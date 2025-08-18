#include "pure_kmest.h"
#include <cmath>
#include <stdexcept>
#include <algorithm>

static double norminv(double p) {
    if (std::abs(p - 0.975) < 1e-6)
        return 1.95996398454005;
    throw std::runtime_error("norminv() only supports 95% CI for now");
}

// Compute CI based on conftype
static std::pair<double, double> compute_ci(double surv, double se, const std::string& type, double z) {
  if (type == "log-log") {
      if (surv <= 0.0 || surv >= 1.0) return {0.0, 1.0};
      double loglog = std::log(-std::log(surv));
      double se_trans = se / (surv * std::log(surv));
      double lower = std::exp(-std::exp(loglog + z * se_trans));
      double upper = std::exp(-std::exp(loglog - z * se_trans));
      return {lower, upper};
  }
  return {0.0, 1.0};  // fallback if unsupported type
}

std::vector<KMPoint> kmest_pure(
  const std::vector<double>& time,
  const std::vector<int>& event,
  const std::string& conftype,
  const double conflev,
  const bool keep_censor
) {
  if (time.size() != event.size())
      throw std::runtime_error("time and event must be the same length");

  int n = static_cast<int>(time.size());
  std::vector<int> indices(n);
  for (int i = 0; i < n; ++i) indices[i] = i;

  // Sort by time ascending, then event descending
  std::sort(indices.begin(), indices.end(), [&](int i, int j) {
      return (time[i] < time[j]) || (time[i] == time[j] && event[i] > event[j]);
  });

  std::vector<double> t_sorted(n);
  std::vector<int> e_sorted(n);
  for (int i = 0; i < n; ++i) {
      t_sorted[i] = time[indices[i]];
      e_sorted[i] = event[indices[i]];
  }

  std::vector<KMPoint> output;
  double surv = 1.0;
  double vcumhaz = 0.0;
  double z = norminv((1.0 + conflev) / 2.0);

  double current_time = -1;
  double nrisk = n, nevent = 0, ncensor = 0;
  bool cache = false;

  for (int i = 0; i < n; ++i) {
      if (t_sorted[i] != current_time) {
          if (i > 0 && cache) {
              surv *= (1.0 - nevent / nrisk);
              if (nrisk > nevent)
                  vcumhaz += nevent / (nrisk * (nrisk - nevent));
              double se = surv * std::sqrt(vcumhaz);
              auto [lower, upper] = compute_ci(surv, se, conftype, z);

              output.push_back({current_time, nrisk, nevent, ncensor, surv, se, lower, upper});
          }

          current_time = t_sorted[i];
          nrisk = n - i;
          nevent = 0;
          ncensor = 0;
          cache = false;
      }

      if (e_sorted[i] == 1) {
          nevent++;
          cache = true;
      } else {
          ncensor++;
          if (keep_censor) cache = true;
      }
  }

  // Final event
  if (cache) {
      surv *= (1.0 - nevent / nrisk);
      if (nrisk > nevent)
          vcumhaz += nevent / (nrisk * (nrisk - nevent));
      double se = surv * std::sqrt(vcumhaz);
      auto [lower, upper] = compute_ci(surv, se, conftype, z);
      output.push_back({current_time, nrisk, nevent, ncensor, surv, se, lower, upper});
  }

  return output;
}