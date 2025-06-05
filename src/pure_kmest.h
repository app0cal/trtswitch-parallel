#ifndef PURE_KMEST_H
#define PURE_KMEST_H

#include <vector>
#include <string>

struct KMPoint {
    double time;
    double nrisk;
    double nevent;
    double ncensor;
    double survival;
    double std_err;
    double lower;
    double upper;
};

std::vector<KMPoint> kmest_pure(
    const std::vector<double>& time,
    const std::vector<int>& event,
    const std::string& conftype = "log-log",
    const double conflev = 0.95,
    const bool keep_censor = false
);

#endif