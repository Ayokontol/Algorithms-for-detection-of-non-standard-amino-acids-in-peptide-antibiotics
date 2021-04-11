#pragma once

#include "data_prepare.h"

#include <vector>
#include <set>
#include <ctime>
#include <cmath>
#include <algorithm>
#include <map>
#include <iostream>

class NoncombinatorialMethod {
public:
    // Noncombinatorial approach with weight
    void ncomb_app_w(std::vector<std::vector<double>>& VP, std::vector<std::vector<double> >& VM, std::vector<double>& D);

    // Noncombinatorial approach with weight
    void ncomb_app(std::vector<std::vector<double>>& VP, std::vector<std::vector<double>>& VM, std::vector<double>& D);
};