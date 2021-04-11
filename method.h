#pragma once

#include "combinatorial_method.h"
#include "data_prepare.h"
#include "noncombinatorial_method.h"

class Method {
    DataPrepare dataPrepare;
    Combinatorial_method combMethod;

    void to_out_file_diff(std::vector<double>& Diff);
public:
    void start_sampling(std::vector<double>& D, std::vector<std::vector<double> >& VM, std::vector<std::vector<double> >& VP);
    void to_out_file(std::vector<double>& D, size_t num_of_current_samples = Num_of_sampl);
};