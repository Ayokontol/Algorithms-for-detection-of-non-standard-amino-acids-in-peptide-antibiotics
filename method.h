#ifndef METHOD_H
#define METHOD_H

#include <thread>
#include <algorithm>

#include "combinatorial_method.h"
#include "data_prepare.h"
#include "noncombinatorial_method.h"

class Method {
    DataPrepare dataPrepare;
    Combinatorial_method combMethod;
    std::vector<std::thread> Pool;

    void to_out_file_diff(std::vector<double>& Diff);
public:
    Method();
    void start_sampling(std::vector<double>& D, std::vector<std::vector<double> >& VM, std::vector<std::vector<double> >& VP);
    void start_sampling_thread(std::vector<double>& D, std::vector<std::vector<double> >& VM, std::vector<std::vector<double> >& VP);
    void to_out_file(std::vector<double>& D, long long num_of_current_samples = Num_of_sampl);
};

#endif