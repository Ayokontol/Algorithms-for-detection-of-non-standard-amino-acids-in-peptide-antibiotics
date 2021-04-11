#pragma once

#include "combinatorial_method.h"
#include "const.h"

#include <fstream>
#include <iostream>

class DataPrepare {
public:
    size_t old_samples = 0;

    std::vector<size_t> Last_change;

    struct m_w {
        double m;           // m/z
        mutable double w;   // weight

        m_w(double m, double w) : m(m), w(w) {}
        m_w() {}

        bool operator<(const m_w& other) const {
            if (abs(m - other.m) < e)
                return false;
            return m < other.m;
        }
    };

    void get_D_from_out_file(std::vector<double>& D);

    // VM -- vector of vector of m/z for all spectra
    // Diff -- result vector of differences
    void make_differences(std::vector<std::vector<double> >& VM, std::vector<double>& Diff);
};
