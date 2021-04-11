#include "data_prepare.h"

// VM -- vector of vector of m/z for all spectra
// Diff -- result vector of differences
void DataPrepare::make_differences(std::vector<std::vector<double> >& VM, std::vector<double>& Diff) {
    Diff.resize(0);
    std::set<m_w> S;

    for (auto M : VM)
        for (size_t i = 0; i < M.size(); ++i)
            for (size_t j = i + 1; j < M.size(); ++j) {
                double m = abs(M[i] - M[j]);

                if (m < 1 - e || m > MAX_DIFF)
                    continue;

                auto it = S.find(m_w(m, 0));
                if (it == S.end()) {
                    S.insert(m_w(m, 0));
                    Diff.push_back(m);
                }
            }
}

void DataPrepare::get_D_from_out_file(std::vector<double>& D) {
    std::ifstream ifs(name_out_file);
    D.resize(D_size);
    Last_change.resize(D_size);
    std::string name;
    double last_change, d;
    size_t delta;
    ifs >> delta;
    old_samples += delta;
    //Mass, last change, name
    ifs >> name;
    ifs >> name;
    ifs >> name;
    ifs >> name;
    for (size_t i = 0; i < D_size; ++i) {
        ifs >> d >> last_change >> name;
        D[i] = d;
        Last_change[i] = last_change;
    }
}

