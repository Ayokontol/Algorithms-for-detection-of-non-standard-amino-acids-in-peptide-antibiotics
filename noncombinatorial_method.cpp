#include "noncombinatorial_method.h"

// Noncombinatorial approach with weight
void NoncombinatorialMethod::ncomb_app_w(std::vector<std::vector<double>>& VP, std::vector<std::vector<double> >& VM, std::vector<double>& D) {

    std::set<DataPrepare::m_w> M_W;

    for (long long k = 0; k < (long long)VP.size(); ++k) {
        for (long long i = 0; i < (long long)VP[k].size(); ++i)
            for (long long j = i + 1; j < (long long)VP[k].size(); ++j) {
                double m = fabs(VM[k][i] - VM[k][j]);
                double w = VP[k][i] * VP[k][j];

                if (m < 1 - e || m > MAX_DIFF)
                    continue;

                auto it = M_W.find(DataPrepare::m_w(m, w));
                if (it == M_W.end())
                    M_W.insert(DataPrepare::m_w(m, w));
                else
                    it->w += w;
            }
    }

    std::vector<DataPrepare::m_w> V(M_W.size());
    std::copy(M_W.begin(), M_W.end(), V.begin());

    std::sort(V.begin(), V.end(), [](DataPrepare::m_w a, DataPrepare::m_w b) {return a.w > b.w; });

    std::ofstream out_w("out_m_w.txt");
    for (auto x : M_W)
        out_w << x.m << " " << x.w << std::endl;

    D.resize(D_size);

    for (long long i = 0; i < D_size; ++i)
        D[i] = V[i].m;

    std::ofstream out_w2("out_m_w_stats.txt");
    for (long long i = 0; i < STATS_SIZE; ++i)
        out_w2 << V[i].m << " " << V[i].w << std::endl;
}

// Noncombinatorial approach with weight
void NoncombinatorialMethod::ncomb_app(std::vector<std::vector<double>>& VP, std::vector<std::vector<double>>& VM, std::vector<double>& D) {

    std::set<DataPrepare::m_w> M_W;

    for (long long k = 0; k < (long long)VP.size(); ++k) {
        for (long long i = 0; i < (long long)VP[k].size(); ++i)
            for (long long j = i + 1; j < (long long)VP[k].size(); ++j) {
                double m = fabs(VM[k][i] - VM[k][j]);

                if (m < 1 - e || m > MAX_DIFF)
                    continue;

                auto it = M_W.find(DataPrepare::m_w(m, 1));
                if (it == M_W.end())
                    M_W.insert(DataPrepare::m_w(m, 1));
                else
                    it->w += 1;
            }
    }

    std::vector<DataPrepare::m_w> V(M_W.size());
    std::copy(M_W.begin(), M_W.end(), V.begin());

    std::sort(V.begin(), V.end(), [](DataPrepare::m_w a, DataPrepare::m_w b) {return a.w > b.w; });

    std::ofstream out_w("out_m.txt");
    for (auto x : M_W)
        out_w << x.m << " " << x.w << std::endl;

    D.resize(STATS_SIZE);

    for (long long i = 0; i < D_size; ++i)
        D[i] = V[i].m;

    std::ofstream out_w2("out_m_stats.txt");
    for (long long i = 0; i < STATS_SIZE; ++i)
        out_w2 << V[i].m << " " << V[i].w << std::endl;
}