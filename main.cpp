#include "method.h"
#include <string>

int main() {
    //std::ifstream ifs("data_glyco.txt"); //2
    std::ifstream ifs("data3.txt");      //1   //peptide data
    //std::ifstream ifs("test.txt");      //2

    long long k = 0;
    std::vector<std::vector<double> > VM;
    std::vector<std::vector<double> > VP;

    int type;
    ifs >> type;

    while (!ifs.eof()) {
        double p, m;

        if (type == 1) {
            long long num;
            long long n;

            ifs >> num >> n;

            if (ifs.eof())
                break;

            std::vector<double> P; // intensity array
            std::vector<double> M; // m/z array

            double sum_p = 0;
            for (long long i = 0; i < n; ++i) {
                ifs >> m >> p;
                P.push_back(p);
                M.push_back(m);
                sum_p += p;
            }

            std::for_each(P.begin(), P.end(), [sum_p](double& x) { x = x / sum_p * 1e5; });

            VM.push_back(M);
            VP.push_back(P);
        }

        //if (ifs.eof())
        //    break;

        if (type == 2) {
            std::vector<double> P; // intensity array
            std::vector<double> M; // m/z array

            double sum_p = 0;
            while (1) {
                int t;
                ifs >> m;
                if (m == -1)
                    break;
                ifs >> p >> t;
                P.push_back(p);
                M.push_back(m);
                sum_p += p;
            }

            std::for_each(P.begin(), P.end(), [sum_p](double& x) { x = x / sum_p; });

            VM.push_back(M);
            VP.push_back(P);
        }

        /*std::vector<double> D1, D2;
        ncomb_app_w(P, M, D1);
        ncomb_app(P, M, D2);

        ofs << num << ") ";
        for (auto x : D1)
            ofs << x << " ";
        ofs << "\n";
        for (auto x : D2)
            ofs << x << " ";
        ofs << "\n\n";*/
    }

    std::vector<double> D;
    //sampling(D, VM, VP);
    std::vector<double> D1, D2;
    //ncomb_app_w(VP, VM, D1);
    //ncomb_app(VP, VM, D2);

    Method method;

    if (IS_PARALLEL)
        method.start_sampling_thread(D, VM, VP);
    else
        method.start_sampling(D, VM, VP);

    /*for (auto x : D)
        ofs << x << " ";
    ofs << std::endl;
    for (auto x : Last_change)
        ofs << x << " ";*/
    method.to_out_file(D);

    return 0;
}

