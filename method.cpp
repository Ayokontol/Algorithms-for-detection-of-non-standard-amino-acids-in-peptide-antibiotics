#include "method.h"

void Method::start_sampling(std::vector<double>& D, std::vector<std::vector<double> >& VM, std::vector<std::vector<double> >& VP) {
    std::srand(unsigned(std::time(0)));

    D.resize(D_size, -1);

    std::vector<double> Diff;

    if (!FROM_FILE) {
        dataPrepare.make_differences(VM, Diff);

        to_out_file_diff(Diff);
    }
    else {
        std::ifstream ifs_diff(name_diff_file);
        size_t num_of_diff;

        ifs_diff >> num_of_diff;

        Diff.resize(num_of_diff);

        for (size_t i = 0; i < num_of_diff; ++i)
            ifs_diff >> Diff[i];
    }

    if (FROM_FILE)
        dataPrepare.get_D_from_out_file(D);
    else {
        //random start
        /*for (size_t i = 0; i < D_size; ++i) {
            D[i] = Diff[rand() % Diff.size()];
        }*/

        //start with noncominatorial algo
        NoncombinatorialMethod noncombMethod;

        noncombMethod.ncomb_app_w(VP, VM, D);
    }

    if (!FROM_FILE)
        dataPrepare.Last_change.resize(D_size, 0);  ///!

    combMethod.prepare_to_sampling(D, VM, VP, Diff, dataPrepare.Last_change, dataPrepare.old_samples);

    for (size_t i = 0; i < Num_of_sampl; ++i) {
        combMethod.sampling(D, VM, VP, Diff, dataPrepare.Last_change, dataPrepare.old_samples, i);
        if (i % 100 == 0 && i > 0) {
            to_out_file(D, i);
        }
    }
}

void Method::to_out_file(std::vector<double>& D, size_t num_of_current_samples) {
    // read amino acids
    std::ifstream amino_acids("amino_acids.txt");

    std::map<std::string, double> Amino_acids_masses;
    size_t num_of_amio_acids;

    amino_acids >> num_of_amio_acids;

    for (size_t i = 0; i < num_of_amio_acids; ++i) {
        std::string name;
        double mass;
        amino_acids >> name >> mass;
        Amino_acids_masses[name] = mass;
    }

    //out
    std::ofstream ofs(name_out_file);
    ofs << dataPrepare.old_samples + num_of_current_samples << std::endl;
    ofs << "Mass, last change, name\n";
    for (size_t i = 0; i < D_size; ++i) {
        ofs << D[i] << " " << dataPrepare.Last_change[i] << " ";

        ofs << "-";
        for (auto x : Amino_acids_masses)
            if (abs(x.second - D[i]) < e) {
                ofs << x.first;
                break;
            }
        ofs << std::endl;
    }
}

void Method::to_out_file_diff(std::vector<double>& Diff) {
    std::ofstream ofs(name_diff_file);

    ofs << Diff.size() << "\n";
    for (auto x : Diff)
        ofs << x << "\n";
}

