#include "method.h"

Method::Method() {
    //Pool.resize(NUM_OF_PROC);
}

void Method::start_sampling(std::vector<double>& D, std::vector<std::vector<double> >& VM, std::vector<std::vector<double> >& VP) {
    std::srand(unsigned(std::time(0)));

    std::cout << "START READING\n";

    D.resize(D_size, -1);

    std::vector<double> Diff;

    if (!FROM_FILE) {
        dataPrepare.make_differences(VM, Diff);

        to_out_file_diff(Diff);
    }
    else {
        std::ifstream ifs_diff(name_diff_file);
        long long num_of_diff;

        ifs_diff >> num_of_diff;

        Diff.resize(num_of_diff);

        for (long long i = 0; i < num_of_diff; ++i)
            ifs_diff >> Diff[i];
    }

    if (FROM_FILE)
        dataPrepare.get_D_from_out_file(D);
    else {
        //random start
        /*for (long long i = 0; i < D_size; ++i) {
            D[i] = Diff[rand() % Diff.size()];
        }*/

        //start with noncominatorial algo
        NoncombinatorialMethod noncombMethod;

        noncombMethod.ncomb_app_w(VP, VM, D);
    }

    if (!FROM_FILE)
        dataPrepare.Last_change.resize(D_size, 0);  ///!

    combMethod.prepare_to_sampling(D, VM, VP, Diff, dataPrepare.Last_change, dataPrepare.old_samples);
    std::cout << "START SAMPLING\n";
    for (long long i = 0; i < Num_of_sampl; ++i) {
        combMethod.sampling(D, VM, VP, Diff, dataPrepare.Last_change, dataPrepare.old_samples, i);
        if (i % 10 == 0 && i > 0) {
            to_out_file(D, i);
        }
    }
}

void Method::start_sampling_thread(std::vector<double>& D, std::vector<std::vector<double> >& VM, std::vector<std::vector<double> >& VP) {
    std::srand(unsigned(std::time(0)));
    D.resize(D_size, -1);

    std::vector<double> Diff;

    if (!FROM_FILE) {
        dataPrepare.make_differences(VM, Diff);

        to_out_file_diff(Diff);
    }
    else {
        std::ifstream ifs_diff(name_diff_file);
        long long num_of_diff;

        ifs_diff >> num_of_diff;

        Diff.resize(num_of_diff);

        for (long long i = 0; i < num_of_diff; ++i)
            ifs_diff >> Diff[i];
    }

    if (FROM_FILE)
        dataPrepare.get_D_from_out_file(D);
    else {
        //random start
        /*for (long long i = 0; i < D_size; ++i) {
            D[i] = Diff[rand() % Diff.size()];
        }*/

        //start with noncominatorial algo
        NoncombinatorialMethod noncombMethod;

        noncombMethod.ncomb_app_w(VP, VM, D);
    }

    if (!FROM_FILE)
        dataPrepare.Last_change.resize(D_size, 0);  ///!

    combMethod.prepare_to_sampling(D, VM, VP, Diff, dataPrepare.Last_change, dataPrepare.old_samples);
    std::cout << "START SAMPLING\n";

    for (long long i = 0; i < Num_of_sampl; ++i) {
        int j = rand() % D_size;
        std::vector<double> All_new_d(NUM_OF_PROC);
        std::vector<double> All_new_lklh(NUM_OF_PROC);
        std::vector< std::vector<std::vector<std::vector<std::vector<bool>>>> > All_VE_z(NUM_OF_PROC);
        std::vector< std::vector<std::vector<std::vector<long long>>> > All_VC(NUM_OF_PROC);
        for (int ii = 0; ii < NUM_OF_PROC; ++ii) {
            Pool.push_back(std::thread(
                &Combinatorial_method::sampling_for_thread,
                &combMethod,
                std::ref(D),
                std::ref(VM),
                std::ref(VP),
                std::ref(Diff),
                std::ref(dataPrepare.Last_change), //remove
                dataPrepare.old_samples,           //remove
                i,                                 //remove
                ii,
                j,
                std::ref(All_new_d[ii]),
                std::ref(All_new_lklh[ii]),
                std::ref(All_VE_z[ii]),
                std::ref(All_VC[ii])));
        }

        for (int ii = 0; ii < NUM_OF_PROC; ++ii)
            if (Pool[Pool.size() - NUM_OF_PROC + ii].joinable())
                Pool[Pool.size() - NUM_OF_PROC + ii].join();

/*        for (auto x : All_new_lklh)
            std::cout << x << " ";
        std::cout << std::endl;
        for (auto x : All_new_d)
            std::cout << x << " ";
        std::cout << std::endl*/;
        int max_arg = std::distance(All_new_lklh.begin(), std::max_element(All_new_lklh.begin(), All_new_lklh.end()));

        if (All_new_lklh[max_arg] > combMethod.Curr_lklh) {
            D[j] = All_new_d[max_arg];
            combMethod.Curr_lklh = All_new_lklh[max_arg];
            combMethod.Curr_VC = std::move(All_VC[max_arg]);
            combMethod.Curr_VE_z = std::move(All_VE_z[max_arg]);

            dataPrepare.Last_change[j] = i * NUM_OF_PROC + dataPrepare.old_samples;
        }


        std::cout << i << "/" << Num_of_sampl << ", iter * NUM_OF_PROC: " << i * NUM_OF_PROC << ", lklh: " << combMethod.Curr_lklh << "\n";
        if (i % 10 == 0 && i > 0) {
            to_out_file(D, i);
        }
    }
}


void Method::to_out_file(std::vector<double>& D, long long num_of_current_samples) {
    // read amino acids
    std::ifstream amino_acids("amino_acids.txt");

    std::map<std::string, double> Amino_acids_masses;
    long long num_of_amio_acids;

    amino_acids >> num_of_amio_acids;

    for (long long i = 0; i < num_of_amio_acids; ++i) {
        std::string name;
        double mass;
        amino_acids >> name >> mass;
        Amino_acids_masses[name] = mass;
    }

    //out
    std::ofstream ofs(name_out_file);
    if (IS_PARALLEL)
        ofs << dataPrepare.old_samples + num_of_current_samples * NUM_OF_PROC << std::endl;
    else
        ofs << dataPrepare.old_samples + num_of_current_samples << std::endl;
    ofs << "Mass, last change, name\n";
    for (long long i = 0; i < D_size; ++i) {
        ofs << D[i] << " " << dataPrepare.Last_change[i] << " ";

        ofs << "-";
        for (auto x : Amino_acids_masses)
            if (fabs(x.second - D[i]) < e) {
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

