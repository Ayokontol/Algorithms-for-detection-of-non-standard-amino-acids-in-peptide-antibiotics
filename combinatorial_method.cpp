#include "combinatorial_method.h"

// find E_z (graph for charge z) for one spectrum l
// D -- alphabet
// M_l -- spectrum
// z -- charge state
// E_z -- edges, size = 0
void Combinatorial_method::make_edges_z(std::vector<double>& D, std::vector<double>& M_l, int z, std::vector<std::vector<bool> >& E_z) {
    size_t n = M_l.size();
    E_z.resize(n, std::vector<bool>(n));

    for (size_t i = 0; i < n; ++i)
        for (size_t j = i + 1; j < n; ++j)
            for (auto d : D)
                if (abs(abs(M_l[i] - M_l[j]) - d / z) < e)
                    E_z[i][j] = true;
}

// DFS to find connected components
// E_z -- graph
// C -- colors
// c -- current color
// v -- start point
void Combinatorial_method::DFS(std::vector<std::vector<bool> >& E_z, std::vector<int>& C, int c, size_t v) {
    C[v] = c;
    for (size_t i = 0; i < E_z[v].size(); ++i)
        if (E_z[v][i] == true && C[i] == -1)
            DFS(E_z, C, c, i);
}

size_t Combinatorial_method::make_comp(std::vector<std::vector<bool> >& E_z, std::vector<int>& C) {
    size_t n = E_z.size();

    C.resize(n, -1);

    size_t c = 0;
    for (size_t i = 0; i < n; ++i)
        if (C[i] == -1) {
            DFS(E_z, C, c, i);
            c++;
        }

    return c; // return num of comp
}

// find P_l_z = Pr(D^l | E^l_z)
// E_z -- graph for charge z
// C -- connected components
// P -- vector of intensity
double Combinatorial_method::Pr_l_z(std::vector<std::vector<bool> >& E_z, std::vector<int>& C, std::vector<double>& P_l) {
    size_t n = E_z.size();

    //std::vector<int> C;

    //size_t c = make_comp(E_z, C);
    size_t c = C.size();

    std::vector<double> P_comp(c, 1);

    double P = 0;

    // E_z is a triangular matrix!!!
    for (size_t i = 0; i < n; ++i)
        for (size_t j = i + 1; j < n; ++j)
            if (E_z[i][j] == true && P_comp[C[i]] == 1)
                P_comp[C[i]] *= P_l[i] * P_l[j];

    for (auto p : P_comp)
        P += p;

    return P;
}

// find Pr(D^l | E^l)
// D -- alphabet
// M_l -- spectrum
// P_l -- vector of intensity for spectrum l
// VE_z -- graphs for all charges
// VC -- connected components for all charges
double Combinatorial_method::Pr_l(
    std::vector<double>& D,
    std::vector<double>& M_l,
    std::vector<double>& P_l,
    std::vector<std::vector<std::vector<bool>>>& VE_z,
    std::vector<std::vector<int>>& VC
) {
    double P = 1;

    for (size_t z = 0; z < Z; ++z) {
        P *= Pr_l_z(VE_z[z], VC[z], P_l);
    }

    return P;
}

double Combinatorial_method::Pr_l_log(
    std::vector<double>& D,
    std::vector<double>& M_l,
    std::vector<double>& P_l,
    std::vector<std::vector<std::vector<bool>>>& VE_z, 
    std::vector<std::vector<int>>& VC
) {
    double P_log = 0;

    for (size_t z = 0; z < Z; ++z) {
        P_log += std::log(Pr_l_z(VE_z[z], VC[z], P_l));
    }

    return P_log;
}

// Likelihood = Ï_l Ï_z Sum_comp Ï_edg p_i * p_j
// D -- alphabet
// M -- vector of vector of m/z
// P -- vector of vector of intensity
double Combinatorial_method::likelihood(std::vector<double>& D, std::vector<std::vector<double> >& M, std::vector<std::vector<double> >& P) {
    size_t n = M.size();
    double Lklh = 1;

    for (size_t i = 0; i < n; ++i)
        Lklh *= Pr_l(D, M[i], P[i], Curr_VE_z[i], Curr_VC[i]);

    return Lklh;
}

// Likelihood = Ï_l Ï_z Sum_comp Ï_edg p_i * p_j
// D -- alphabet
// M -- vector of vector of m/z
// P -- vector of vector of intensity
double Combinatorial_method::likelihood_log(std::vector<double>& D, std::vector<std::vector<double> >& M, std::vector<std::vector<double> >& P) {
    size_t n = M.size();
    double Lklh = 1;

    for (size_t i = 0; i < n; ++i)
        Lklh += Pr_l_log(D, M[i], P[i], Curr_VE_z[i], Curr_VC[i]);

    return Lklh;
}


// check P(D | D[i] = new_d) == 1
// D -- alphabet
// index -- index of D[i]
// new_d -- D[i] = new_d
bool Combinatorial_method::check(std::vector<double>& D, size_t index, double new_d) {
    std::vector<double> Scale = { 3, 3 / 2, 2 / 3, 2, 1 / 3, 1 / 2 };

    if (new_d < 1 - e || new_d > MAX_DIFF)
        return false;

    for (size_t i = 0; i < D.size(); ++i) {
        if (i == index)
            continue;

        if (abs(D[i] - new_d) < 0.5)
            return false;

        for (size_t j = 0; j < Scale.size(); ++j)
            if (abs(1 - D[i] / new_d * Scale[j]) < e)
                return false;
    }
    return true;
}

// sampling with scales
// D -- alphabet
// VM -- vector of vector of m/z
// VP -- vector of vector of intensity
// i -- index of D[i]
bool Combinatorial_method::sampling_2(std::vector<double>& D, std::vector<std::vector<double> >& VM, std::vector<std::vector<double> >& VP, size_t i) {
    std::srand(unsigned(std::time(0)));

    std::vector<double> Scale = { 3, 3 / 2, 2 / 3, 2, 1 / 3, 1 / 2 };

    double d = D[i] * Scale[rand() % Scale.size()];
    double d_old = D[i];

    size_t count = 0;
    while (!check(D, i, d)) {
        d = D[i] * Scale[rand() % Scale.size()];
        if (count > 20)
            return false;
        count++;
    }

    double lklh_old = Curr_lklh;
    auto Old_VE_z = std::move(Curr_VE_z);
    auto Old_VC = std::move(Curr_VC);

    Curr_VE_z.resize(VM.size(), std::vector<std::vector<std::vector<bool>>>(3));
    Curr_VC.resize(VM.size(), std::vector<std::vector<int>>(3));

    D[i] = d;

    for (size_t l = 0; l < VM.size(); ++l)
        for (int z = 0; z < Z; ++z)
            make_edges_z(D, VM[l], z + 1, Curr_VE_z[l][z]);
    for (size_t l = 0; l < VM.size(); ++l)
        for (int z = 0; z < Z; ++z)
            make_comp(Curr_VE_z[l][z], Curr_VC[l][z]);

    Curr_lklh = likelihood_log(D, VM, VP);

    if (lklh_old >= Curr_lklh) {
        D[i] = d_old;
        Curr_VE_z = std::move(Old_VE_z);
        Curr_VC = std::move(Old_VC);
        Curr_lklh = lklh_old;
        return false;
    }
    return true;
}

//sampling
// D -- alphabet
// VM -- vector of vector of m/z
// VP -- vector of vector of intensity
// Diff -- vector of differences
// i -- index of D[i]
bool Combinatorial_method::sampling_1(std::vector<double>& D, std::vector<std::vector<double> >& VM, std::vector<std::vector<double> >& VP, std::vector<double>& Diff, size_t i) {
    std::srand(unsigned(std::time(0)));

    double d = Diff[rand() % Diff.size()];
    double d_old = D[i];

    while (!check(D, i, d)) {
        d = Diff[rand() % Diff.size()];
    }

    double lklh_old = Curr_lklh;
    auto Old_VE_z = std::move(Curr_VE_z);
    auto Old_VC = std::move(Curr_VC);

    Curr_VE_z.resize(VM.size(), std::vector<std::vector<std::vector<bool>>>(3));
    Curr_VC.resize(VM.size(), std::vector<std::vector<int>>(3));

    D[i] = d;

    for (size_t l = 0; l < VM.size(); ++l)
        for (int z = 0; z < Z; ++z)
            make_edges_z(D, VM[l], z + 1, Curr_VE_z[l][z]);

    //check if the same
    //bool fl_the_same = true;
    //for (size_t l = 0; l < VM.size(); ++l)
    //    for (int z = 0; z < Z; ++z)
    //        for (size_t i = 0; i < Curr_VE_z[l][z].size(); ++i)
    //            for (size_t j = 0; j < Curr_VE_z[l][z][i].size(); ++j) {
    //                //if (Curr_VE_z[l][z][i][j] == true || Old_VE_z[l][z][i][j] == true)
    //                //    std::cout << "not empty\n";
    //                if (Curr_VE_z[l][z][i][j] != Old_VE_z[l][z][i][j]) {
    //                    fl_the_same = false;
    //                    std::cout << "not the same\n";
    //                    break;
    //                }
    //            }
    //if (fl_the_same)
    //    std::cout << "the same\n";

    for (size_t l = 0; l < VM.size(); ++l)
        for (int z = 0; z < Z; ++z)
            make_comp(Curr_VE_z[l][z], Curr_VC[l][z]);

    Curr_lklh = likelihood_log(D, VM, VP);

    if (lklh_old >= Curr_lklh) {
        //if (lklh_old == Curr_lklh)
        //    std::cout << "the same lklhd\n";
        D[i] = d_old;
        Curr_VE_z = std::move(Old_VE_z);
        Curr_VC = std::move(Old_VC);
        Curr_lklh = lklh_old;
        return false;
    }
    return true;
}

void Combinatorial_method::prepare_to_sampling(
    std::vector<double>& D,
    std::vector<std::vector<double> >& VM,
    std::vector<std::vector<double> >& VP,
    std::vector<double>& Diff,
    std::vector<size_t>& Last_change,
    size_t old_samples) {

    Curr_VE_z.resize(VM.size(), std::vector<std::vector<std::vector<bool>>>(3));
    Curr_VC.resize(VM.size(), std::vector<std::vector<int>>(3));
    for (size_t l = 0; l < VM.size(); ++l)
        for (int z = 0; z < Z; ++z)
            make_edges_z(D, VM[l], z + 1, Curr_VE_z[l][z]);
    for (size_t l = 0; l < VM.size(); ++l)
        for (int z = 0; z < Z; ++z)
            make_comp(Curr_VE_z[l][z], Curr_VC[l][z]);

    Curr_lklh = likelihood_log(D, VM, VP);
}

// D -- alphabet (random or from file)
// VM -- all spectra
// VP -- all intensity
void Combinatorial_method::sampling(
    std::vector<double>& D,
    std::vector<std::vector<double> >& VM,
    std::vector<std::vector<double> >& VP,
    std::vector<double>& Diff,
    std::vector<size_t>& Last_change,
    size_t old_samples,
    size_t i) {

    bool tmp;
    size_t j = rand() % D_size;
    tmp = sampling_1(D, VM, VP, Diff, j);
    if (tmp) {
        Last_change[j] = i + old_samples;
        //std::cout << "1st change ";
    }
    std::cout << i << "\n";
    j = rand() % D_size;
    tmp = sampling_2(D, VM, VP, j);
    if (tmp) {
        Last_change[j] = i + old_samples;
        //std::cout << "2nd change\n";
    }
}
