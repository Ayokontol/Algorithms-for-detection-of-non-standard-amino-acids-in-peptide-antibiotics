#include "combinatorial_method.h"

// find E_z (graph for charge z) for one spectrum l
// D -- alphabet
// M_l -- spectrum
// z -- charge state
// E_z -- edges, size = 0
void Combinatorial_method::make_edges_z(
    std::vector<double> const & D,
    std::vector<double> const & M_l,
    int z,
    std::vector<std::vector<bool> >& E_z) const {
    long long n = M_l.size();
    E_z.resize(n, std::vector<bool>(n));

    for (long long i = 0; i < n; ++i)
        for (long long j = i + 1; j < n; ++j)
            for (auto d : D)
                if (fabs(fabs(M_l[i] - M_l[j]) - d / z) < e)
                    E_z[i][j] = true;
}

// DFS to find connected components
// E_z -- graph
// C -- colors
// c -- current color
// v -- start point
void Combinatorial_method::DFS(std::vector<std::vector<bool> > const & E_z, std::vector<long long>& C, long long c, long long v) const {
    C[v] = c;
    for (long long i = 0; i < (long long)E_z[v].size(); ++i)
        if (E_z[v][i] == true && C[i] == -1)
            DFS(E_z, C, c, i);
}

long long Combinatorial_method::make_comp(std::vector<std::vector<bool> > const & E_z, std::vector<long long>& C) const {
    long long n = E_z.size();

    C.resize(n, -1);

    long long c = 0;
    for (long long i = 0; i < n; ++i)
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
double Combinatorial_method::Pr_l_z(
    std::vector<std::vector<bool>> const & E_z,
    std::vector<long long> const & C,
    std::vector<double> const & P_l
) const {
    long long n = E_z.size();

    //std::vector<int> C;

    //long long c = make_comp(E_z, C);
    long long c = C.size();

    std::vector<double> P_comp(c, 1);

    double P = 0;

    // E_z is a triangular matrix!!!
    for (long long i = 0; i < n; ++i)
        for (long long j = i + 1; j < n; ++j)
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
    std::vector<double> const & D,
    std::vector<double> const & M_l,
    std::vector<double> const & P_l,
    std::vector<std::vector<std::vector<bool>>> const & VE_z,
    std::vector<std::vector<long long>> const & VC
) const {
    double P = 1;

    for (long long z = 0; z < Z; ++z) {
        P *= Pr_l_z(VE_z[z], VC[z], P_l);
    }

    return P;
}

double Combinatorial_method::Pr_l_log(
    std::vector<double> const & D,
    std::vector<double> const & M_l,
    std::vector<double> const & P_l,
    std::vector<std::vector<std::vector<bool>>> const & VE_z,
    std::vector<std::vector<long long>> const & VC
) const {
    double P_log = 0;

    for (long long z = 0; z < Z; ++z) {
        P_log += std::log(Pr_l_z(VE_z[z], VC[z], P_l));
    }

    return P_log;
}

// Likelihood = Ï_l Ï_z Sum_comp Ï_edg p_i * p_j
// D -- alphabet
// M -- vector of vector of m/z
// P -- vector of vector of intensity
double Combinatorial_method::likelihood(
    std::vector<double>& D, 
    std::vector<std::vector<double> >& M, 
    std::vector<std::vector<double> >& P,
    std::vector<std::vector<std::vector<std::vector<bool>>>>& VE_z,
    std::vector<std::vector<std::vector<long long>>>& VC
) const {
    long long n = M.size();
    double Lklh = 1;

    for (long long i = 0; i < n; ++i)
        Lklh *= Pr_l(D, M[i], P[i], VE_z[i], VC[i]);

    return Lklh;
}

// Likelihood = Ï_l Ï_z Sum_comp Ï_edg p_i * p_j
// D -- alphabet
// M -- vector of vector of m/z
// P -- vector of vector of intensity
double Combinatorial_method::likelihood_log(
    std::vector<double> const & D,
    std::vector<std::vector<double> > const & M,
    std::vector<std::vector<double> > const & P,
    std::vector<std::vector<std::vector<std::vector<bool>>>> const & VE_z,
    std::vector<std::vector<std::vector<long long>>> const & VC
) const {
    long long n = M.size();
    double Lklh = 1;

    for (long long i = 0; i < n; ++i)
        Lklh += Pr_l_log(D, M[i], P[i], VE_z[i], VC[i]);

    return Lklh;
}

// check P(D | D[i] = new_d) == 1
// D -- alphabet
// index -- index of D[i]
// new_d -- D[i] = new_d
bool Combinatorial_method::check(std::vector<double> const & D, long long index, double new_d) const {
    std::vector<double> Scale = { 3, 3 / 2, 2 / 3, 2, 1 / 3, 1 / 2 };

    if (new_d < 1.0 - e || new_d > MAX_DIFF)
        return false;

    for (long long i = 0; i < (long long)D.size(); ++i) {
        if (i == index)
            continue;

        if (fabs(D[i] - new_d) < 0.5)
            return false;
        for (long long j = 0; j < (long long)Scale.size(); ++j)
            if (fabs(1 - D[i] / new_d * Scale[j]) < e) {
                return false;
            }
    }
    return true;
}

// sampling with scales
// D -- alphabet
// VM -- vector of vector of m/z
// VP -- vector of vector of intensity
// i -- index of D[i]
bool Combinatorial_method::sampling_2(
    std::vector<double>& D,
    std::vector<std::vector<double> > const & VM,
    std::vector<std::vector<double> > const & VP,
    long long i
) {
    std::srand(unsigned(std::time(0)));

    std::vector<double> Scale = { 3, 3 / 2, 2 / 3, 2, 1 / 3, 1 / 2 };

    double d = D[i] * Scale[rand() % Scale.size()];
    double d_old = D[i];

    long long count = 0;
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
    Curr_VC.resize(VM.size(), std::vector<std::vector<long long>>(3));

    D[i] = d;

    for (long long l = 0; l < (long long)VM.size(); ++l)
        for (int z = 0; z < Z; ++z)
            make_edges_z(D, VM[l], z + 1, Curr_VE_z[l][z]);
    for (long long l = 0; l < (long long)VM.size(); ++l)
        for (int z = 0; z < Z; ++z)
            make_comp(Curr_VE_z[l][z], Curr_VC[l][z]);

    Curr_lklh = likelihood_log(D, VM, VP, Curr_VE_z, Curr_VC);

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
bool Combinatorial_method::sampling_1(
    std::vector<double>& D,
    std::vector<std::vector<double> > const & VM,
    std::vector<std::vector<double> > const & VP,
    std::vector<double> const & Diff,
    long long i
) {
    std::srand((std::time(NULL)));

    double d = Diff[rand() % Diff.size()];
    double d_old = D[i];

    long long count = 0;
    while (!check(D, i, d)) {
        int t = rand() % Diff.size();
        //std::cout << Diff[t] << "\n";
        d = Diff[t];
        if (count > 20)
            return false;
        count++;
    }

    double lklh_old = Curr_lklh;
    auto Old_VE_z = std::move(Curr_VE_z);
    auto Old_VC = std::move(Curr_VC);

    Curr_VE_z.resize(VM.size(), std::vector<std::vector<std::vector<bool>>>(3));
    Curr_VC.resize(VM.size(), std::vector<std::vector<long long>>(3));

    D[i] = d;

    for (long long l = 0; l < (long long)VM.size(); ++l)
        for (int z = 0; z < Z; ++z)
            make_edges_z(D, VM[l], z + 1, Curr_VE_z[l][z]);

    for (long long l = 0; l < (long long)VM.size(); ++l)
        for (int z = 0; z < Z; ++z)
            make_comp(Curr_VE_z[l][z], Curr_VC[l][z]);

    Curr_lklh = likelihood_log(D, VM, VP, Curr_VE_z, Curr_VC);

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

//sampling
// D -- alphabet
// VM -- vector of vector of m/z
// VP -- vector of vector of intensity
// Diff -- vector of differences
// i -- index of D[i]
bool Combinatorial_method::sampling_1_for_thread(
    std::vector<double> const & D, 
    std::vector<std::vector<double> > const & VM, 
    std::vector<std::vector<double> > const & VP, 
    std::vector<double> const & Diff, 
    long long i,
    int cur_proc,
    double& new_d,
    double& new_lklh,
    std::vector<std::vector<std::vector<std::vector<bool>>>>& New_VE_z,
    std::vector<std::vector<std::vector<long long>>>& New_VC
) const {
    std::srand(std::time(NULL) + cur_proc);

    new_lklh = -INF;

    new_d = Diff[rand() % Diff.size()];
    std::vector<double> New_D(D);

    long long count = 0;
    while (!check(D, i, new_d)) {
        int t = rand() % Diff.size();
        //std::cout << Diff[t] << "\n";
        new_d = Diff[t];
        if (count > 20)
            return false;
        count++;
    }

    //std::cout << new_d << "\n";

    New_VE_z.resize(VM.size(), std::vector<std::vector<std::vector<bool>>>(3));
    New_VC.resize(VM.size(), std::vector<std::vector<long long>>(3));

    New_D[i] = new_d;

    for (long long l = 0; l < (long long)VM.size(); ++l)
        for (int z = 0; z < Z; ++z)
            make_edges_z(New_D, VM[l], z + 1, New_VE_z[l][z]);

    for (long long l = 0; l < (long long)VM.size(); ++l)
        for (int z = 0; z < Z; ++z)
            make_comp(New_VE_z[l][z], New_VC[l][z]);

    new_lklh = likelihood_log(New_D, VM, VP, New_VE_z, New_VC);
    return true;
}

// sampling with scales
// D -- alphabet
// VM -- vector of vector of m/z
// VP -- vector of vector of intensity
// i -- index of D[i]
bool Combinatorial_method::sampling_2_for_thread(
    std::vector<double>& D,
    std::vector<std::vector<double> > const& VM,
    std::vector<std::vector<double> > const& VP,
    long long i,
    int cur_proc,
    double& new_d,
    double& new_lklh,
    std::vector<std::vector<std::vector<std::vector<bool>>>>& New_VE_z,
    std::vector<std::vector<std::vector<long long>>>& New_VC
) const {
    std::srand(unsigned(std::time(0)));

    new_lklh = -INF;

    std::vector<double> Scale = { 3, 3 / 2, 2 / 3, 2, 1 / 3, 1 / 2 };

    new_d = D[i] * Scale[rand() % Scale.size()];
    std::vector<double> New_D(D);

    long long count = 0;
    while (!check(D, i, new_d)) {
        new_d = D[i] * Scale[rand() % Scale.size()];
        if (count > 20)
            return false;
        count++;
    }

    New_VE_z.resize(VM.size(), std::vector<std::vector<std::vector<bool>>>(3));
    New_VC.resize(VM.size(), std::vector<std::vector<long long>>(3));

    New_D[i] = new_d;

    for (long long l = 0; l < (long long)VM.size(); ++l)
        for (int z = 0; z < Z; ++z)
            make_edges_z(New_D, VM[l], z + 1, New_VE_z[l][z]);

    for (long long l = 0; l < (long long)VM.size(); ++l)
        for (int z = 0; z < Z; ++z)
            make_comp(New_VE_z[l][z], New_VC[l][z]);

    new_lklh = likelihood_log(New_D, VM, VP, New_VE_z, New_VC);
    return true;
}


void Combinatorial_method::prepare_to_sampling(
    std::vector<double> const & D,
    std::vector<std::vector<double>> const & VM,
    std::vector<std::vector<double>> const & VP,
    std::vector<double>& Diff,
    std::vector<long long>& Last_change,
    long long old_samples
) {

    Curr_VE_z.resize(VM.size(), std::vector<std::vector<std::vector<bool>>>(3));
    Curr_VC.resize(VM.size(), std::vector<std::vector<long long>>(3));
    for (long long l = 0; l < (long long)VM.size(); ++l)
        for (int z = 0; z < Z; ++z)
            make_edges_z(D, VM[l], z + 1, Curr_VE_z[l][z]);
    for (long long l = 0; l < (long long)VM.size(); ++l)
        for (int z = 0; z < Z; ++z)
            make_comp(Curr_VE_z[l][z], Curr_VC[l][z]);

    Curr_lklh = likelihood_log(D, VM, VP, Curr_VE_z, Curr_VC);
}

// D -- alphabet (random or from file)
// VM -- all spectra
// VP -- all intensity
void Combinatorial_method::sampling(
    std::vector<double>& D,
    std::vector<std::vector<double> >& VM,
    std::vector<std::vector<double> >& VP,
    std::vector<double>& Diff,
    std::vector<long long>& Last_change,
    long long old_samples,
    long long i
) {
    bool tmp;
    long long j = rand() % D_size;
    std::cout << "START 1st SAMPL\n";
    tmp = sampling_1(D, VM, VP, Diff, j);
    if (tmp) {
        Last_change[j] = i + old_samples;
        //std::cout << "1st change ";
    }
    std::cout << i << "\n";
    j = rand() % D_size;
    std::cout << "START 2nd SAMPL\n";
    tmp = sampling_2(D, VM, VP, j);
    if (tmp) {
        Last_change[j] = i + old_samples;
        //std::cout << "2nd change\n";
    }
}

void Combinatorial_method::sampling_for_thread(
    std::vector<double>& D,
    std::vector<std::vector<double> >& VM,
    std::vector<std::vector<double> >& VP,
    std::vector<double>& Diff,
    std::vector<long long>& Last_change,
    long long old_samples, 
    long long i, // step of sampl
    int cur_proc,
    long long j, // element from D
    double& new_d,
    double& new_lklh,
    std::vector<std::vector<std::vector<std::vector<bool>>>>& New_VE_z,
    std::vector<std::vector<std::vector<long long>>>& New_VC) 
{
    std::vector<std::vector<std::vector<std::vector<bool>>>> New_VE_z_1;
    std::vector<std::vector<std::vector<long long>>> New_VC_1;
    double new_d_1, new_lklh_1;
    bool tmp;
    tmp = sampling_1_for_thread(D, VM, VP, Diff, j, cur_proc, new_d_1, new_lklh_1, New_VE_z_1, New_VC_1);
    tmp |= sampling_2_for_thread(D, VM, VP, j, cur_proc, new_d, new_lklh, New_VE_z, New_VC);

    if (new_lklh_1 > new_lklh) {
        new_lklh = new_lklh_1;
        new_d = new_d_1;
        New_VE_z = std::move(New_VE_z_1);
        New_VC = std::move(New_VC_1);
    }
}
