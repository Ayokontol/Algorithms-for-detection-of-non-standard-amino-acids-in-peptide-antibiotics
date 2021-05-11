#ifndef COMB_METHOD_H
#define COMB_METHOD_H

#include <vector>
#include <set>
#include <ctime>
#include <cmath>
#include <algorithm>
#include <map>
//#include <iostream>

#include "const.h"

class Combinatorial_method {
private:
    // find E_z (graph for charge z) for one spectrum l
    // D -- alphabet
    // M_l -- spectrum
    // z -- charge state
    // E_z -- edges, size = 0
    void make_edges_z(std::vector<double> const & D, std::vector<double> const & M_l, int z, std::vector<std::vector<bool> >& E_z) const;

    // DFS to find connected components
    // E_z -- graph
    // C -- colors
    // c -- current color
    // v -- start point
    void DFS(std::vector<std::vector<bool> > const & E_z, std::vector<long long>& C, long long c, long long v) const;

    long long make_comp(std::vector<std::vector<bool> > const & E_z, std::vector<long long>& C) const;

    // find P_l_z = Pr(D^l | E^l_z)
    // E_z -- graph for charge z
    // C -- connected components
    // P -- vector of intensity
    double Pr_l_z(std::vector<std::vector<bool>> const & E_z, std::vector<long long> const & C, std::vector<double> const & P_l) const;

    // find Pr(D^l | E^l)
    // D -- alphabet
    // M_l -- spectrum
    // P_l -- vector of intensity for spectrum l
    // VE_z -- graphs for all charges
    // VC -- connected components for all charges
    double Pr_l(
        std::vector<double> const & D,
        std::vector<double> const & M_l,
        std::vector<double> const & P_l,
        std::vector<std::vector<std::vector<bool>>> const & VE_z,
        std::vector<std::vector<long long>> const & VC
    ) const;

    double Pr_l_log(
        std::vector<double> const & D,
        std::vector<double> const & M_l,
        std::vector<double> const & P_l,
        std::vector<std::vector<std::vector<bool>>> const & VE_z,
        std::vector<std::vector<long long>> const & VC
    ) const;

    // Likelihood = Ï_l Ï_z Sum_comp Ï_edg p_i * p_j
    // D -- alphabet
    // M -- vector of vector of m/z
    // P -- vector of vector of intensity
    double likelihood(
        std::vector<double>& D,
        std::vector<std::vector<double> >& M,
        std::vector<std::vector<double> >& P,
        std::vector<std::vector<std::vector<std::vector<bool>>>>& VE_z,
        std::vector<std::vector<std::vector<long long>>>& VC
    ) const;

    // Likelihood = Ï_l Ï_z Sum_comp Ï_edg p_i * p_j
    // D -- alphabet
    // M -- vector of vector of m/z
    // P -- vector of vector of intensity
    double likelihood_log(
        std::vector<double> const & D,
        std::vector<std::vector<double> > const & M,
        std::vector<std::vector<double> > const & P,
        std::vector<std::vector<std::vector<std::vector<bool>>>> const & VE_z,
        std::vector<std::vector<std::vector<long long>>> const & VC
    ) const;


    // check P(D | D[i] = new_d) == 1
    // D -- alphabet
    // index -- index of D[i]
    // new_d -- D[i] = new_d
    bool check(std::vector<double> const & D, long long index, double new_d) const;

public:
    double Curr_lklh;

    std::vector<std::vector<std::vector<std::vector<bool>>>> Curr_VE_z; // VE_z[l][z] -- graph for scpectrum l, charge state z
    std::vector<std::vector<std::vector<long long>>> Curr_VC;           // Curr_VC[l][z] -- connected components for scpectrum l, charge state z

    // sampling with scales
    // D -- alphabet
    // VM -- vector of vector of m/z
    // VP -- vector of vector of intensity
    // i -- index of D[i]
    bool sampling_2(
        std::vector<double>& D,
        std::vector<std::vector<double> > const & VM,
        std::vector<std::vector<double> > const & VP,
        long long i);

    //sampling
    // D -- alphabet
    // VM -- vector of vector of m/z
    // VP -- vector of vector of intensity
    // Diff -- vector of differences
    // i -- index of D[i]
    bool sampling_1(
        std::vector<double>& D,
        std::vector<std::vector<double> > const & VM,
        std::vector<std::vector<double> > const & VP,
        std::vector<double> const & Diff,
        long long i);

    bool sampling_1_for_thread(
        std::vector<double> const & D,
        std::vector<std::vector<double>> const & VM,
        std::vector<std::vector<double>> const & VP,
        std::vector<double> const & Diff,
        long long i,
        int cur_proc,
        double& new_d,
        double& new_lklh,
        std::vector<std::vector<std::vector<std::vector<bool>>>>& New_VE_z,
        std::vector<std::vector<std::vector<long long>>>& New_VC) const;

    bool sampling_2_for_thread(
        std::vector<double>& D,
        std::vector<std::vector<double> > const& VM,
        std::vector<std::vector<double> > const& VP,
        long long i,
        int cur_proc,
        double& new_d,
        double& new_lklh,
        std::vector<std::vector<std::vector<std::vector<bool>>>>& New_VE_z,
        std::vector<std::vector<std::vector<long long>>>& New_VC
    ) const;

    void sampling_for_thread(
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
        std::vector<std::vector<std::vector<long long>>>& New_VC);

    void prepare_to_sampling(
        std::vector<double> const & D,
        std::vector<std::vector<double>> const & VM,
        std::vector<std::vector<double>> const & VP,
        std::vector<double>& Diff,
        std::vector<long long>& Last_change,
        long long old_samples);
    
    void sampling(
        std::vector<double>& D,
        std::vector<std::vector<double> >& VM,
        std::vector<std::vector<double> >& VP,
        std::vector<double>& Diff,
        std::vector<long long>& Last_change,
        long long old_samples,
        long long i);

    //void sampling_for_thread(
    //    std::vector<double>& D,
    //    std::vector<std::vector<double> >& VM,
    //    std::vector<std::vector<double> >& VP,
    //    std::vector<double>& Diff,
    //    std::vector<long long>& Last_change,
    //    long long old_samples);
};

#endif
