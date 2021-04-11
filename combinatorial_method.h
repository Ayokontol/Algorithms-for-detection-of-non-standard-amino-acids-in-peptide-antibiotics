#pragma once

#include <vector>
#include <set>
#include <ctime>
#include <cmath>
#include <algorithm>
#include <map>
#include <iostream>

#include "const.h"

class Combinatorial_method {
private:
    double Curr_lklh;

    std::vector<std::vector<std::vector<std::vector<bool>>>> Curr_VE_z; // VE_z[l][z] -- graph for scpectrum l, charge state z
    std::vector<std::vector<std::vector<int>>> Curr_VC;                 // Curr_VC[l][z] -- connected components for scpectrum l, charge state z

    // find E_z (graph for charge z) for one spectrum l
    // D -- alphabet
    // M_l -- spectrum
    // z -- charge state
    // E_z -- edges, size = 0
    void make_edges_z(std::vector<double>& D, std::vector<double>& M_l, int z, std::vector<std::vector<bool> >& E_z);

    // DFS to find connected components
    // E_z -- graph
    // C -- colors
    // c -- current color
    // v -- start point
    void DFS(std::vector<std::vector<bool> >& E_z, std::vector<int>& C, int c, size_t v);

    size_t make_comp(std::vector<std::vector<bool> >& E_z, std::vector<int>& C);

    // find P_l_z = Pr(D^l | E^l_z)
    // E_z -- graph for charge z
    // C -- connected components
    // P -- vector of intensity
    double Pr_l_z(std::vector<std::vector<bool> >& E_z, std::vector<int>& C, std::vector<double>& P_l);

    // find Pr(D^l | E^l)
    // D -- alphabet
    // M_l -- spectrum
    // P_l -- vector of intensity for spectrum l
    // VE_z -- graphs for all charges
    // VC -- connected components for all charges
    double Pr_l(
        std::vector<double>& D,
        std::vector<double>& M_l,
        std::vector<double>& P_l,
        std::vector<std::vector<std::vector<bool>>>& VE_z,
        std::vector<std::vector<int>>& VC
    );

    double Pr_l_log(
        std::vector<double>& D,
        std::vector<double>& M_l,
        std::vector<double>& P_l,
        std::vector<std::vector<std::vector<bool>>>& VE_z,
        std::vector<std::vector<int>>& VC
    );

    // Likelihood = Ï_l Ï_z Sum_comp Ï_edg p_i * p_j
    // D -- alphabet
    // M -- vector of vector of m/z
    // P -- vector of vector of intensity
    double likelihood(std::vector<double>& D, std::vector<std::vector<double> >& M, std::vector<std::vector<double> >& P);

    // Likelihood = Ï_l Ï_z Sum_comp Ï_edg p_i * p_j
    // D -- alphabet
    // M -- vector of vector of m/z
    // P -- vector of vector of intensity
    double likelihood_log(std::vector<double>& D, std::vector<std::vector<double> >& M, std::vector<std::vector<double> >& P);


    // check P(D | D[i] = new_d) == 1
    // D -- alphabet
    // index -- index of D[i]
    // new_d -- D[i] = new_d
    bool check(std::vector<double>& D, size_t index, double new_d);

public:
    // sampling with scales
    // D -- alphabet
    // VM -- vector of vector of m/z
    // VP -- vector of vector of intensity
    // i -- index of D[i]
    bool sampling_2(std::vector<double>& D, std::vector<std::vector<double> >& VM, std::vector<std::vector<double> >& VP, size_t i);

    //sampling
    // D -- alphabet
    // VM -- vector of vector of m/z
    // VP -- vector of vector of intensity
    // Diff -- vector of differences
    // i -- index of D[i]
    bool sampling_1(std::vector<double>& D, std::vector<std::vector<double> >& VM, std::vector<std::vector<double> >& VP, std::vector<double>& Diff, size_t i);

    void prepare_to_sampling(
        std::vector<double>& D,
        std::vector<std::vector<double> >& VM,
        std::vector<std::vector<double> >& VP,
        std::vector<double>& Diff,
        std::vector<size_t>& Last_change,
        size_t old_samples);
    
    void sampling(
        std::vector<double>& D,
        std::vector<std::vector<double> >& VM,
        std::vector<std::vector<double> >& VP,
        std::vector<double>& Diff,
        std::vector<size_t>& Last_change,
        size_t old_samples,
        size_t i);
};
