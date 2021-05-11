#ifndef CONST_H
#define CONST_H

#include <string>
#include <fstream>
#include <iostream>
#include <thread>

//typedef unsigned long long long;

//std::ofstream cout("log.txt");

const double e = 0.01;            // epsilon
const long long D_size = 20;       // alphabet size
const long long Num_of_sampl = 2000; // number of sampling interation
const double Z = 3;               // max charge state
const int MAX_DIFF = 250;
const int NUM_OF_PROC = std::thread::hardware_concurrency();
const double INF = 1000000;
const bool IS_PARALLEL = true;

const bool FROM_FILE = 1;

const long long STATS_SIZE = 30;

const std::string name_data = "peptide_data";

const std::string name_out_file = "out_" + name_data + ".txt";
const std::string name_diff_file = "diff_" + name_data + ".txt";

#endif