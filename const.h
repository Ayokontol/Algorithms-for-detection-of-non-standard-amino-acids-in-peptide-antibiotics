#pragma once

const double e = 0.01;            // epsilon
const size_t D_size = 20;       // alphabet size
const size_t Num_of_sampl = 490; // number of sampling interation
const double Z = 3;               // max charge state
const int MAX_DIFF = 250;

const bool FROM_FILE = 1;

const size_t STATS_SIZE = 30;

const std::string name_data = "glyco_start_with_noncomb_w";

const std::string name_out_file = "out_" + name_data + ".txt";
const std::string name_diff_file = "diff_" + name_data + ".txt";
