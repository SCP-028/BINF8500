#pragma once

#include <string>
#include <vector>
#include <sstream> // std::stringstream
#include <fstream> // std::ifstream
#include <cmath>   // pow, sqrt, log, isnan
#include <limits>  // std::numeric_limits<float>::quiet_NaN();
#include <cstdio>
#include <iostream>

typedef std::vector<std::vector<float>> Matrix;

namespace nw
{
std::string read_fasta(const std::string input_file);
void strip_non_alphabetic(std::string &s);
void to_upper_case(std::string &s);
void reverse_string(std::string &s);
void prepare_sequence(std::string &s);
void initialize_score_matrix(Matrix &m, float gap_score);
float score_top_left(Matrix &m,
                     size_t i, size_t j,
                     std::string &seq1, std::string &seq2,
                     float match_score, float mismatch_score);
float score_left(Matrix &m,
                 size_t i, size_t j, float gap_score);
float score_top(Matrix &m,
                size_t i, size_t j, float gap_score);
float max_score(float, float, float);
} // namespace nw
