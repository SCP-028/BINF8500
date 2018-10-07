#pragma once

#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <cassert>

using Svec = std::vector<std::string>;
template <typename T>
using Matrix = std::vector<std::vector<T>>;

namespace preprocess
{
Svec read_motifs(const std::string &s);     // flat file
std::string read_DNA(const std::string &s); // fasta file
void normalize_sequence(std::string &s);
void to_upper_case(std::string &s);
} // namespace preprocess

class PSSM
{
private:
  size_t m_motif_len;
  Matrix<float> m_score_matrix; // columns A, C, G, T
  Matrix<float> m_reverse_matrix;

public:
  std::vector<float> m_background_prob;
  PSSM(const Svec &motifs, const std::string &DNA);
  std::vector<float> count_char_occurrence(const std::string &s);
  void print_score_matrix(const unsigned rounding);
  void construct_freq_matrix(const Svec &v);
  void add_pseudocount(const float pc);
  void convert_to_score_matrix();
  void generate_reverse_matrix();
  float calc_score_for_forward(const std::string &s, size_t i);
  float calc_score_for_reverse(const std::string &s, size_t i);
  std::string generate_reverse_strand(const std::string &s);
};
