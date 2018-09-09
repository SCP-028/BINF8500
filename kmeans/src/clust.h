#pragma once

#include <string>
#include <vector>
#include <sstream> // std::istringstrem
#include <fstream> // std::ifstream, ofstream
#include <cmath>   // std::pow, sqrt
#include <limits>  // std::numeric_limits<float>::max()

class Sample
{
private:
  size_t sample_id;            // sample index from the original matrix
  std::string sample_name;     // name of the sample
  std::vector<float> features; // rest of the row in the input
  size_t cluster_id;           // the cluster that the sample belongs to
  float distance_to_centroid;  // should be updated with cluster_id

public:
  Sample(const std::vector<std::string> &, size_t sample_id);
  void set_cluster_id(size_t);
  void set_distance_to_centroid(std::vector<float> &, bool initialize);
  size_t get_sample_id();
  std::vector<float> get_features();
  size_t get_cluster_id();
  float get_distance_to_centroid();
};

class Cluster
{
private:
  size_t cluster_id;           // 0 is reserved for not belonging to any cluster.
  size_t cluster_size;         // avoid empty clusters
  std::vector<float> centroid; // same length as Sample::features
  std::vector<Sample> samples; // samples within this cluster

public:
  Cluster(int cluster_id, Sample &);
  size_t get_cluster_id();
  std::vector<float> get_centroid();
  void add_sample(Sample &);
  void remove_sample(Sample &);
  std::vector<Sample> get_samples();
  void update_centroid();
};

// helper functions
namespace kmeans
{
std::vector<std::string> split(const std::string &string_to_split,
                               std::vector<std::string> &output_vector,
                               const char delimiter);

std::vector<std::vector<std::string>> read_table(std::vector<std::vector<std::string>> &matrix,
                                                 const std::string input_file);

float col_mean(std::vector<Sample> &, size_t column_num);

void scale_features(std::vector<Sample> &matrix);

float distance(std::vector<float> &v1, std::vector<float> &v2);

void print_matrix_row(std::vector<std::vector<std::string>> &matrix, size_t n);

void print_cluster_result(std::vector<Cluster> &res,
                          std::vector<std::vector<std::string>> &matrix);

void reset_samples(std::vector<Sample> &samples);
} // namespace kmeans
