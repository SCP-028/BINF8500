#pragma once

#include <string>
#include <vector>
#include <sstream> // std::istringstrem
#include <fstream> // std::ifstream, ofstream
#include <cmath>   // std::pow, sqrt

class Sample
{
private:
  std::string sample_name;     // name of the sample
  std::vector<float> features; // rest of the row in the input
  size_t cluster_id = 0;       // the cluster that the sample belongs to
  float distance_to_centroid;  // should be updated with cluster_id

public:
  Sample(const std::vector<std::string> &);
  void set_cluster_id(size_t);
  void set_distance_to_centroid(std::vector<float> &, bool initialize);
  std::string get_sample_name();
  std::vector<float> get_features();
  size_t get_cluster_id();
  float get_distance_to_centroid();
};

class Cluster
{
private:
  size_t cluster_id;           // assigned to the samples. 0 is reserved for not belonging to any cluster.
  std::vector<float> centroid; // same length as Sample::features
  std::vector<Sample> samples; // samples within this cluster

public:
  Cluster(int cluster_id, Sample &);
  size_t get_cluster_id();
  std::vector<float> get_centroid();
  void add_sample(Sample &);
  void remove_sample(Sample &);
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
} // namespace kmeans
