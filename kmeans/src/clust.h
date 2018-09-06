#pragma once

#include <string>
#include <vector>
#include <sstream> // std::istringstrem
#include <fstream> // std::ifstream, ofstream
#include <cmath>   // std::pow, sqrt

class Sample
{
public:
  std::string sample_name;     // name of the sample
  std::vector<float> features; // rest of the row in the input
  int cluster_id;              // the cluster that the sample belongs to
  float distance_to_centroid;  // should be updated with cluster_id

  Sample(const std::vector<std::string> &);
};

class Cluster
{
public:
  int cluster_id;                   // assigned to the samples
  std::vector<float> centroid;      // same length as Sample::features
  std::vector<string> sample_names; // samples within this cluster

  Cluster(int cluster_id, Sample &);
  void addSample(Sample);
  void removeSample(Sample);
  void updateCentroid();
};

inline std::vector<std::string>
split(const std::string &string_to_split,
      std::vector<std::string> &output_vector,
      const char delimiter);

inline std::vector<std::vector<std::string>>
read_table(std::vector<std::vector<std::string>> &matrix,
           const std::string input_file);

inline void
scale_features(std::vector<Sample> &matrix, const size_t ncol);
