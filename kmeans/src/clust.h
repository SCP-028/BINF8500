#pragma once

#include <cassert>
#include <cmath>    // pow, sqrt, log, isnan
#include <cstdio>   // printf
#include <fstream>  // std::ifstream
#include <limits>   // std::numeric_limits<float>::quiet_NaN();
#include <random>
#include <sstream>  // std::istringstrem
#include <string>
#include <vector>

class Sample {
   private:
    size_t sample_id;             // sample index from the original matrix
    std::string sample_name;      // name of the sample
    std::vector<float> features;  // rest of the row in the input
    size_t cluster_id;            // the cluster that the sample belongs to
    float distance_to_centroid;   // should be updated with cluster_id

   public:
    Sample();
    Sample(const std::vector<std::string> &, size_t sample_id);
    void set_cluster_id(size_t);
    void set_feature(size_t i, float feature);
    void set_distance_to_centroid(std::vector<float> &);
    size_t get_sample_id();
    std::string get_sample_name();
    std::vector<float> get_features();
    size_t get_cluster_id();
    float get_distance_to_centroid();
    void reset_sample();
};

class Cluster {
   private:
    size_t cluster_id;    // 0 is reserved for not belonging to any cluster.
    size_t cluster_size;  // avoid empty clusters
    std::vector<float> centroid;  // same length as Sample::features
    std::vector<Sample> samples;  // samples within this cluster

   public:
    Cluster(size_t cluster_id);
    void add_sample(Sample &);
    void remove_sample(Sample &);
    size_t get_cluster_id();
    std::vector<float> get_centroid();
    std::vector<Sample> get_samples();
    size_t get_cluster_size();
    void update_centroid();
};

// helper functions
namespace kmeans {
std::vector<std::string> split(const std::string &string_to_split,
                               std::vector<std::string> &output_vector,
                               const char delimiter);

std::vector<std::vector<std::string>> read_table(
    std::vector<std::vector<std::string>> &matrix,
    const std::string input_file);

float col_mean(std::vector<Sample> &, size_t column_num);

void scale_features(std::vector<Sample> &matrix);

float distance(std::vector<float> &v1, std::vector<float> &v2);

void initialize_clusters(std::vector<Cluster> &clusters,
                         std::vector<Sample> &samples, size_t k);

Sample furthest_sample_in_clusters(std::vector<Cluster> &);
}  // namespace kmeans
