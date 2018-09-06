#include "clust.h"
using std::vector, std::string;

/***************
 * class Sample *
 ****************/
Sample::Sample(const vector<string> &v)
{
    // first item should be the name, all the others should be converted to floats
    this->sample_name = v[0];
    for (size_t i = 1; i < v.size(); i++)
    {
        this->features.push_back(std::stof(v[i].c_str()));
    };
}

void Sample::set_cluster_id(size_t i)
{
    this->cluster_id = i;
}

void Sample::set_distance_to_centroid(vector<float> &centroid, bool initialize = false)
{
    if (initialize)
    {
        distance_to_centroid = 0.0;
    }
    else
    {
        float ans;
        ans = kmeans::distance(centroid, features);
        distance_to_centroid = ans;
    }
}

std::string Sample::get_sample_name()
{
    return this->sample_name;
}

std::vector<float> Sample::get_features()
{
    return this->features;
}

size_t Sample::get_cluster_id()
{
    return this->cluster_id;
}
float Sample::get_distance_to_centroid()
{
    return this->distance_to_centroid;
}

/*****************
 * class Cluster *
 *****************/
Cluster::Cluster(int cluster_id, Sample &s)
{
    this->cluster_id = cluster_id;
    this->centroid = s.get_features();
    add_sample(s);
    s.set_distance_to_centroid(centroid, true);
}

size_t Cluster::get_cluster_id()
{
    return this->cluster_id;
}
std::vector<float> Cluster::get_centroid()
{
    return this->centroid;
}

void Cluster::add_sample(Sample &s)
{
    samples.push_back(s);
    s.set_cluster_id(this->cluster_id);
}
void Cluster::remove_sample(Sample &s)
{
    for (size_t i = 0; i <= samples.size(); i++)
    {
        if (samples[i].get_sample_name() == s.get_sample_name())
        {
            samples.erase(samples.begin() + i);
            s.set_cluster_id(0);
        }
    }
}
void Cluster::update_centroid()
{
    const size_t NCOL = samples[0].get_features().size();
    centroid.resize(NCOL);
    for (size_t j = 0; j < NCOL; j++)
    {
        centroid[j] = kmeans::col_mean(samples, j);
    }
}

/********************
 * helper functions *
 ********************/
namespace kmeans
{
vector<string>
split(const string &line, vector<string> &v, const char delimiter = '\t')
{
    // Split a string and store into a vector of strings, removing the given delimiter.
    std::istringstream iss(line);
    string token;
    while (getline(iss, token, delimiter))
    {
        v.push_back(token);
    }
    return v;
}

vector<vector<string>>
read_table(vector<vector<string>> &matrix, const string input_file)
{
    // Read a file into a vector of vectors (representing a matrix).
    std::ifstream fin(input_file);

    if (fin.is_open())
    {
        string line;
        while (getline(fin, line))
        {
            vector<string> row;
            row = split(line, row);
            matrix.push_back(row);
        }
        matrix.shrink_to_fit();
        fin.close();
    }
    else
    {
        printf("Cannot open file: %s\n", input_file.c_str());
        exit(1);
    }

    return matrix;
}

float col_mean(vector<Sample> &v, size_t ncol)
{
    float feature_mean = 0.0;
    for (auto item : v)
    {
        feature_mean += item.get_features()[ncol];
    }
    return feature_mean;
}

void scale_features(vector<Sample> &v)
{
    // Scale the matrix such that each column has mean 0 and sd 1.
    const size_t NCOL = v[0].get_features().size(),
                 NROW = v.size();
    for (size_t j = 0; j < NCOL; j++)
    {
        float feature_mean, feature_sd = 0.0;
        // calculate mean
        feature_mean = col_mean(v, j);
        // calculate standard deviation
        for (size_t i = 0; i < NROW; i++)
        {
            feature_sd += std::pow(v[i].get_features()[j] - feature_mean, 2);
        }
        feature_sd = std::sqrt(feature_sd / NROW);
        // scale feature
        for (size_t i = 0; i < NROW; i++)
        {
            v[i].get_features()[j] = (v[i].get_features()[j] - feature_mean) / feature_sd;
        }
    };
}

float distance(vector<float> &v1, vector<float> &v2)
{
    float ans = 0.0;
    for (size_t i = 0; i < v1.size(); i++)
    {
        ans += std::pow(v1[i] - v2[i], 2);
    }
    ans = std::sqrt(ans);
    return ans;
}
} // namespace kmeans
