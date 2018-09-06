#include "clust.h"
using std::vector, std::string;

Sample::Sample(const vector<string> &v)
{
    // first item should be the name, all the others should be converted to floats
    this->sample_name = v[0];
    for (size_t i = 1; i < v.size(); i++)
    {
        this->features.push_back(std::stof(v[i].c_str()));
    };
}

Cluster::Cluster(int cluster_id, Sample &s)
{
    this->cluster_id = cluster_id;
    this->centroid = s.features;
    this->sample_names.push_back(s.sample_name);
};
//   void addSample(Sample);
//   void removeSample(Sample);
//   void updateCentroid();

inline vector<string>
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

inline vector<vector<string>>
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

inline void
scale_features(vector<Sample> &v, const size_t n)
{
    // Scale the matrix such that each column has mean 0 and sd 1.
    for (size_t j = 0; j < n; j++)
    {
        float feature_mean = 0.0, feature_sd = 0.0;
        const size_t feature_num = v.size();
        // calculate mean
        for (size_t i = 0; i < feature_num; i++)
        {
            feature_mean += v[i].features[j];
        }
        feature_mean /= feature_num;
        // calculate standard deviation
        for (size_t i = 0; i < feature_num; i++)
        {
            feature_sd += std::pow(v[i].features[j] - feature_mean, 2);
        }
        feature_sd = std::sqrt(feature_sd / feature_num);
        // scale feature
        for (size_t i = 0; i < feature_num; i++)
        {
            v[i].features[j] = (v[i].features[j] - feature_mean) / feature_sd;
        }
    };
}
