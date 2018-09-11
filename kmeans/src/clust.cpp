#include "clust.h"
using std::string;
using std::vector;

/***************
 * class Sample *
 ****************/

Sample::Sample()
{
    this->sample_name = "EMPTY_SAMPLE";
    this->cluster_id = 0;
    this->distance_to_centroid = std::numeric_limits<float>::quiet_NaN();
}

Sample::Sample(const vector<string> &v, size_t sample_id)
{
    this->sample_id = sample_id;
    // first item should be the name, all the others should be converted to floats
    this->sample_name = v[0];
    for (size_t i = 1; i < v.size(); i++)
    {
        this->features.push_back(std::stof(v[i].c_str()));
    };
    this->cluster_id = 0;
    this->distance_to_centroid = std::numeric_limits<float>::quiet_NaN();
}

void Sample::set_cluster_id(size_t i)
{
    this->cluster_id = i;
}

void Sample::set_distance_to_centroid(vector<float> &centroid)
{
    this->distance_to_centroid = kmeans::distance(centroid, features);
}

size_t Sample::get_sample_id()
{
    return this->sample_id;
}

string Sample::get_sample_name()
{
    return this->sample_name;
}

vector<float> Sample::get_features()
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

void Sample::reset_sample()
{
    this->cluster_id = std::numeric_limits<unsigned>::quiet_NaN();
    this->distance_to_centroid = std::numeric_limits<float>::quiet_NaN();
}

/*****************
 * class Cluster *
 *****************/
Cluster::Cluster(size_t cluster_id)
{
    this->cluster_id = cluster_id;
    this->cluster_size = 0;
}

void Cluster::add_sample(Sample &s)
{
    s.set_cluster_id(this->cluster_id);
    if (this->centroid.size() == 0)
    {
        this->centroid = s.get_features();
    }
    s.set_distance_to_centroid(this->centroid);
    this->samples.push_back(s);
    this->cluster_size++;
}
void Cluster::remove_sample(Sample &s)
{
    for (size_t i = 0; i <= samples.size(); i++)
    {
        if (samples[i].get_sample_id() == s.get_sample_id())
        {
            s.reset_sample();
            samples.erase(samples.begin() + i);
            this->cluster_size--;
        }
    }
}
size_t Cluster::get_cluster_id()
{
    return this->cluster_id;
}
vector<float> Cluster::get_centroid()
{
    return this->centroid;
}

vector<Sample> Cluster::get_samples()
{
    return this->samples;
}

size_t Cluster::get_cluster_size()
{
    return this->samples.size();
}

void Cluster::update_centroid()
{
    const size_t NCOL = this->centroid.size();
    // centroid.resize(NCOL);
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
    for (auto &item : v)
    {
        feature_mean += item.get_features()[ncol];
    }
    feature_mean /= v.size();
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
            feature_sd += pow(v[i].get_features()[j] - feature_mean, 2.0);
        }
        feature_sd = sqrt(feature_sd / NROW);
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
        ans += pow(v1[i] - v2[i], 2);
    }
    return ans;
}

void initialize_clusters(vector<Cluster> &clusters, vector<Sample> &samples, size_t k)
{
    for (size_t cluster_id = 2; cluster_id <= k; cluster_id++)
    {
        vector<float> last_centroid = clusters.back().get_centroid();
        float max_distance = 0.0;
        size_t new_centroid_sample_id = 0;
        for (auto &sample : samples)
        {
            if (sample.get_cluster_id() == 0) // doesn't have a cluster yet
            {
                vector<float> _features = sample.get_features();
                float _distance = kmeans::distance(_features, last_centroid);
                if (_distance > max_distance)
                {
                    max_distance = _distance;
                    new_centroid_sample_id = sample.get_sample_id();
                }
            }
        }
        // create new cluster with this sample as the centroid
        clusters.emplace_back(cluster_id);
        clusters.back().add_sample(samples[new_centroid_sample_id - 1]);
    }
}

Sample furthest_sample_in_clusters(vector<Cluster> &v)
{
    Sample furthest_sample;
    for (auto &c : v)
    {
        if (c.get_cluster_size() != 0)
        {
            vector<Sample> ss = c.get_samples();
            for (auto &s : ss)
            {
                if (std::isnan(furthest_sample.get_distance_to_centroid()) ||
                    s.get_distance_to_centroid() > furthest_sample.get_distance_to_centroid())
                {
                    furthest_sample = s;
                }
            }
        }
    }
    return furthest_sample;
}
} // namespace kmeans
