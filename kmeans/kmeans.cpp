#include <cstdio>   // printf
#include <cstdlib>  // exit
#include <fstream>  // std::ifstream, ofstream
#include <iterator> // std::ostream_iterator
#include <math.h>   //std::pow, sqrt
#include <sstream>  // std::istringstrem
#include <string>   // std::string
#include <vector>   // std::vector, swap
using namespace std;

class Sample
{
  public:
    string name;
    vector<float> features;
    int cluster;

    Sample(const vector<string> &v)
    {
        // first item should be the name, all the others should be converted to floats
        this->name = v[0];
        for (size_t i = 1; i < v.size(); i++)
        {
            this->features.push_back(stof(v[i].c_str()));
        };
    }
};

vector<string> split(const string &line, vector<string> &v, const char delimiter = '\t')
{
    /* Split a string into a vector of strings, removing the given delimiter. */
    istringstream iss(line);
    string token;
    while (getline(iss, token, delimiter))
    {
        v.push_back(token);
    }
    return v;
}

vector<vector<string>> read_table(vector<vector<string>> &matrix, const string input_file)
{
    /* Read a file into a vector of vectors (representing a matrix). */
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

void scale_features(vector<Sample> &v, const size_t n)
{
    for (size_t j = 0; j < n; j++)
    {
        float feature_mean, feature_sd = 0.0;
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
            feature_sd += pow(v[i].features[j] - feature_mean, 2);
        }
        feature_sd = sqrt(feature_sd / feature_num);
        // scale feature
        for (size_t i = 0; i < feature_num; i++)
        {
            v[i].features[j] = (v[i].features[j] - feature_mean) / feature_sd;
        }
    };
}

int main(int argc, char **argv)
{
    if (argc != 2)
    {
        printf("[Error] %s takes 1 arguement, but %d were given.\n\n"
               "Usage: %s <input_file>\n\n"
               "This program implements the k-means clustering algorithm.\n\n"
               "The appropriate value for k is chosen automatically using the Bayesian information criterion.\n",
               argv[0], argc - 1, argv[0]);
        exit(1);
    }
    // Read file into vector of vectors
    vector<vector<string>> matrix;
    matrix.reserve(1000);
    matrix = read_table(matrix, argv[1]);
    const size_t NROW = matrix.size();
    const size_t NCOL = matrix[0].size() - 1; // first column is used as rowname later
    // Convert rows into objects for better readability
    vector<Sample> samples;
    for (size_t i = 1; i < NROW; i++) // skip header row
    {
        Sample g(matrix[i]);
        samples.push_back(g);
    }
    // Normalize data before clustering -> mean 0 and sd 1 on each column
    scale_features(samples, NCOL);

    // Pick 1st centroid randomly, then pick the point that's furthest,
    // Assign each point to the closest centroid
    // Update the location of the centroid by averaging all the points within the cluster

    // Check the Optimal Growth Temperature in the output (GC content & temp)
    for (auto g : samples)
    {
        printf("%s,", g.name.c_str());
        for (auto i : g.features)
        {
            printf("%.2f,", i);
        }
        printf("\n");
    }
    return 0;
}
