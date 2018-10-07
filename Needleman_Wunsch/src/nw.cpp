#include "nw.h"

using std::string;
using std::vector;

namespace nw
{
string read_fasta(const string &input_file)
{
    std::ifstream fin(input_file);
    std::stringstream buffer;
    string firstline;
    if (fin.is_open())
    {
        getline(fin, firstline); // skip description line
        buffer << fin.rdbuf();
    }
    else
    {
        printf("Cannot open file: %s\n", input_file.c_str());
        exit(1);
    }
    return buffer.str();
}
void strip_non_alphabetic(string &s)
{
    size_t j = 0;
    for (auto &c : s)
    {
        if (c >= 'A' && c <= 'Z')
        {
            s[j] = c;
            j++;
        }
    }
    s = s.substr(0, j);
}
void to_upper_case(string &s)
{
    // using ASCII difference between upper and lower cases
    for (auto &c : s)
    {
        if (c >= 'a' && c <= 'z')
        {
            c -= 32;
        }
    }
}
void reverse_string(string &s)
{
    const size_t l = s.length();
    for (size_t i = 0; i < l / 2; i++)
    {
        std::swap(s[i], s[l - i - 1]);
    }
}
void prepare_sequence(string &s)
{
    to_upper_case(s);
    strip_non_alphabetic(s);
}
void initialize_score_matrix(Matrix &m, float gap_score)
{
    m[0][0] = 0;
    for (size_t i = 1; i < m.size(); i++)
    {
        m[i][0] = m[i - 1][0] + gap_score;
    }
    for (size_t j = 1; j < m[0].size(); j++)
    {
        m[0][j] = m[0][j - 1] + gap_score;
    }
}
float score_top_left(Matrix &m, size_t i, size_t j,
                     const string &seq1, const string &seq2,
                     float match_score, float mismatch_score)
{
    float score;
    if (seq1[i - 1] == seq2[j - 1])
    {
        score = m[i - 1][j - 1] + match_score;
    }
    else
    {
        score = m[i - 1][j - 1] + mismatch_score;
    }
    return score;
}
float score_left(Matrix &m, size_t i, size_t j, float gap_score)
{
    return m[i][j - 1] + gap_score;
}
float score_top(Matrix &m, size_t i, size_t j, float gap_score)
{
    return m[i - 1][j] + gap_score;
}
float max_score(const float a, const float b, const float c)
{
    float max = a;
    if (b > max)
    {
        max = b;
    }
    if (c > max)
    {
        max = c;
    }
    return max;
}
size_t count_gap(const std::string &s)
{
    size_t num = 0;
    for (auto &c : s)
    {
        if (c == '-')
        {
            num += 1;
        }
    }
    return num;
}
} // namespace nw
