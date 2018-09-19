#include "nw.h"

using std::string;
using std::vector;

namespace nw
{
string read_fasta(const string input_file)
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
void strip_newline(string &s)
{
    for (size_t i = s.find('\n'); i != string::npos; i = s.find('\n'))
    {
        s.erase(i);
    }
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
    strip_newline(s);
    to_upper_case(s);
    reverse_string(s);
}
void initialize_score_matrix(vector<vector<float>> &m, float gap_score)
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
float score_top_left(vector<vector<float>> &m,
                     size_t i, size_t j,
                     string &seq1, string &seq2,
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

float score_left(std::vector<std::vector<float>> &m,
                 size_t i, size_t j, float gap_score)
{
    return m[i][j - 1] + gap_score;
}

float score_top(std::vector<std::vector<float>> &m,
                size_t i, size_t j, float gap_score)
{
    return m[i - 1][j] + gap_score;
}
float max_score(float a, float b, float c)
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
void print_score_matrix(vector<vector<float>> &m,
                        string &seq1, string &seq2)
{
    printf(",Seq2");
    for (auto &c : seq2)
    {
        printf(",%c", c);
    }
    printf("\nSeq1");
    for (auto &score : m[0])
    {
        printf(",%.1f", score);
    }
    for (size_t i = 1; i <= seq1.length(); i++)
    {
        printf("\n%c", seq1[i - 1]);
        for (auto &score : m[i])
        {
            printf(",%.1f", score);
        }
    }
}
} // namespace nw
