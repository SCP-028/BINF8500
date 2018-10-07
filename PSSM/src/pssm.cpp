#include "pssm.h"

using std::string;
using std::vector;

namespace preprocess
{
string read_DNA(const string &DNA_path)
{
    std::ifstream fin(DNA_path);
    std::stringstream buffer;
    string firstline;
    if (fin.is_open())
    {
        getline(fin, firstline);
        buffer << fin.rdbuf();
        fin.close();
        return buffer.str();
    }
    else
    {
        printf("Cannot open file: %s\n", DNA_path.c_str());
        exit(1);
    }
}

Svec read_motifs(const string &motif_path)
{
    std::ifstream fin(motif_path);

    if (fin.is_open())
    {
        Svec motifs;
        motifs.reserve(1000);
        string line;
        while (getline(fin, line))
        {
            motifs.emplace_back(line);
        }
        fin.close();
        motifs.shrink_to_fit();
        return motifs;
    }
    else
    {
        printf("Cannot open file:%s\n", motif_path.c_str());
        exit(1);
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
void normalize_sequence(string &s)
{
    size_t i = 0;
    for (auto &c : s)
    {
        if (c == 'A' || c == 'C' || c == 'G' || c == 'T')
        {
            s[i] = c;
            i++;
        }
    }
    s = s.substr(0, i);
}
} // namespace preprocess

PSSM::PSSM(const Svec &motifs, const string &DNA)
{
    const size_t motif_num = motifs.size();
    assert(motif_num > 0);
    m_motif_len = motifs[0].size();
    m_score_matrix.resize(m_motif_len);
    for (size_t j = 0; j < m_motif_len; j++)
    {
        string col(motif_num, 'x');
        for (size_t i = 0; i < motif_num; i++)
        {
            col[i] = motifs[i][j];
        }
        m_score_matrix[j] = count_char_occurrence(col);
    }
    m_background_prob = count_char_occurrence(DNA);
    const size_t DNA_len = DNA.size();
    assert(DNA_len > 0);
    for (auto &score : m_background_prob)
    {
        score /= DNA_len;
    }
}

vector<float> PSSM::count_char_occurrence(const string &s)
{
    vector<float> ans{0, 0, 0, 0}; // 4 columns: A C G T
    for (auto &c : s)
    {

        switch (c)
        {
        case 'A':
            ans[0] += 1;
            break;
        case 'C':
            ans[1] += 1;
            break;
        case 'G':
            ans[2] += 1;
            break;
        case 'T':
            ans[3] += 1;
            break;

        default:
            break;
        }
    }
    return ans;
}

void PSSM::print_score_matrix(const unsigned rounding)
{
    printf("\n%10s%10c%10c%10c%10c\n", "Position", 'A', 'C', 'G', 'T');
    for (size_t i = 0; i < m_motif_len; i++)
    {
        printf("%10zu", i + 1);
        for (auto &score : m_score_matrix[i])
        {
            printf("%10.*f", rounding, score);
        }
        printf("\n");
    }
}

void PSSM::add_pseudocount(const float pc)
{
    for (auto &row : m_score_matrix)
    {
        for (auto &score : row)
        {
            score += pc;
        }
    }
}

void PSSM::convert_to_score_matrix()
{
    for (auto &row : m_score_matrix)
    {
        float total_score = 0.0;
        for (auto &freq_pseudo : row)
        {
            total_score += freq_pseudo;
        }
        const size_t row_size = row.size();
        for (size_t i = 0; i < row_size; i++)
        {
            row[i] = std::log2(row[i] / (total_score * m_background_prob[i]));
        }
    }
}

void PSSM::generate_reverse_matrix()
{
    const size_t ncol = m_score_matrix[0].size();

    // get m_reverse_matrix to a proper size
    m_reverse_matrix.resize(m_motif_len);
    for (auto &row : m_reverse_matrix)
    {
        row.resize(ncol);
    }

    // A <-> T and C <-> G
    for (size_t i = 0; i < m_motif_len; i++)
    {
        for (unsigned j = 0; j < ncol; j++)
        {
            m_reverse_matrix[m_motif_len - i - 1][j] = m_score_matrix[i][ncol - j - 1];
        }
    }
}

float PSSM::calc_score_for_forward(const std::string &s, size_t istart)
{
    float total_score = 0.0;
    for (size_t i = 0; i < m_motif_len; i++)
    {
        switch (s[i + istart])
        {
        case 'A':
            total_score += m_score_matrix[i][0];
            break;
        case 'C':
            total_score += m_score_matrix[i][1];
            break;
        case 'G':
            total_score += m_score_matrix[i][2];
            break;
        case 'T':
            total_score += m_score_matrix[i][3];
            break;

        default:
            break;
        };
    }
    return total_score;
}

float PSSM::calc_score_for_reverse(const std::string &s, size_t istart)
{
    float total_score = 0.0;
    for (size_t i = 0; i < m_motif_len; i++)
    {
        switch (s[i + istart])
        {
        case 'A':
            total_score += m_reverse_matrix[i][0];
            break;
        case 'C':
            total_score += m_reverse_matrix[i][1];
            break;
        case 'G':
            total_score += m_reverse_matrix[i][2];
            break;
        case 'T':
            total_score += m_reverse_matrix[i][3];
            break;

        default:
            break;
        };
    }
    return total_score;
}

string PSSM::generate_reverse_strand(const std::string &s)
{
    const size_t s_size = s.size() - 1;
    string ans(s_size + 1, 'x');
    for (size_t i = 0; i <= s_size; i++)
    {
        switch (s[i])
        {
        case 'A':
            ans[s_size - i] = 'T';
            break;
        case 'C':
            ans[s_size - i] = 'G';
            break;
        case 'G':
            ans[s_size - i] = 'C';
            break;
        case 'T':
            ans[s_size - i] = 'A';
            break;

        default:
            break;
        }
    }
    return ans;
}
