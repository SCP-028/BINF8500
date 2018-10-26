#pragma once

#include <cmath>
#include <cstdio>
#include <fstream>
#include <random>
#include <sstream>
#include <string>
#include <vector>
using namespace std;

const static float PSEUDOCOUNT = 0.25;

template <typename T>
using Matrix = vector<vector<T>>;

namespace infiles {
inline Matrix<string> read_fasta(const string &input_file) {
    std::ifstream fin(input_file);
    std::stringstream buffer;
    if (fin.is_open()) {
        buffer << fin.rdbuf();
    } else {
        printf("Cannot open file: %s\n", input_file.c_str());
        exit(1);
    }
    string line;
    getline(buffer, line);  // first sequence
    vector<string> seq{line, ""};
    Matrix<string> ans;
    ans.reserve(30);
    while (getline(buffer, line)) {
        if (line[0] == '>') {
            ans.emplace_back(seq);
            seq[0] = line;
            seq[1] = "";
        } else {
            seq[1] += line;
        }
    }
    ans.emplace_back(seq);
    ans.shrink_to_fit();
    return ans;
}
}  // namespace infiles

namespace pssm {
inline Matrix<float> init_freq_matrix(const float pc, const float motif_len) {
    // Creates a motif_len by 4 matrix filled with pseudo-counts
    // as initial values
    Matrix<float> ans(motif_len, vector<float>(4, pc));
    return ans;
}
inline void count_char_occurrence(const string &s, vector<float> &bg_freq) {
    for (auto &c : s) {
        switch (c) {
            case 'A':
                bg_freq[0] += 1;
                break;
            case 'C':
                bg_freq[1] += 1;
                break;
            case 'G':
                bg_freq[2] += 1;
                break;
            case 'T':
                bg_freq[3] += 1;
                break;

            default:
                break;
        }
    }
    for (auto &freq : bg_freq) {
        freq /= s.size();
    }
}
inline void count_char_occurrence(const string &s, Matrix<float> &freq_matrix) {
    //     pos      A  C  G  T
    //      0       pc pc pc pc
    //      1       pc pc pc pc
    //     ...      ...........
    // motif_len-1  pc pc pc pc
    const size_t motif_len = s.size();
    for (size_t i = 0; i < motif_len; i++) {
        switch (s[i]) {
            case 'A':
                freq_matrix[i][0] += 1;
                break;
            case 'C':
                freq_matrix[i][1] += 1;
                break;
            case 'G':
                freq_matrix[i][2] += 1;
                break;
            case 'T':
                freq_matrix[i][3] += 1;
                break;

            default:
                break;
        }
    }
}
inline void convert_to_pssm(const vector<float> &bg_freq,
                            Matrix<float> &freq_matrix) {
    // not taking log because can't have negative scores
    for (auto &row : freq_matrix) {
        float total_freq = 0.0;
        for (auto &ele : row) {
            total_freq += ele;
        }
        for (int i = 0; i < 4; i++) {
            row[i] /= total_freq * bg_freq[i];
        }
    }
}
inline float calc_score(const Matrix<float> &pssm, const string &s,
                        const size_t start_pos, const size_t motif_len) {
    float total_score = 0.0;
    for (size_t i = 0; i < motif_len; i++) {
        switch (s[i + start_pos]) {
            case 'A':
                total_score += pssm[i][0];
                break;
            case 'C':
                total_score += pssm[i][1];
                break;
            case 'G':
                total_score += pssm[i][2];
                break;
            case 'T':
                total_score += pssm[i][3];
                break;

            default:
                break;
        };
    }
    return total_score;
}
inline Matrix<float> build_pssm(const Matrix<string> &fasta,
                                const size_t motif_len, const size_t seq_i,
                                vector<size_t> &motif_positions) {
    vector<float> bg_freq(4, 0.0);
    Matrix<float> score_matrix = pssm::init_freq_matrix(PSEUDOCOUNT, motif_len);
    const size_t fasta_size = fasta.size();
    for (size_t i = 0; i < fasta_size; i++) {
        // use left-out sequence to calculate background
        // frequency, and use other motifs to build PSSM
        if (i == seq_i) {
            pssm::count_char_occurrence(fasta[i][1], bg_freq);
        } else {
            string motif_seq =
                fasta[i][1].substr(motif_positions[i], motif_len);
            pssm::count_char_occurrence(motif_seq, score_matrix);
        }
    }
    pssm::convert_to_pssm(bg_freq, score_matrix);
    return score_matrix;
}
}  // namespace pssm

namespace gibbs {
inline vector<size_t> init_motif_positions(const Matrix<string> &seqs,
                                           const size_t motif_len,
                                           mt19937 seed) {
    vector<size_t> motif_positions;
    motif_positions.reserve(seqs.size());
    for (auto &seq : seqs) {
        uniform_int_distribution<size_t> range(0, seq[1].size() - motif_len);
        motif_positions.emplace_back(range(seed));
    }
    return motif_positions;
}
vector<float> calc_scores_in_seq(const Matrix<string> &fasta,
                                 const size_t seq_i, const float motif_len,
                                 const Matrix<float> &score_matrix) {
    const size_t motif_num = fasta[seq_i][1].size() - motif_len;
    vector<float> motif_scores(motif_num, 0.0);
    for (size_t _pos = 0; _pos < motif_num; _pos++) {
        motif_scores[_pos] =
            pssm::calc_score(score_matrix, fasta[seq_i][1], _pos, motif_len);
    }
    return motif_scores;
}
size_t update_position(const Matrix<string> &fasta, const size_t motif_len,
                       mt19937 seed, vector<size_t> &motif_positions,
                       vector<float> &fasta_scores) {
    size_t num_changes = 0;
    const size_t seq_num = fasta.size();
    for (size_t i = 0; i < seq_num; i++) {
        // Use all words except seq[i] to build the PSSM
        Matrix<float> score_matrix =
            pssm::build_pssm(fasta, motif_len, i, motif_positions);

        // Use PSSM to assign score to each position in seq 1
        vector<float> motif_scores =
            gibbs::calc_scores_in_seq(fasta, i, motif_len, score_matrix);

        // Select new location in seq 1 with probability proportional to the
        // previous scores using a weighted probability distribution
        std::discrete_distribution<size_t> generator(motif_scores.begin(),
                                                     motif_scores.end());
        size_t current_pos = motif_positions[i], new_pos = generator(seed);
        if (new_pos != current_pos) {
            num_changes++;
            motif_positions[i] = new_pos;
            fasta_scores[i] =
                pssm::calc_score(score_matrix, fasta[i][1], new_pos, motif_len);
        }
    }
    return num_changes;
}
inline void update_score(const vector<float> &fasta_scores,
                         const vector<size_t> &motif_positions,
                         float &max_score, vector<size_t> &final_positions) {
    float total_score = 0.0;
    for (auto &ele : fasta_scores) {
        total_score += ele;
    }
    if (total_score > max_score) {
        max_score = total_score;
        final_positions = motif_positions;
    }
}
}  // namespace gibbs
