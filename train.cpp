#include <iostream>
#include <vector>
#include <string>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <cstring>
#include "hmm.h"
using namespace std;

static void train(HMM* hmm, int iteration, const char* seq_file_name);
static void dump_model(HMM* hmm, const char* dump_file_name);

int main(int argc, const char** argv) {
  if (argc != 5) {
    cerr << "Usage: ./train iteration model_init.txt seq_model_x.txt model_x.txt\n"; exit(-1);
  }
  int iteration = strtol(argv[1], nullptr, 10);
  if (iteration <= 0) {
    cerr << "Iteration has to be greater than zero.\n"; exit(-1);
  }

  HMM* hmm_initial = new HMM();
  hmm_initial->loadHMM(argv[2]);
  
  train(hmm_initial, iteration, argv[3]);
  dump_model(hmm_initial, argv[4]);

  return 0;
}

void train(HMM* hmm, int iteration, const char* seq_file_name) {
  vector<string> seq_data;
  read_seq_data(seq_data, seq_file_name);
  int T = seq_data.at(0).length(), N = hmm->get_N();

  cout << "Training with sequence: " << seq_file_name << "...\n";
  for (int it = 0; it < iteration; ++it) {
    cout << "Iteration: " << it+1 << "\n";
    // calculate alpha, beta, gamma, epsilon sequence by sequence
    // then accumulate epsion and gamma
    // then finally update parameter of hmm
    vector<double> pi_numer(N, 0.0);
    vector<vector<double>> epsilon_accumulation(N, vector<double>(N, 0.0));
    vector<double> gamma_accumulation(N, 0.0);
    vector<vector<double>> b_numer(N, vector<double>(MAX_OBSERV, 0.0));

    for (int i = 0, n = seq_data.size(); i < n; ++i) {
      vector<vector<double>> alpha(T, vector<double>(N, 0.0));
      hmm->calculate_alpha(alpha, seq_data[i]);
      vector<vector<double>> beta(T, vector<double>(N, 0.0));
      hmm->calculate_beta(beta, seq_data[i]);
      vector<vector<double>> gamma(T, vector<double>(N, 0.0));
      hmm->calculate_gamma(gamma, alpha, beta);
      vector<vector<vector<double>>> epsilon(T-1, vector<vector<double>>(N, vector<double>(N, 0.0)));
      hmm->calculate_epsilon(epsilon, alpha, beta, seq_data[i]);

      hmm->accumulate_params(gamma, epsilon, pi_numer, epsilon_accumulation,
                             gamma_accumulation, b_numer, seq_data[i]);
    }

    hmm->update_model(pi_numer, seq_data.size(), epsilon_accumulation, gamma_accumulation, b_numer);
  }
}

static void dump_model(HMM* hmm, const char* dump_file_name) {
  FILE* fp = fopen(dump_file_name, "w+");
  hmm->dumpHMM(fp);
  fclose(fp);
}
