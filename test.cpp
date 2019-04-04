#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include "hmm.h"
using namespace std;

int main(int argc, const char** argv) {
  bool acc = false;
  if (argc < 4) {
    cerr << "Usage: ./test modellist.txt testing_data.txt result.txt [-a]\n"; exit(-1);
  }
  if (argc == 5 && strcmp(argv[4], "-a") == 0) {
    cout << "Computing accuracy with testing_answer.txt.\n"; acc = true;
  }

  /* initialize hmms and get modellist files*/
  vector<HMM*> hmms(5);
  for (int i = 0; i < 5; ++i) {
    hmms[i] = new HMM();
  }
  vector<string> model_files;
  load_models(argv[1], hmms, model_files);

  /* read testing data */
  vector<string> seq_data;
  read_seq_data(seq_data, argv[2]);

  /* compute best fit sequences */
  ofstream result_file(argv[3]);
  for (int i = 0; i < seq_data.size(); ++i) {
    double max_prob = 0.0;
    int best_model_index = -1;
    for (int j = 0; j < hmms.size(); ++j) {
      double res = hmms[j]->calculate_prob(seq_data[i]);
      if (res > max_prob) {
        max_prob = res;
        best_model_index = j;
      }
    }

    result_file << model_files[best_model_index] << " " << max_prob << endl;
  }

  /* computing accuracy */
  if (acc) {
    ofstream acc_file("acc.txt");
    ifstream res_file(argv[3]);
    ifstream answer_file("testing_answer.txt");
    if (!answer_file || !res_file) return 0;

    string buf1, buf2;
    double numer = 0, deno = 0;
    while (getline(res_file, buf1)) {
      getline(answer_file, buf2);
      size_t p = buf1.find(' ');
      string f1 = buf1.substr(0, p);
      if (f1 == buf2) {   // sequence model match
        numer += 1;
      }
      deno += 1;
    }
    acc_file << numer / deno;
  }

  return 0;
}
