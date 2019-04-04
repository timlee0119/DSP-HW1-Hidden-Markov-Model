#ifndef HMM_HEADER_
#define HMM_HEADER_

#include <vector>
#include <string>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
using namespace std;

#ifndef MAX_STATE
#	define MAX_STATE	10
#endif

#ifndef MAX_OBSERV
#	define MAX_OBSERV	26
#endif

#ifndef MAX_SEQ
#	define	MAX_SEQ		200
#endif

#ifndef MAX_LINE
#	define MAX_LINE 	256
#endif

#ifndef o
#  define o(t) seq[t] - 'A'
#endif

#ifndef Pi
#  define Pi(i) initial[i]
#endif

#ifndef A
#  define A(i, j) transition[i][j]
#endif

#ifndef B
#  define B(i, j) observation[j][i]
#endif

static FILE *open_or_die(const char*, const char*);

class HMM
{
private:
   /* data */
   char *model_name;
   int state_num;					//number of state
   int observ_num;					//number of observation
   double initial[MAX_STATE];			//initial prob.
   double transition[MAX_STATE][MAX_STATE];	//transition prob.
   double observation[MAX_OBSERV][MAX_STATE];	//observation prob.

public:
   HMM() {}
   ~HMM() {}
   int get_N() { return state_num; }
   void loadHMM(const char*);
   void dumpHMM(FILE*);
   void calculate_alpha(vector<vector<double>>&, const string&);
   void calculate_beta(vector<vector<double>>&, const string&);
   void calculate_gamma(vector<vector<double>>&,
                        const vector<vector<double>>&,
                        const vector<vector<double>>&);
   void calculate_epsilon(vector<vector<vector<double>>>&,
                          const vector<vector<double>>&,
                          const vector<vector<double>>&,
                          const string&);
   void accumulate_params(const vector<vector<double>>&,
                          const vector<vector<vector<double>>>&,
                          vector<double>&,
                          vector<vector<double>>&,
                          vector<double>&,
                          vector<vector<double>>&,
                          const string&);
   void update_model(const vector<double>& pi_numer,
                     const int sample_number,
                     const vector<vector<double>>& epsilon_accumulation,
                     const vector<double>& gamma_accumulation,
                     const vector<vector<double>>& b_numer);
   double calculate_prob(const string& seq) const;
};

void HMM::loadHMM(const char *filename)
{
   int i, j;
   FILE *fp = open_or_die( filename, "r");

   this->model_name = (char *)malloc( sizeof(char) * (strlen( filename)+1));
   strcpy( this->model_name, filename );

   char token[MAX_LINE] = "";
   while( fscanf( fp, "%s", token ) > 0 )
   {
      if( token[0] == '\0' || token[0] == '\n' ) continue;

      if( strcmp( token, "initial:" ) == 0 ){
         fscanf(fp, "%d", &this->state_num );

         for( i = 0 ; i < this->state_num ; i++ )
            fscanf(fp, "%lf", &( this->initial[i] ) );
      }
      else if( strcmp( token, "transition:" ) == 0 ){
         fscanf(fp, "%d", &this->state_num );

         for( i = 0 ; i < this->state_num ; i++ )
            for( j = 0 ; j < this->state_num ; j++ )
               fscanf(fp, "%lf", &( this->transition[i][j] ));
      }
      else if( strcmp( token, "observation:" ) == 0 ){
         fscanf(fp, "%d", &this->observ_num );

         for( i = 0 ; i < this->observ_num ; i++ )
            for( j = 0 ; j < this->state_num ; j++ )
               fscanf(fp, "%lf", &( this->observation[i][j]) );
      }
   }
}

void HMM::dumpHMM(FILE *fp)
{
   int i, j;

   //fprintf( fp, "model name: %s\n", this->model_name );
   fprintf( fp, "initial: %d\n", this->state_num );
   for( i = 0 ; i < this->state_num - 1; i++ )
      fprintf( fp, "%.5lf ", this->initial[i]);
   fprintf(fp, "%.5lf\n", this->initial[ this->state_num - 1 ] );

   fprintf( fp, "\ntransition: %d\n", this->state_num );
   for( i = 0 ; i < this->state_num ; i++ ){
      for( j = 0 ; j < this->state_num - 1 ; j++ )
         fprintf( fp, "%.5lf ", this->transition[i][j] );
      fprintf(fp,"%.5lf\n", this->transition[i][this->state_num - 1]);
   }

   fprintf( fp, "\nobservation: %d\n", this->observ_num );
   for( i = 0 ; i < this->observ_num ; i++ ){
      for( j = 0 ; j < this->state_num - 1 ; j++ )
         fprintf( fp, "%.5lf ", this->observation[i][j] );
      fprintf(fp,"%.5lf\n", this->observation[i][this->state_num - 1]);
   }
}

void HMM::calculate_alpha(vector<vector<double>>& alpha, const string& seq) {
   int T = seq.length(), N = state_num;
   /* initialize alpha */
   for (int i = 0; i < N; ++i) {
      // alpha[0][i] = initial[i] * observation[o(0)][i];
      alpha[0][i] = Pi(i) * B(i, o(0));
   }
   /* induction */
   for (int t = 0; t < T-1; ++t) {
      for (int j = 0; j < N; ++j) {
         double prev = 0.0;
         for (int i = 0; i < N; ++i) {
            prev += ( alpha[t][i] * A(i, j) );
         }
         alpha[t+1][j] = prev * B(j, o(t+1));
      }
   }
}

void HMM::calculate_beta(vector<vector<double>>& beta, const string& seq) {
   int T = seq.length(), N = state_num;
   /* initialize beta */
   for (int i = 0; i < N; ++i) {
      beta[T-1][i] = 1;
   }
   /* induction */
   for (int t = T-2; t >= 0; --t) {
      for (int i = 0; i < N; ++i) {
         for (int j = 0; j < N; ++j) {
            beta[t][i] += ( A(i, j) * B(j, o(t+1)) * beta[t+1][j] );
         }
      }
   }
}

void HMM::calculate_gamma(vector<vector<double>>& gamma,
                          const vector<vector<double>>& alpha,
                          const vector<vector<double>>& beta) {
   int T = gamma.size(), N = state_num;
   for (int t = 0; t < T; ++t) {
      double sum = 0.0;
      for (int i = 0; i < N; ++i) {
         sum += alpha[t][i] * beta[t][i];
      }
      for (int i = 0; i < N; ++i) {
         gamma[t][i] = (alpha[t][i] * beta[t][i]) / sum;
      }
   }
}

void HMM::calculate_epsilon(vector<vector<vector<double>>>& epsilon,
                            const vector<vector<double>>& alpha,
                            const vector<vector<double>>& beta,
                            const string& seq) {
   for (int t = 0; t < epsilon.size(); ++t) {
      double denominator = 0.0;
      for (int i = 0; i < state_num; ++i) {
         for (int j = 0; j < state_num; ++j) {
            denominator += alpha[t][i] * A(i, j) * B(j, o(t+1)) * beta[t+1][j];
         }
      }
      for (int i = 0; i < state_num; ++i) {
         for (int j = 0; j < state_num; ++j) {
            epsilon[t][i][j] = alpha[t][i] * A(i, j) * B(j, o(t+1)) * beta[t+1][j]
                               / denominator;
         }
      }
   }
}

void HMM::accumulate_params(const vector<vector<double>>& gamma,
                            const vector<vector<vector<double>>>& epsilon,
                            vector<double>& pi_numer,
                            vector<vector<double>>& epsilon_accumulation,
                            vector<double>& gamma_accumulation,
                            vector<vector<double>>& b_numer,
                            const string& seq) {
   int N = state_num;
   /* for pi_numer, gamma_accumulation and b_numer */
   for (int i = 0; i < N; ++i) {
      pi_numer[i] += gamma[0][i];
      for (int t = 0; t < gamma.size(); ++t) {
         gamma_accumulation[i] += gamma[t][i];
         b_numer[i][o(t)] += gamma[t][i];
      }
   }
   /* for epsilon_accumulation */
   for (int i = 0; i < N; ++i) {
      for (int j = 0; j < N; ++j) {
         for (int t = 0; t < epsilon.size(); ++t) {
            epsilon_accumulation[i][j] += epsilon[t][i][j];
         }
      }
   }
}

void HMM::update_model(const vector<double>& pi_numer,
                       const int sample_number,
                       const vector<vector<double>>& epsilon_accumulation,
                       const vector<double>& gamma_accumulation,
                       const vector<vector<double>>& b_numer) {
   /* Update pi(i) and a(i, j) */
   for (int i = 0; i < state_num; ++i) {
      Pi(i) = pi_numer[i] / sample_number;
      for (int j = 0; j < state_num; ++j) {
         A(i, j) = epsilon_accumulation[i][j] / gamma_accumulation[i];
      }
   }
   /* Update b(i, k)*/
   for (int i = 0; i < state_num; ++i) {
      for (int k = 0; k < observ_num; ++k) {
         B(i, k) = b_numer[i][k] / gamma_accumulation[i];
      }
   }
}

double HMM::calculate_prob(const string& seq) const{
   int T = seq.size(), N = state_num;
   vector<vector<double>> delta(T, vector<double>(N, 0.0));
   /* Initialization */
   for (int i = 0; i < N; ++i) {
      delta[0][i] = Pi(i) * B(i, o(0));
   }
   /* Recursion */
   for (int t = 1; t < T; ++t) {
      for (int j = 0; j < N; ++j) {
         double max = 0.0;
         for (int i = 0; i < N; ++i) {
            double res = delta[t-1][i] * A(i, j);
            if (res > max) max = res;
         }
         delta[t][j] = max * B(j, o(t));
      }
   }
   /* Termination */
   double prob = 0.0;
   for (int i = 0; i < N; ++i) {
      if (delta[T-1][i] > prob) prob = delta[T-1][i];
   }

   return prob;
}

/* ======== Local static functions ========*/
FILE *open_or_die( const char *filename, const char *ht )
{
   FILE *fp = fopen( filename, ht );
   if( fp == NULL ){
      perror( filename);
      exit(1);
   }

   return fp;
}

void read_seq_data(vector<string>& seq_data, const char* seq_file_name) {
  ifstream seq_file(seq_file_name);
  if (!seq_file) {
    cerr << "Cannot open " << seq_file_name << ", aborting...\n"; exit(-1);
  }
  string buf;
  while (getline(seq_file, buf)) {
    seq_data.push_back(buf);
  }
}

static int load_models( const char *listname, vector<HMM*> &hmms, vector<string>& model_files )
{
   FILE *fp = open_or_die( listname, "r" );

   int count = 0;
   char filename[MAX_LINE] = "";
   while( fscanf(fp, "%s", filename) == 1 ){
      hmms[count]->loadHMM(filename);
      model_files.push_back(filename);
      count ++;

      if( count >= hmms.size() ){
         // return count;
         break;
      }
   }
   fclose(fp);

   return count;
}

static void dump_models(const vector<HMM*> &hmms)
{
   for (auto it = hmms.begin(); it != hmms.end(); ++it) {
      (*it)->dumpHMM(stderr);
   }
}
/* ======== End ===========================*/

#endif