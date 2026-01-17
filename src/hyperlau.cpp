#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

// Use the Rcpp namespace to access Rcpp functions
using namespace Rcpp;

// Use the arma namespace to access Armadillo functions
using namespace arma;

#include <cmath>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <armadillo>
#include <stdio.h>
#include <numeric>
#include <fstream>
#include <unordered_map>
#include <chrono>

using namespace std;
using namespace arma;
using namespace std::chrono;
using std::ofstream;
using Matrix = std::vector<std::vector<double>>;

int globalL;
string globalstr;
vector<string> globalR;

void myexit(int code)
{
  Rcpp::stop("exiting");
}

/* this function was taken from HyperHMM https://github.com/StochasticBiology/hypercube-hmm
 A function for converting a number into a binary string.
 Input variables:
 - int n: The number you want to convert from int to binary string
 - int L: The total length of the binary string. This number needs to be greater then log_2(n)
 Output:
 - A string of length L representing the integer n as a binary string
 */
string number2binary(int n, int L){
  string binary = "";
  
  //Create the binary number (without added 0's to the left)
  for (int v = pow(2,L-1); v>=1; v /=2){
    if (n >= v){
      n -= v;
      binary = binary + "1";
    }else{
      binary = binary + "0";
    }
  }
  //Add zeros to the left such that it has length L
  while (binary.length() < L){
    binary = "0" + binary;
  }
  return binary;
}

/*this function was taken from HyperHMM https://github.com/StochasticBiology/hypercube-hmm
 A function for converting a binary string into a number.
 Input variables:
 - string bin: The binary string you want to convert to a number
 - int L: The total length of the binary string.
 Output:
 - An integer form of the binary string
 */
int binary2int(string bin, int L){
  int number = 0;
  
  //Loop through the binary number, checking at each position if the entry is a 1, if yes, add the corresponding value to the output
  for (int i = 0; i<L;i++){
    if (bin[i] == '1'){
      number +=pow(2,L-i-1);
    }
  }
  return number;
}


/* this function was taken from HyperHMM https://github.com/StochasticBiology/hypercube-hmm
 A function for creating numerous lists with information about the possible states to go to from a given state.
 Input variables (All of the vectors needs to be empty as input):
 - vector<int> n_partners: The number of possible states to go to from the state represented by the index.
 - vector<int> cumulative_partners: A vector where the element at index i represents the index where t=i starts in the vector possible
 - vector<string> partners: A vector containing all the indices (in CRS format) for states in correct order.
 - int L: The total length of the binary string.
 Output:
 - Updated versions of n_partners, cumulative_partners, and partners.
 */
void possible_transitions(vector<int>& n_partners, vector<int>& cumulative_partners, vector<int>& partners, int L){
  int end_node_int;
  string end_node;
  //Loop through all the states
  for (int i = 0; i < pow(2,L); i++){
    string vertex = number2binary(i,L);
    n_partners.push_back(0);
    for (int j = 0; j < L; j ++){
      if (vertex[j] == '0'){
        //if we find a 0 in the state it is possible to add a 1. Hence it is possible to go to the state with all equal elements except this place.
        end_node = vertex;
        end_node[j] = '1';
        end_node_int = binary2int(end_node, L);
        partners.push_back(end_node_int);
        n_partners[i]++;
      }
    }
    //Add the correct number to the cumulative vector
    if ( i == 0){
      cumulative_partners.push_back(0);
    }else{
      int c = cumulative_partners[i-1] + n_partners[i -1];
      cumulative_partners.push_back(c);
    }
  }
}

/* This function removes duplications in the input data but counts and stores how often they occurred in the original dataset
 Input variables: 
 - vector<string> original data: contains the rows of the input dataset as they are
 - vector<int> frequencies: empty vector in which the frequencies can be stored
 Output:
 - vector<string> unique_data: contains the rows of the input dataset after removing duplications
 - vector<int> frequencies: contains the number indicating how often the corresponding data point was contained in the original dataset, stored in the appropriate order to the order of the data points in the vector unique_data
 */

vector<string> red_data(vector<string> original_data, vector<int>& frequencies) {
  vector<string> unique_data;
  unordered_map<string, int> frequency;
  
  for (const auto& str : original_data) {
    if (frequency.find(str) == frequency.end()) {
      frequency[str] = 1; // First occurrence of the number
      unique_data.push_back(str);
    } else {
      ++frequency[str]; // Increment the frequency
    }
  }
  
  for (string str : unique_data) {
    frequencies.push_back(frequency[str]);
  }
  
  return unique_data;
}


/*creates a uniform rate matrix (stored as a vector), as an initial guess, with a base rate of 1 for every feature and 0 for all other pairwise interactions. Only for the models 1,2,3 and 4.
 Input:
 - vector<double> rate_matrix: zero vector of length corresponding to the number of parameters (L^model)
 - int model: the integer indicating which model is used (for the models 1,2,3 and 4)
 - int L: integer indicating the number of features under consideration (length of the binary strings)
 Output: 
 - vector<double> rate_matrix: contains a 1 in all entries corresponding to the base rate of a feature and a 0 in all other entries.
 */
void uniform_rate_matrix(vector<double>& rate_matrix, int model, int L){
  
  if (model == 1){
    for (int i = 0; i < L; i++){
      rate_matrix[i] = 1;
    }
  }else if (model == 2){
    for (int i = 0; i < L; i++){
      rate_matrix[i*L] = 1;
    }
  }else if (model == 3){
    for (int i = 0; i < L; i++){
      rate_matrix[i*L*L] = 1;
    }
  }else if (model == 4){
    for (int i = 0; i < L; i++){
      rate_matrix[i*L*L*L] = 1;
    }
  }else if (model == -1){
    Rprintf("There is something wrong with building the initial matrix for model -1 \n");
    myexit(1);
  }else{
    Rprintf("Invalid value for model \n");
    myexit(1);
  }
  
}


/* This function was taken from HyperHMM https://github.com/StochasticBiology/hypercube-hmm
 A  function that calculates the uniformized values >0 for the transition matrix, with the corresponding information about their position in the matrix.
 Input variables: (All of the vectors need to be empty as input):
 - vector A_val: contains the values (>0) that occur in the uniform transition matrix
 - vector A_row_ptr: stores the information about the corresponding row in the transition matrix
 - vector A_col_idx: stores the information about the corresponding column in the transition matrix
 - int L: the total length of the binary string
 Output:
 - Updated versions of A_val, A_row_ptr and A_col_idx
 */
void uniform_transition_matrix(arma::vec& A_val, arma::vec& A_row_ptr, arma::vec& A_col_idx, int L){
  vector<int> n_partners;
  vector<int> c_partners;
  vector<int> partners;
  possible_transitions(n_partners, c_partners, partners, L);
  int k = 0;
  int c = 0;
  //Loop through all the transition probabilities (>0) and update the specific values to be uniform
  for (int i = 0; i<pow(2,L); i++){
    int n_end_vertices = n_partners[i];
    int r = A_row_ptr(i);
    A_row_ptr(i+1) = n_end_vertices +r;
    c = 0;
    for (int j=0; j<n_end_vertices; j++){
      int j2 = partners[k];
      A_val(r+c) = 1./n_end_vertices;
      A_col_idx(r+c) = j2;
      k++;
      c++;
    }
  }
}

/*function that calculates the rate for an transition from node j to node i in the hypercube
 Input: 
 - int i: integer representation of the node in the hypercube that is the destination of the transition
 - int j: integer representation of the node in the hypercube that is the origin of the transition
 - vector<double> x_current: vector that contains the transition rates used for model 1,2,3 and 4
 - int L: integer indicating the number of features under consideration (length of the binary strings)
 - int model: the integer indicating which model is used (for the models 1,2,3 and 4)
 Output:
 - double rate: rate of the considered transition
 */
double rate(int i, int j, vector<double> x_current, int L, int model){
  
  double rate = 0;
  string i_str = number2binary(i,L);
  string j_str = number2binary(j,L);
  int n_diff = 0;
  for (int k = 0; k < L; k++){
    if (i_str[k] != j_str[k]){
      n_diff = n_diff + 1;
    }
  }
  
  int locus;
  if (n_diff != 1){
    
    rate = 0;
  }else {
    for (int k = 0; k< L; k++){
      if (i_str[k] != j_str[k]){
        locus = k;
        
        if (model == 1){
          rate = x_current[locus];
        }else if (model == 2){
          rate = x_current[locus*L + locus];
          for (int l= 0; l<L; l++){
            int digit = j_str[l]- '0';
            rate = rate + digit*x_current[l*L + locus];
          }
        }else if (model == 3){
          rate = x_current[locus*L*L+ locus*L + locus];
          for (int l = 0; l<L; l++){
            for (int m = l; m<L; m++){
              int digit_l = j_str[l] - '0';
              int digit_m = j_str[m] - '0';
              rate = rate + digit_l*digit_m*x_current[l*L*L + m*L + locus];
            }
          }
        }else if (model == 4){
          rate = x_current[locus*L*L*L+locus*L*L + locus*L + locus];
          for (int l = 0; l<L; l++){
            for (int m = l; m<L;m++){
              for (int n = m; n<L; n++){
                int digit_l = j_str[l] - '0';
                int digit_m = j_str[m] - '0';
                int digit_n = j_str[n] - '0';
                rate = rate + digit_l*digit_m*digit_n*x_current[l*L*L*L + m*L*L + n*L + locus];
              }
            }
          }
        }else{
          myexit(1);
        }
      }
    }
    
    rate =  exp(rate);
    
  }
  
  return rate;
}

/*transforms an vector with transition rates into the form of a matrix containing transition probabilities
 Input:
 - int L: integer indicating the number of features under consideration (length of the binary strings)
 - vector<double> x_current: vector that contains the transition rates used for model 1,2,3 and 4
 - int model: the integer indicating which model is used (for the models 1,2,3 and 4)
 - mat mat_temp: zero matrix of the size 2^L x 2^L
 Output:
 - mat mat_temp: matrix containing transition probabilities based on the given transition rates
 */
void build_trans_matrix(int L, vector<double> x_current, int model, mat& mat_temp){
  for (int j = 0; j< pow(2,L); j++){
    for (int i = 0; i < pow(2,L); i++){
      if (j< i){
        mat_temp(i,j) = rate(i,j,x_current,L,model);
        if (mat_temp(i,j) < 0){
          mat_temp(i,j) = 0;
        }
      }
    }	
    double sum = 0;
    for (int i = 0; i<pow(2,L); i++){
      sum = sum + mat_temp(i,j);
    }
    for (int i = 0; i<pow(2,L); i++){
      mat_temp(i,j) = mat_temp(i,j)/sum;
    }
  }
  for (int k = 0; k < pow(2,L); k++){
    mat_temp(k,pow(2,L)-1) = 0;
  }
}

/* calculates the log-likelihood of seeing the input data given a certain rate vector for the models 1,2,3 and 4
 Input: 
 -vector<string> before: vector containing all ancestor states (in string format) of the input data
 -vector<string> after: vector containing all descendant states (in string format) of the input data (in appropriate order to the vector "before")
 -vector<int> frequ: vector containing the frequencies of the occurrence of a corresponding ancestor - descendant combinations in the original dataset (in appropriate order to the vectors "before" and "after".) 
 -int L: integer indicating the number of features under consideration (length of the binary strings)
 -vector<double> x_current: vector that contains the transition rates used for model 1,2,3 and 4
 -int model: the integer indicating which model is used (for the models 1,2,3 and 4)
 Output:
 - double that indicates the log-likelihood of seeing the input data given the rate vector x_current
 */
double loglh(vector<string> before, vector<string> after, vector<int> frequ, int L, vector<double> x_current, int model){
  
  vec loglh(before.size(), fill::zeros);
  
  // 2.) Use the transition matrix state vector approach to simulate the system over timesteps t
  mat P_desc(pow(2,L),L+1, fill::zeros);
  P_desc(0,0) = 1;
  mat mat_temp(pow(2,L),pow(2,L), fill::zeros);
  build_trans_matrix(L,x_current,model,mat_temp);
  for (int t = 1; t <= L; t++){
    P_desc.col(t) = mat_temp * P_desc.col(t-1);
  }
  
  // 1.) Write down the sets of all the compatible states for the before and the after datapoint
  for (int l = 0; l< before.size(); l++){
    
    
    string DP_b = before[l];
    
    vec c_before(pow(2,L), fill::zeros);
    
    for (int i = 0; i<pow(2,L); i++){
      for (int k = 0; k<L;){
        string node = globalR[i];
        if((DP_b[k] == '?')&&(k<L-1)){
          k++;
        }else if((DP_b[k] == '?')&&(k=L-1)){
          c_before(i) = 1;
          break;
        }else if((DP_b[k] == node[k])&&(k<L-1)){
          k++;
        }else if((DP_b[k]==node[k])&&(k=L-1)){
          c_before(i) = 1;
          break;
        }else{
          break;
        }
        
      }
    }
    
    string DP_a = after[l];
    
    vec c_after(pow(2,L), fill::zeros);
    
    for (int i = 0; i<pow(2,L); i++){
      for (int k = 0; k<L;){
        string node = globalR[i];
        if((DP_a[k] == '?')&&(k<L-1)){
          k++;
        }else if((DP_a[k] == '?')&&(k=L-1)){
          c_after(i) = 1;
          break;
        }else if((DP_a[k] == node[k])&&(k<L-1)){
          k++;
        }else if((DP_a[k]==node[k])&&(k=L-1)){
          c_after(i) = 1;
          break;
        }else{
          break;
        }
      }
    }
    
    
    mat P = P_desc;
    
    // 3.) Loop through t from 0 to L
    
    vector<double> comp_probs;
    mat Q(pow(2,L),L+1, fill::zeros);
    
    for (int t = 0; t<=L; t++){
      vec P_0 = P.col(t);
      
      vec P_0_comp(pow(2,L), fill::zeros);
      for(int i = 0; i< c_before.size(); i++){
        P_0_comp[i] = c_before[i]* P_0[i];
      }
      
      if (sum(P_0_comp) != 0){
        vec Q_0_comp(pow(2,L), fill::zeros);
        vec Q_col_0 = P_0_comp;
        
        mat Q(pow(2,L),L+1, fill::zeros);
        Q.col(t) = Q_col_0;
        for (int s = t+1; s<=L; s++){
          Q.col(s) = mat_temp * Q.col(s-1);
        }
        
        for(int s = 0; s<=L; s++){
          vec Q_0 = Q.col(s);
          for(int i = 0; i< c_after.size(); i++){
            Q_0_comp[i] = c_after[i]*Q_0[i];
          }
          
          double sum_0 = sum(Q_0_comp);
          if (sum_0 != 0){
            comp_probs.push_back(sum_0);
          }
          
        }
      }
    }
    
    double lh = 0;
    for (int i = 0; i< comp_probs.size(); i++){
      lh = lh + comp_probs[i];
    }
    
    double factor = (L+1)*(L+2);
    double normlh = 2/factor*lh;
    loglh[l] = log(normlh)*frequ[l];
    
  }
  
  return sum(loglh);
  
}

/* calculates the log-likelihood of seeing the input data given a certain transition matrix for the model F (-1)
 Input: 
 -vector<string> before: vector containing all ancestor states (in string format) of the input data
 -vector<string> after: vector containing all descendant states (in string format) of the input data (in appropriate order to the vector "before")
 -vector<int> frequ: vector containing the frequencies of the occurrence of a corresponding ancestor - descendant combinations in the original dataset (in appropriate order to the vectors "before" and "after") 
 -int L: integer indicating the number of features under consideration (length of the binary strings)
 - mat x_current: matrix of size 2^L x 2^L containing transition probabilities
 Output:
 - double that indicates the log-likelihood of seeing the input data given the transition matrix x_current
 */
double loglh_minus_1(vector<string> before, vector<string> after, vector<int> frequ, int L, mat x_current){
  
  //Calculate the new likelihood given the new matrix
  vec loglh(before.size(), fill::zeros);
  
  // 2.) Use the transition matrix state vector approach to simulate the system over timesteps t
  mat P_desc(pow(2,L),L+1, fill::zeros);
  P_desc(0,0) = 1;
  for (int t = 1; t<=L; t++){
    P_desc.col(t) = x_current * P_desc.col(t-1);
  }
  
  // 1.) Write down the sets of all the compatible states for the before and the after datapoint
  for (int l = 0; l< before.size(); l++){
    
    
    string DP_b = before[l];
    
    vec c_before(pow(2,L), fill::zeros);
    
    for (int i = 0; i<pow(2,L); i++){
      for (int k = 0; k<L;){
        string node = globalR[i];
        if((DP_b[k] == '?')&&(k<L-1)){
          k++;
        }else if((DP_b[k] == '?')&&(k=L-1)){
          c_before(i) = 1;
          break;
        }else if((DP_b[k] == node[k])&&(k<L-1)){
          k++;
        }else if((DP_b[k]==node[k])&&(k=L-1)){
          c_before(i) = 1;
          break;
        }else{
          break;
        }
        
      }
    }
    
    
    string DP_a = after[l];
    
    vec c_after(pow(2,L), fill::zeros);
    
    for (int i = 0; i<pow(2,L); i++){
      for (int k = 0; k<L;){
        string node = globalR[i];
        if((DP_a[k] == '?')&&(k<L-1)){
          k++;
        }else if((DP_a[k] == '?')&&(k=L-1)){
          c_after(i) = 1;
          break;
        }else if((DP_a[k] == node[k])&&(k<L-1)){
          k++;
        }else if((DP_a[k]==node[k])&&(k=L-1)){
          c_after(i) = 1;
          break;
        }else{
          break;
        }
      }
    }
    
    
    mat P = P_desc;
    
    
    // 3.) Loop through t from 0 to L
    
    vector<double> comp_probs;
    mat Q(pow(2,L),L+1, fill::zeros);
    
    for (int t = 0; t<=L; t++){
      vec P_0 = P.col(t);
      
      vec P_0_comp(pow(2,L), fill::zeros);
      for(int i = 0; i< c_before.size(); i++){
        P_0_comp[i] = c_before[i]* P_0[i];
      }
      
      if (sum(P_0_comp) != 0){
        vec Q_0_comp(pow(2,L), fill::zeros);
        vec Q_col_0 = P_0_comp;
        
        mat Q(pow(2,L),L+1, fill::zeros);
        Q.col(t) = Q_col_0;
        for (int s = t+1; s<=L; s++){
          Q.col(s) = x_current * Q.col(s-1);
        }
        
        for(int s = 0; s<=L; s++){
          vec Q_0 = Q.col(s);
          for(int i = 0; i< c_after.size(); i++){
            Q_0_comp[i] = c_after[i]*Q_0[i];
          }
          
          double sum_0 = sum(Q_0_comp);
          if (sum_0 != 0){
            comp_probs.push_back(sum_0);
          }
          
        }
      }
    }
    
    double lh = 0;
    for (int i = 0; i< comp_probs.size(); i++){
      lh = lh + comp_probs[i];
    }
    
    double factor = (L+1)*(L+2);
    double normlh = 2/factor*lh;
    loglh[l] = log(normlh)*frequ[l];
    if(normlh <=0){
      // IGJ: this error gets thrown when we work with the loaded package and the TB L=10 data
      // but not when we run that from within the dev environment??
      Rprintf("Formal error in dataset!: normlh is %.3e\n", normlh);
      Rprintf("this is for l = %i: before %s after %s\n", l, before[l].c_str(), after[l].c_str());
      myexit(1);
    }
    
  }
  return sum(loglh);
  
  
}

/*performs the simulated annealing process for model 1 - 4 (based on a rate vector)
 Input: 
 - vector<double> x_initial: vector that contains the initial transition rates with which the simulated annealing process should be started
 - vector<double> best_mat: vector that contains the best transition rates found so far (at this point usually the same as x_initial)
 - int L: integer indicating the number of features under consideration (length of the binary strings)
 - vector<string> before: vector containing all ancestor states (in string format) of the input data
 - vector<string> after: vector containing all descendant states (in string format) of the input data (in appropriate order to the vector "before")
 - vector<int> frequ: vector containing the frequencies of the occurrence of a corresponding ancestor - descendant combinations in the original dataset (in appropriate order to the vectors "before" and "after")
 - int model: the integer indicating which model is used (for the models 1,2,3 and 4)
 - vector<double> prog_best_lik: empty vector
 - double denom: simulating annealing rate. At the end of every simulated annealing loop, the temperature is divided by this number
 Output:
 - vector<double> best_mat: vector containing the best transition rates that were found in order to optimize the likelihood function
 - vector<double> prog_best_lik: vector containing the log-likelihoods for every simulated annealing loop
 */
void simulated_annealing(vector<double> x_initial, vector<double>& best_mat, int L, vector<string> before, vector<string> after, vector<int> frequ, int model, vector<double>& prog_best_lik, double denom){
  double temp = 1;
  vector<double> x_old = x_initial;
  vector<double> x_current = x_initial;
  vector<double> x_current_temp = x_initial;
  
  double lik_initial = loglh(before, after, frequ, L, x_initial, model);
  
  double best_lik = lik_initial;
  double new_lik = lik_initial;
  double old_lik = lik_initial;
  double thresh = 0.000001;
  double up = log(thresh/temp);
  double down = log(1/denom);
  double num_it = up/down;
  while (temp > thresh){
    
    
    auto start = high_resolution_clock::now();
    
    
    //Small perturbations to get a new matrix
    for (int i = 0; i < x_old.size(); i++){
      x_current[i] = x_old[i] + (drand48() - 0.5)*0.05;
    }
    
    x_current_temp = x_current;
    
    new_lik = loglh(before, after, frequ, L, x_current, model);
    
    if (new_lik > old_lik || exp(-(old_lik -new_lik)/temp) > drand48()){
      old_lik = new_lik;
      x_old = x_current;
    }
    
    if (new_lik > best_lik){
      best_lik = new_lik;
      best_mat = x_current;
    }
    
    prog_best_lik.push_back(best_lik);
    
    if (temp == 1){
      auto stop = high_resolution_clock::now();
      auto duration = duration_cast<milliseconds>(stop - start);
      Rprintf("Estimated time: %.3e minutes\n", num_it*duration.count()/60000);
      
    }
    
    temp = temp/denom;
  }
}

/*performs the simulated annealing process for model F (-1) (based on a transition matrix)
 Input: 
 - mat x_initial: matrix that contains the initial transition probabilities with which the simulated annealing process should be started
 - mat best_mat: matrix that contains the best transition probabilities found so far (at this point usually the same as x_initial)
 - int L: integer indicating the number of features under consideration (length of the binary strings)
 - vector<string> before: vector containing all ancestor states (in string format) of the input data
 - vector<string> after: vector containing all descendant states (in string format) of the input data (in appropriate order to the vector "before")
 - vector<int> frequ: vector containing the frequencies of the occurrence of a corresponding ancestor - descendant combinations in the original dataset (in appropriate order to the vectors "before" and "after")
 - vector<double> prog_best_lik: empty vector
 - double denom: simulating annealing rate. At the end of every simulated annealing loop, the temperature is divided by this number
 Output:
 - mat best_mat: matrix containing the best transition probabilities that were found in order to optimize the likelihood function
 - vector<double> prog_best_lik: vector containing the log-likelihoods for every simulated annealing loop
 */
void simulated_annealing_minus_1(mat x_initial, mat& best_mat, int L, vector<string> before, vector<string> after, vector<int> frequ, vector<double>& prog_best_lik, double denom){
 /* for(int i = 0; i < pow(2,L); i++){
    for( int j = 0; j < pow(2,L)-1; j++){
      if(x_initial(i,j) != 0){
      x_initial(i,j) = 1;
      }
    }
  }
    for(int i = 0; i < pow(2,L); i++){
    for( int j = 0; j < pow(2,L)-1; j++){
      x_initial(i,j) /= sum(x_initial(j));
    }
  }*/
  
  double temp = 1;
  mat x_old = x_initial;
  mat x_current = x_initial;
  mat x_current_temp = x_initial;
  
  double lik_initial = loglh_minus_1(before, after, frequ, L, x_initial);
  
  double best_lik = lik_initial;
  double new_lik = lik_initial;
  double old_lik = lik_initial;
  
  double thresh = 0.000001;
  double up = log(thresh/temp);
  double down = log(1/denom);
  double num_it = up/down;
  int iteration = 0;
  
  while (temp > thresh){
    
    auto start = high_resolution_clock::now();

    //Small perturbations to get a new matrix
    for(int i = 0; i < pow(2,L); i++){
      for( int j = 0; j < pow(2,L); j++){
        if(x_old(i,j) == 0){
          x_current(i,j) = x_old(i,j);
        }else{
          x_current(i,j) = x_old(i,j) + (drand48()-0.5)*0.05;
          if( x_current(i,j) < 0){
            x_current(i,j) = 0.0001;
          }
        }
      }
    }
    
    
    x_current_temp = x_current;
    
    for(int i = 0; i < pow(2,L); i++){
      for( int j = 0; j < pow(2,L)-1; j++){
        x_current(i,j) /= sum(x_current_temp.col(j));
      }
    }
    
    
    new_lik = loglh_minus_1(before, after, frequ, L, x_current);
    
    if (new_lik > old_lik || exp(-(old_lik -new_lik)/temp) > drand48()){
      old_lik = new_lik;
      x_old = x_current;
    }
    
    if (new_lik > best_lik){
      best_lik = new_lik;
      best_mat = x_current;
    }
    
    prog_best_lik.push_back(best_lik);
    
    if (temp == 1){
      auto stop = high_resolution_clock::now();
      auto duration = duration_cast<milliseconds>(stop - start);
      Rprintf("Estimated time: %.3e minutes\n", (num_it*duration.count())/60000);
      
    }
    
    temp = temp/denom;
  }
}

List HyperLAU(NumericMatrix obs,
              Nullable<NumericMatrix> initialstates,
              NumericVector nboot,
              NumericVector model,
              NumericVector rate);

//' HyperLAU computation (Armadillo)
//'
//' @param obs Numeric matrix of observations
//' @param initialstates optional
//' @param model optional
//' @param rate optional
//' @param nboot optional
//' 
//' @return A fitted HyperLAU model, including a dataframe of transitions, vector of likelihoods, and the dimensionality of the problem
//' @examples
//' data = matrix(c(0,0,1, 0,1,1, 1,1,1), ncol=3, nrow=3)
//' fit = HyperLAU(data)
//' @export
// [[Rcpp::export]]
List HyperLAU(NumericMatrix obs,
               Nullable<NumericMatrix> initialstates = R_NilValue,
               NumericVector nboot = 0,
               NumericVector model = -1,
               NumericVector rate = 1.001)
 {
   
   Rprintf("Starting HyperLAU\n");
   
   List l;
   
   //default values for parameters
   int _bootstrap = nboot[0];
   
   int _model = model[0];
   double denom = rate[0];
   
   int n = obs.nrow();
   int m = obs.ncol();
   
   Rprintf("Found observation matrix %i x %i\n", n, m);
   
   vector<string> data(n);
   
   for (int i = 0; i < n; i++) {
     Rprintf("  getting row %i\n", i);
     
     std::string s, s1;

     for (int j = 0; j < m; j++) {
       if(obs(i, j) == 2) {
         s += "?";
       } else {
         s += std::to_string(static_cast<int>(obs(i, j)));
       }
     }
    // data[2*i+1] = s;
     if(initialstates.isUsable()) {
       NumericMatrix _initialstates(initialstates);
    //   std::string s1;
       for (int j = 0; j < m; j++) {
         if(_initialstates(i, j) == 2) {
           s1 += "?";
         } else {
           s1 += std::to_string(static_cast<int>(_initialstates(i, j)));
         }
       }
     } else {
         s1.assign(m, '0');
       }
       //data[2*i+0] = s1;
     
     data[i] = s1 + " " + s;
  }
   
   globalL = m;
   
   // store the input data in different vectors that can be used by the algorithm
   for (int i = 0; i<pow(2,globalL); i++){
     string x = number2binary(i,globalL);
     globalR.push_back(x);
   }
   
   NumericVector boot_v;
   NumericVector from_v;
   NumericVector to_v;
   NumericVector prob_v;
   NumericVector flux_v;
   
   NumericVector boot_lik_v;
   NumericVector best_lik_v;
   
   std::random_device rd;
   std::mt19937 gen(rd());
   uniform_int_distribution<int> distribution(0, n-1);
   
   for(int boot = 0; boot <= _bootstrap; boot++) {
     Rprintf("Bootstrap %i\n", boot);
   Rprintf("Checking duplicates\n");
   // check for duplications 
   vector<int> frequ;
   vector<string> reduced_data;
   
   if(boot == 0) {
     reduced_data = red_data(data,frequ);
   } else {
     string d_bt;
     vector<string> data_bt;
     for (int j = 0; j<n; j++){
       int ri = distribution(gen);
       d_bt = data[ri];
       data_bt.push_back(d_bt);
     }
     reduced_data = red_data(data_bt,frequ);
   }
   
   // setting L to the string length (= number of considered features)
   string first_row = reduced_data[0];
   
   Rprintf("Extracting words\n");
   
   string S;
   vector<string> all;
   
   for(int i =0; i < reduced_data.size(); i++){
     S = reduced_data[i];
     stringstream ss(S);
     string word;
     while (ss >> word) { // Extract word from the stream.
       all.push_back(word);
     }
   }
   
   Rprintf("Building dataset\n");
   
   string next;
   vector<string> before;
   vector<string> after;
   for(int j = 0; j < all.size(); j++){
     next = all[j];
     if (j % 2 == 0){
       before.push_back(next);
     } else{
       after.push_back(next);
     }
     Rprintf("Datapoint %i: %s\n", j, next.c_str());
   }
   
   Rprintf("Initialising matrix\n");
   
   //initialising the rate matrix
   int num_param = pow(globalL,_model);
   vector<double> rate_vec(num_param, 0.0);
   arma::mat rate_matrix(pow(2,globalL),pow(2,globalL), arma::fill::zeros);
   
   if (_model != -1){
     uniform_rate_matrix(rate_vec, _model, globalL);
   }else{
     arma::vec A_val(pow(2,globalL-1)*globalL, arma::fill::zeros);
     arma::vec A_row_ptr(pow(2,globalL)+1, arma::fill::zeros);
     arma::vec A_col_idx(pow(2,globalL-1)*globalL, arma::fill::zeros);
     
     //Calculation of the values in the transition matrix and corresponding information about their placement
     uniform_transition_matrix(A_val, A_row_ptr, A_col_idx, globalL);
     
     
     //Fitting these values into the shape of a matrix
     int c = 0;
     int j = 1;
     for (int i = 0; i < pow(2,globalL-1)*globalL; i++){
       int r = A_col_idx(i);
       if (i >= A_row_ptr(j)){
         c++;
         j++;
       }
       rate_matrix(r,c) = A_val(i);
     }
   }
   
   // transform the rate matrix into the needed form. Models 1-4 deal with a vector form of the matrix, model F (-1) deals with a matrix
   vector<double> x_initial_vec;
   vector<double> best_vec;
   mat x_initial;
   mat best_mat;
   if (_model != -1){
     x_initial_vec = rate_vec;
     best_vec = rate_vec;
   }else{
     x_initial = rate_matrix;
     best_mat = rate_matrix;
   }	
   
   //Simulated annealing 
   vector<double> prog_best_lik;
   Rprintf("Simulated annealing started \n");
   if (_model != -1){
     simulated_annealing(x_initial_vec, best_vec, globalL, before, after, frequ, _model, prog_best_lik, denom);
   }else{
     simulated_annealing_minus_1(x_initial, best_mat, globalL, before, after, frequ, prog_best_lik, denom);
   }
   
   // get the final transition matrix
   mat final_trans_matrix(pow(2,globalL), pow(2,globalL), fill::zeros);
   if (_model !=-1){
     build_trans_matrix(globalL, best_vec, _model, final_trans_matrix);
   }else{
     final_trans_matrix = best_mat;
   }
   
   // get the final probability distribution
   mat P_fin(pow(2,globalL),globalL+1, fill::zeros);
   P_fin(0,0) = 1;
   for (int t = 1; t<=globalL; t++){
     P_fin.col(t) = final_trans_matrix * P_fin.col(t-1);
   }
   
   //Calculating the fluxes
   mat fluxes(pow(2,globalL),pow(2,globalL), fill::zeros);
   for(int i = 0; i<pow(2,globalL); i++){
     for(int j = 0; j < pow(2,globalL); j++){
       double flux = sum(P_fin.row(j))* final_trans_matrix(i,j);
       fluxes(i,j) = flux;
     }
   }
   
   //creating outputs
   // 	ofstream outdata;
   // 	outdata.open("transitions_" + out_name + ".txt");
   // 	outdata << "From " << "To " << "Probability " << "Flux"<< endl;
   for(int i = 0; i<pow(2,globalL); i++){
     for(int j = 0; j < pow(2,globalL); j++){
       if (final_trans_matrix(j,i) != 0){
         boot_v.push_back(boot);
         from_v.push_back(i);
         to_v.push_back(j);
         prob_v.push_back(final_trans_matrix(j,i));
         flux_v.push_back(fluxes(j,i));
         
         // 				outdata << i << " " << j << " " << final_trans_matrix(j,i) << " " << fluxes(j,i) << endl;
       }
     }
   }
   //  	outdata.close();

   //  	ofstream liklh;
   // 	liklh.open("best_likelihood_" + out_name + ".txt");
   for(int i = 0; i < prog_best_lik.size();i++){
     boot_lik_v.push_back(boot);
     best_lik_v.push_back(prog_best_lik[i]);
     //	liklh << prog_best_lik[i] << endl;
   }
   //  	liklh.close();
   }
   
   DataFrame rdf = DataFrame::create(Named("Bootstrap") = boot_v,
                                     Named("From") = from_v,
                                     Named("To") = to_v,
                                     Named("Probability") = prob_v,
                                     Named("Flux") = flux_v);
   DataFrame rl = DataFrame::create(Named("Bootstrap") = boot_lik_v,
                                    Named("Likelihood") = best_lik_v);
   List rlist = List::create(Named("Dynamics") = rdf,
                             Named("Likelihood") = rl,
                             Named("L") = m);
   
   return(rlist);
 }
 
 
