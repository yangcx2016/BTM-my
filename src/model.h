/**
 * Biterm topic model(BTM) with Gbbis sampling 
 * Author: Xiaohui Yan(xhcloud@gmail.com)
 * 2012-9-25
 */
#ifndef _MODEL_H
#define _MODEL_H

#include <vector>
#include <fstream>
#include "biterm.h"
#include "doc.h"
#include "pvec.h"
#include "pmat.h"

using namespace std;

class Model {
public:
  //vector<Biterm> bs;
  vector<Doc> blogs;

protected:
  int W;				// vocabulary size
  int K;				// number of topics
  int H;                //number of hashtags
  int n_iter;			// maximum number of iteration of Gibbs Sampling
  int save_step;

  double alpha;			// hyperparameters of p(z)
  double beta;			// hyperparameters of p(w|z)
  double gamma;
  
  // sample recorders
  Pvec<int> nb_z;	// n(b|z), size K*1
  Pmat<int> nwz;	  // n(w,z), size K*W
  Pmat<int> nhz;        //n(h,z), size K*H
  Pmat<int> ndz;        //n(d,z), size D*K
  Pvec<int> nh;         //size K*1, the number of hashtags belong to topic k
  Pvec<int> hz;         //size H*1

  Pvec<double> pw_b;   // the background word distribution  

  // If true, the topic 0 is set to a background topic that 
  // equals to the emperiacal word dsitribution. It can filter 
  // out common words
  bool has_background; 
  
  int UN;
  Pmat<int> nu_z;

public:
  Model(int K, int W, double a, double b, double g, int n_iter, int save_step,
		bool has_b = false): 
	K(K), W(W), alpha(a), beta(b), gamma(g),
	n_iter(n_iter), has_background(has_b),
	save_step(save_step) {
	pw_b.resize(W);
	nwz.resize(K, W);
    nh.resize(K);
	nb_z.resize(K);
    UN = 0;
    H = 0;
  }
  
  // run estimate procedures
  void run(string docs_pt, string res_dir, string doc_user, string doc_ht);
  
private:
  // intialize memeber varibles and biterms
  void model_init();		// load from docs
  void load_docs(string docs_pt, string user_pt, string ht_pt);
  
  // update estimate of a biterm
  void update_biterm(Biterm& bi);
  
  // reset topic proportions for biterm b
  void reset_biterm_topic(Biterm& bi, int d);
  
  // assign topic proportions for biterm b
  void assign_biterm_topic(Biterm& bi, int k, int d);
  
  // compute condition distribution p(z|b)
  void compute_pz_b(Biterm& bi, Pvec<double>& p, int d);

  void save_res(string res_dir);
  void save_pz(string pt);
  void save_pw_z(string pt);
  void save_pu_z(string dir);
  void save_ph_z(string dir);
  
  void update_docs(int d);
  void assign_hashtag_topic(int ht, int k);
  void reset_hashtag_topic(int ht);
  void compute_pz_h(int h, Pvec<double>& pz, int d);
};

#endif
