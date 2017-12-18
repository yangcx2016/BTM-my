#ifndef _DOC_H
#define _DOC_H

#include <string>
#include <vector>
#include <list>
#include <cassert>
#include <sstream>

#include "biterm.h"

using namespace std;

class Doc {
private:
  vector<int> ws;	// word sequence
  vector<Biterm> bts; //biterms
  vector<int> hts; //hashtags
  int user;
  
public: 
  Doc(const string& s) {read_doc(s);}

  Doc(const string& wds, const string& hts, int u) {
    read_doc(wds);
    read_ht(hts);
    user = u;
  }

  int size() const {return ws.size();}
  
  const vector<int>& get_ws() const {return ws;}
  vector<Biterm>& get_bts() {return bts;}
  vector<int>& get_hts() {return hts;}
  
  const int get_w(int i) const {
	assert(i < ws.size());
	return ws[i];
  }

  const int get_ht_num() const {return hts.size();}
  const int get_size() const {return ws.size();}

  /**
  * Extract biterms from a document
  *   `win`:  window size for biterm extraction
  *   `bs`: the output biterms
  */
  void gen_biterms(int win = 15) {
	if (ws.size() < 2) return;
	
	for (int i = 0; i < ws.size()-1; ++i) 
	  for (int j = i+1; j < min(i + win, int(ws.size())); ++j) 
		bts.push_back( Biterm(ws[i], ws[j], user) );
  }

private:
  void read_doc(const string& s) {
    istringstream iss(s);
	int w;
    while (iss >> w)  ws.push_back(w);
  }
  void read_ht(const string& s) {
    istringstream iss(s);
	int w;
    while (iss >> w)  hts.push_back(w);
  }
};
  
#endif

