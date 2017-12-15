#ifndef _BITERM_H
#define _BITERM_H

#include <cmath>
#include <string>
#include <sstream>
#include <cassert>

using namespace std;

class Biterm {
private:
  int wi;
  int wj;
  int z; // topic assignment
  int user;
  
public: 
  Biterm(int w1, int w2, int u): z(-1) {
	wi = min(w1, w2);
	wj = max(w1, w2);
        user = u;
  }

  // s format:wi   wj    z
  Biterm(string s) {
	istringstream iss(s);
	iss >> wi >> wj >> z;
  }
  
  int get_wi() const {return wi;}
  int get_wj() const {return wj;}
  
  int get_z() const {return z;}
  void set_z(int k) {z = k;}
  void reset_z() {z = -1;}

  int get_user() const {return user;}
  void set_usr(int u) { user =u;}
  
  string str() const {
	ostringstream os;
	os << wi << '\t' << wj << '\t' << '\t' << z;
	return os.str();
  }  
};

#endif
