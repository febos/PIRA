#include <string>
#include <vector>
#include <iostream>
#include <algorithm>

// https://doi.org/10.1073/pnas.77.11.6309

using namespace std;

pair<float,string> Nussinov(string seq, int minloop = 3) {


	return make_pair(1.0, "a");
	
}


int main() {

	string seq = "GGGCGGCUAGCUCAGCGGAAGAGCGCUCGCCUCACACGCGAGAGGUCGUAGGUUCAAGUCCUACGCCGCCCACCA";
	string dbn = "(((((((..((((....[..)))).(((((.......))))).....(((((..]....))))))))))))....";

	pair<float,string> scorepred = Nussinov(seq);
	float score = scorepred.first;
	string pred = scorepred.second;
	
	cout << seq << endl;
	cout << dbn << endl;
	cout << pred << endl;
	cout << score << endl;
	
	return 0;
}

