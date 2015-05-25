
#include "headers.h"
#include "config.h"
#include "IO.h"

/*
 *
 * FILTER REF
 *
 */

void writeNewRef(Config &config, string file, string file2, vector<int> &coverage, string &part_genome, string &whole_genome, vector<int> &positions) {

	ofstream out(file.c_str());
	ofstream out2(file2.c_str());
	vector<int> starts;
	vector<int> stops;
	unsigned int start = 0;
	unsigned int stop = 0;
	unsigned int position = 0;
	for(unsigned int i = 0; i < coverage.size(); i++) {
		if(coverage[i] > config.COV_THRES) {
			start = i;
			stop = i;
			position = i;
			break;
		}
	}
	for(unsigned int i = position; i < coverage.size(); i++) {
		if(coverage[i] > config.COV_THRES) {
			if(stop + 32000 > i) {
				stop = i;
			}
			else {
				starts.push_back(start-config.COV_PADDING);
				stops.push_back(stop + config.COV_PADDING);
				start = i;
				stop = i;
			}
		}
	}
	starts.push_back(start);
	stops.push_back(stop);

	for(unsigned int i = 0; i < starts.size(); i++) {
		int pos3 = starts[i]-32000;
		int pos4 = stops[i]+1;
		out2 << "> " << starts[i] << "-" << stops[i] << endl;
		for(int j = 0; j < config.COV_GAPLEN; j++) {
			positions.push_back(pos3);
			part_genome += "N";
			out2 << "N";
			pos3++;
		}
		for(long j =  starts[i]; j <= stops[i]; j++) {
			positions.push_back(j);
			out2 << whole_genome[j];
			part_genome.push_back(whole_genome[j]);
		}
		for(int j = 0; j < config.COV_GAPLEN; j++) {
			positions.push_back(pos4);
			out2 << "N";
			part_genome.push_back('N');
			pos4++;
		}
		out2 << endl;
	}

}

/*
 *
 * COVERAGE
 *
 */


void calculateCoverageFromCig(vector<int> &coverage, string infile) {
	ifstream res(infile.c_str());
	string line;
	string linestr;
	while(getline(res, line)) {
		getline(res, line);
		getline(res, line);
		getline(res, line);
		getline(res, line);
		vector<long> array;
		gapArrayFromString(array, linestr);
		for(unsigned int k = 0; k < array.size(); k+=2) {
			int start = array[k];
			int stop = array[k+1];
			for(int l = start; l <= stop; l++) {
				coverage[l] = coverage[l]+1;
			}
		}

	getline(res, line);
	}
}


void calculateCoverageFromResults(vector<int> &coverage, string infile) {
	ifstream res(infile.c_str());
	string line;
	string linestr;
	while(getline(res, line)) {
		getline(res, line);
		getline(res, line);
		linestr = line.substr(0, line.size()-2);
		vector<long> array;
		gapArrayFromString(array, linestr);
		for(unsigned int k = 0; k < array.size(); k+=2) {
			int start = array[k];
			int stop = array[k+1];
			for(int l = start; l <= stop; l++) {
				coverage[l] = coverage[l]+1;
			}
		}

	getline(res, line);
	}
}

void calculateCoverageFromVector(Config &config, vector<int> &coverage, vector<Result> &results) {
	for(unsigned int i = 0; i < results.size(); i++) {
		Result r = results[i];
		if(r.gapArray.size() == 0) {
			r.gapArray.push_back(r.start);
			r.gapArray.push_back(r.stop);
		}
		/*if(r.br % 50 == 0) {
			cout << "read: " << r.br << endl;
		}*/
		for(unsigned int k = 0; k < r.gapArray.size(); k+=2) {
			int start;
			int stop;
			if(config.IS_SECOND_PHASE) {
				start = config.POSITIONS_INDEX[r.gapArray[k]];
				stop = config.POSITIONS_INDEX[r.gapArray[k+1]];
			}
			else {
				start = r.gapArray[k];
				stop = r.gapArray[k+1];
			}
			for(int l = start; l <= stop; l++) {
				coverage[l] = coverage[l]+1;
			}
		}
	}
}

void writeCoverageFiltered(string path, vector <int> &coverage) {
	ofstream outs(path.c_str());
	for(unsigned int i = 870000; i < 880000; i++) {
	    if(coverage[i] > -1) outs << i << "\t";
	    if(coverage[i] > -1) outs << coverage[i] << endl;
	}

}

void writeCoverage(string path, vector <int> &coverage) {
	ofstream outs(path.c_str());
	for(unsigned int i = 0; i < coverage.size(); i++) {
	    outs << coverage[i] << "\t";
	}
	outs << endl;
}

void checkCoverage(vector<int> &coverage) {
	for(unsigned int i = 0; i < coverage.size(); i++) {
	    if(coverage[i] > 0) cout << coverage[i] << " ";
	}
	cout << endl;
}

void evaluateCoverage(vector<int> &final_coverage, vector<int> &correct_coverage, string infile) {
	ofstream streamout(infile.c_str());
	//cout << "final coverage: " << final_coverage.size() << endl;
	//cout << "correct coverage: " << correct_coverage.size() << endl;
	int sum = 0;
//	cout << "starting" << endl;
	for(unsigned int i = 0; i < final_coverage.size(); i++) {
		sum += abs(final_coverage[i] - correct_coverage[i]);
		//cout << i;
	}
	cout << "Coverage difference: " << sum << endl;
	streamout << "Coverage difference: " << sum << endl;
}

