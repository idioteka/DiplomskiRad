
#include "headers.h"
#include "config.h"
#include "statistics.h"
#include "helpers.h"
#include "IO.h"

void findSecondPhaseReads(Config &config, vector<Read> &alignedReads, vector<Read> &secondPhaseReads,
		vector<Result> &results, vector<Result> &doneResults) {

	map<int, Result> resultsMap;
	for(unsigned int i = 0; i < results.size(); i++) {
		Result r = results[i];
		std::map<int,Result>::iterator it = resultsMap.begin();
		resultsMap.insert(it, pair<int,Result>(r.br, results[i]));
	}
	int num = 0;

	for(unsigned int i = 0; i < alignedReads.size(); i++) {
		Read r = alignedReads[i];
		map<int, Result>::iterator pos = resultsMap.find(r.br);
		Result r2 = pos->second;
  		if(r2.maxScore * config.PRECISION_CUTOFF > r2.score) {
			secondPhaseReads.push_back(r);
			num++;
		//	TOTAL_BASE_NUM += r.content.size();
		}
  		else {
  			doneResults.push_back(r2);
  		//	TOTAL_BASE_NUM += r.content.size();
  			config.ALIGNED_BASE_NUM += r.content.size();
  		}
	}
}

void analyseSecondPhaseReads(Info &info, vector<Read> &allReads, vector<Read> &unalignedReads,
		vector<Read> &newAlignedReads, vector<Result> &results,
		vector<Result> &secondResults, vector<Result> &product, vector<Result> &product2) {


	map<int, Result> resultsMap;
	for(unsigned int i = 0; i < results.size(); i++) {
		Result r = results[i];
		std::map<int,Result>::iterator it = resultsMap.begin();
		resultsMap.insert(it, pair<int,Result>(r.br, results[i]));
	}

	map<int, Result> secondResultsMap;
	for(unsigned int i = 0; i < secondResults.size(); i++) {
		Result r = secondResults[i];
		std::map<int,Result>::iterator it = secondResultsMap.begin();
		secondResultsMap.insert(it, pair<int,Result>(r.br, secondResults[i]));
	}

	int better = 0;
	int worse = 0;
	int unfonud = 0;
	int same = 0;

	for(unsigned int i = 0; i < newAlignedReads.size(); i++) {
		Read r = newAlignedReads[i];
		map<int, Result>::iterator it2 = secondResultsMap.find(r.br);
		Result secR = it2->second;
		map<int, Result>::iterator it3 = resultsMap.find(r.br);
		Result firR = it3->second;
		if(secR.score > firR.score) {
			product2.push_back(secR);
			better++;
		}
		else if(secR.score < firR.score) {
			product.push_back(firR);
			worse++;
		}
		else {
			product.push_back(firR);
			same++;
		}
	}
	for(unsigned int i = 0; i < unalignedReads.size(); i++) {
		Read r = unalignedReads[i];
		map<int, Result>::iterator it2 = resultsMap.find(r.br);
		Result r2 = it2->second;
		product.push_back(r2);
		unfonud++;
	}

	info.better = better;
	info.worse = worse;
	info.same = same;
	info.unfound = unfonud;

	cout << "Better: " << better << endl;
	cout << "Worse: " << worse << endl;
	cout << "Same: " << same << endl;
	cout << "Unfound: " << unfonud << endl;

}
