
#include "config.h"

void findSecondPhaseReads(Config &config, vector<Read> &alignedReads, vector<Read> &secondPhaseReads,
		vector<Result> &results, vector<Result> &doneResults);

void analyseSecondPhaseReads(Info &info, vector<Read> &allReads, vector<Read> &unalignedReads,
		vector<Read> &newAlignedReads, vector<Result> &results,
		vector<Result> &secondResults, vector<Result> &product, vector<Result> &product2);
