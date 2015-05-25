
#include "score.h"
#include "config.h"

/*
 *
 * MAKE MATCH STRING
 *
 */

int scoreNoIndelsAndMakeMatchString(string &read, string &ref, long refStart, string &match);
long makeGref(string &ref, vector<long> &gaps, long refStartLoc, long refEndLoc, string &gref);
void fillUnlimited(string &read, string &ref, long refStartLoc, long refEndLoc,
		long max[], vector<vector<vector<long> > > &packed);
void traceback(string &read, string &ref, long refStartLoc, long refEndLoc, int row,
		int col, int state, Result &r, vector<vector<vector<long> > > &packed);
void fillLimited(string &read, string &ref, long refStartLoc, long refEndLoc, int score,
		vector<long> &gaps, long max[], Result &r);
void makeMatchStringForSite(Config &config, SiteScore ss, string &read, long *sizes,
		long *sites, Result &r, string &whole_genome, int threadId);
void makeMatchString(vector<SiteScore> &results, string &read, string &read_reverse,
		long *sizes, long *sites, vector<Result> &resultsFinal, string &whole_genome,
		int threadId, int br, int maxScore);
void filterBadReads(Config &config, int score, vector<SiteScore> &results, string &read, string &read_reverse,
		long *sizes, long *sites, vector<vector<Result> > &databaseResults, vector<Result> &resultsFinal, string &whole_genome,
		int threadId, int br, int maxScore);
