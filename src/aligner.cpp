#include "headers.h"
#include "genome_index.h"
#include "heap.h"
#include "const.h"
#include "key_manager.h"
#include "helpers.h"
#include "score.h"
#include "match_string.h"
#include "prescan.h"
#include "coverage.h"
#include "statistics.h"
#include "statistics2.h"
#include "parameters.h"
#include "IO.h"
#include "second_phase.h"
#include "config.h"

bool isNull(SiteScore &ss) {
	if(ss.start == -1 && ss.stop == -1 && ss.score == -1) return true;
	else return false;
}

bool isNull(Triplet &t){
	if(t.column == -1 && t.row == -1 && t.site == -1) return true;
	else return false;
}

int calcApproxHitsCutoff(Config &config, int keys, int max_hits, int current_cutoff, bool perfect) {

	int reduction = min( max( (max_hits)/config.HIT_REDUCTION_DIV, config.MAX_HITS_REDUCTION2), max(config.MAXIMUM_MAX_HITS_REDUCTION, keys/8));
	int approxCuttof = max(max(config.MIN_APPROX_HITS_TO_KEEP, current_cutoff), max_hits-reduction);
	return approxCuttof;
}

void makeGapArray(vector<long> &locArray, long minLoc, long minGap, vector<long> &gapArray){
	int gaps = 0;
	bool doSort=false;

	if(locArray[0] < 0) { locArray[0] = minLoc; }

	for(unsigned int i = 1; i < locArray.size(); i++){
		if(locArray[i] < 0){ locArray[i] = locArray[i-1] + 1; }
		else { locArray[i] += i; }
		if(locArray[i] < locArray[i-1]) {
			doSort = true;
		}
	}

	long maxsite = 0;
	for(unsigned int i = locArray.size()-1; i > 0; i--) {
		if(locArray[i] > -1) {
			maxsite = locArray[i];
			break;
		}
	}

	long minsite = 0;
	for(unsigned int i = 0; i < locArray.size(); i++) {
		if(locArray[i] > -1) {
			minsite = locArray[i];
			break;
		}
	}

	if(doSort) {
		sort(locArray.begin(), locArray.begin() + locArray.size());
	}

	for(unsigned int i = 1; i < locArray.size(); i++){
		long dif = locArray[i] - locArray[i-1];
		if(dif>minGap){
			gaps++;
		}
	}

	if(gaps < 1) {
		locArray.clear();
		return;
	}

	gapArray.push_back(minsite);

	for(unsigned int i = 1, j = 1; i < locArray.size(); i++){
		if(locArray[i-1] < minsite) continue;
		long dif = locArray[i] - locArray[i-1];
		if(dif > minGap) {
			if(locArray[i-1] > maxsite) {
				break;
			}
			else if(locArray[i] > maxsite && locArray[i-1] < maxsite) {
				gapArray.push_back(locArray[i-1]);
				break;
			}
			else if(locArray[i-1] == maxsite) {
				break;
			}
			gapArray.push_back(locArray[i-1]);
			gapArray.push_back(locArray[i]);
			j+=2;
		}
	}

	gapArray.push_back(maxsite);
}

/*
 *
 * ALIGNER
 *
 */

void printLocationArray(vector<int> locArray) {
	for(unsigned int i = 0; i < locArray.size(); i++) {
		cout << locArray[i] << " ";
	}
	cout << endl;
}

int align(Config &config, int bestScores[], vector<int> &keys, string &read, vector<int> &offsets, long *sizes,
		long *sites, vector<SiteScore> &results, bool all_bases_covered,
		int max_score, int max_quick_score, bool fully_defined, int strand, string &whole_genome) {
	vector<long> starts;
	vector<long> stops;
	int get_hits = getHits(config, keys, config.MAX_LEN, starts, stops, sizes);

	bool filter_by_qscore = (get_hits >= 5);
	int minScore = (int) (config.MIN_SCORE_MULT * max_score);
	int minQuickScore = (int)(config.MIN_QSCORE_MULT * max_quick_score);

	int current_top_score = bestScores[0];
	int bestqscore = bestScores[3];
	int max_hits = bestScores[1];
	int perfects_found = bestScores[5];

	int cutoff = max(minScore, (int)(current_top_score * config.DYNAMIC_SCORE_THRESH));
	int qcutoff = max(bestScores[2], minQuickScore);

	int approx_hits_cutoff = calcApproxHitsCutoff(config, get_hits, max_hits, config.MIN_APPROX_HITS_TO_KEEP, current_top_score >= max_score);
	if(approx_hits_cutoff > get_hits){
		return current_top_score;
	}
	bool short_circuit = (all_bases_covered && get_hits == (int)offsets.size() && filter_by_qscore);

	if(current_top_score >= max_score){
		qcutoff = max(qcutoff, (int)(max_quick_score*config.DYNAMIC_QSCORE_THRESH_PERFECT));
	}

	Heap heap = Heap();

	vector<long> sites_tmp;
	heap.clear();
	heap.initHeap(offsets.size());

	for(int i = 0; i < get_hits; i++){
		long start = starts[i];

		long a = sites[start];
		long a2 = a - offsets[i];

		Triplet t(i, start, a2);
		sites_tmp.push_back(a2);

		heap.add(t);
	}

	vector<long> loc_array;
	SiteScore prevSS;

	int numberOfSites = 0;

	while(!heap.isEmpty()) {
		Triplet t = heap.peek();
		long site = t.site;
		int center_index = t.column;
		long max_nearby_site = site;
		int approx_hits = 0;
		long minsite = site - config.MAX_INDEL;
		long maxsite = site + config.MAX_INDEL2;
		for(int column = 0, chances = get_hits - approx_hits_cutoff; column < get_hits && chances >= 0; column++) {
			long x = sites_tmp[column];
			if(x >= minsite && x <= maxsite) {
				if(x > max_nearby_site) {
					max_nearby_site = x;
				}
				approx_hits++;
			} else {
				chances--;
			}
		}

		//cout << "approx hits: " << approx_hits << endl;
		//cout << "approx_hits_cutoff: " << approx_hits_cutoff << endl;

		if(approx_hits >= approx_hits_cutoff){

			int score;
			int qscore = quickScore(config, sites_tmp, center_index, offsets, approx_hits, get_hits);
			int scoreZ = scoreZ2(config, sites_tmp, center_index, offsets, approx_hits, get_hits);
			qscore+=scoreZ;

			long mapStart = site;
			long mapStop = max_nearby_site;

			if(qscore < qcutoff){
				score=-1;
			}else{
				if(short_circuit && qscore == max_quick_score){
					score = max_score;
				}else{
					score = extendScore(config, read, offsets, sites_tmp, center_index, loc_array, get_hits, approx_hits, whole_genome);

					long min = std::numeric_limits<long>::max();
					long max = std::numeric_limits<long>::min();
					for(unsigned int i = 0; i < loc_array.size(); i++){
						long x = loc_array[i];
						if(x >- 1){
							if(x < min) {
								min = x;
							}
							if(x > max) {
								max = x;
							}
						}
					}
					if(min<0 || max<0){
						score=-99999;
					}
					mapStart = min;
					mapStop = max;
				}
				if(score == max_score){
					qcutoff = max(qcutoff, (int)(max_quick_score * config.DYNAMIC_QSCORE_THRESH_PERFECT));
					approx_hits_cutoff = calcApproxHitsCutoff(config, get_hits, max_hits, config.MIN_APPROX_HITS_TO_KEEP, current_top_score >= max_score);
				}

				if(score >= cutoff){
					qcutoff = max(qcutoff, (int)(qscore * config.DYNAMIC_QSCORE_THRESH));
					bestqscore = max(qscore, bestqscore);
				}
			}
			if(score>=cutoff){
				//cout << "cutoff " << cutoff << "\n";

				if(score > current_top_score){
					max_hits = max(approx_hits, max_hits);
					approx_hits_cutoff = calcApproxHitsCutoff(config, get_hits, max_hits, approx_hits_cutoff, current_top_score >= max_score);
					cutoff = max(cutoff, (int)(score * config.DYNAMIC_SCORE_THRESH));

					if(score >= max_score){
						cutoff = max(cutoff, (int)(score*0.95f));
					}
					current_top_score = score;
				}

				long site2 = mapStart;
				long site3 = mapStop + read.size() - 1;

				vector<long> gapArray;

				//printLocationArray(loc_array);
				if(site3-site2 >= config.MINGAP + (int)read.size()){

					makeGapArray(loc_array, site2, config.MINGAP, gapArray);
					//cout << "score:" << score << endl;
					/*if(gapArray.size() != 0){
						gapArray[0] = min(gapArray[0], site2);
						gapArray[gapArray.size()-1] = max(gapArray[gapArray.size()-1], site3);
					}*/
					if(gapArray.size() != 0) {
						site2 = gapArray[0];
						site3 = gapArray[gapArray.size()-1];
					}
				}

				SiteScore ss;
				bool perfect1 = score == max_score && fully_defined;
				bool inbounds = (site2 >= 0 && site3 < (long)whole_genome.size());

				if(inbounds && (int) gapArray.size() == 0 && isNull(prevSS) == false
						&& strand == prevSS.strand && overlap(prevSS.start, prevSS.stop, site2, site3)){

					int betterScore = max(score, prevSS.score);
					int minStart = min(prevSS.start, site2);
					int maxStop = max(prevSS.stop, site3);
					bool perfect2 = prevSS.score == max_score && fully_defined;

					bool shortEnough = (maxStop-minStart < 2 * (int)read.size());

					if(prevSS.start == site2 && prevSS.stop == site3){
						results[results.size()-1].score = betterScore;
						results[results.size()-1].perfect = (prevSS.perfect || perfect1 || perfect2);
						//if(prevSS.perfect){prevSS.semiperfect=true;}

					} else if(shortEnough && prevSS.start == site2){
						if(perfect2){
							//do nothing
						}else if(perfect1){
							results[results.size()-1].stop = site3;
							if(!prevSS.perfect) { perfects_found++; }
							results[results.size()-1].perfect = true;
							results[results.size()-1].score = betterScore;
						}else{
							results[results.size()-1].stop = maxStop;
							// prevSS.setPerfect(bases); TODO add perfect set
						}
						results[results.size()-1].score = betterScore;
					}else if(shortEnough && prevSS.stop == site3){
						if(perfect2){
							//do nothing
						}else if(perfect1){
							results[results.size()-1].start = site2;
							if(!prevSS.perfect) { perfects_found++; }
							results[results.size()-1].perfect = true;
							results[results.size()-1].score = betterScore;
						}else{
							results[results.size()-1].start = minStart;
							//prevSS.setPerfect(bases); TODO add perfect set
						}
						results[results.size()-1].score = betterScore;
					}else if(shortEnough && (maxStop - minStart <= (int) read.size() + config.MAX_SUBSUMPTION_LENGTH) && !perfect1 && !perfect2){
						results[results.size()-1].start = minStart;
						results[results.size()-1].stop = maxStop;
						results[results.size()-1].score = betterScore;
						//prevSS.setPerfect(bases); TODO add perfect set
					}else{
						ss = SiteScore(site2, site3, score, approx_hits, perfect1, strand, gapArray);
						//printLocationArray(loc_array);
						//if(!perfect1){ss.setPerfect(bases);} TODO add perfect set
					}

				}else if(inbounds){
					ss = SiteScore(site2, site3, score, approx_hits, perfect1, strand, gapArray);
					//printLocationArray(loc_array);
					//if(!perfect1){ss.setPerfect(bases);} TODO add perfect setter
				}

				if(isNull(ss) == false){
					results.push_back(ss);
					numberOfSites++;
					if(ss.perfect){
						if(isNull(prevSS) == true || !prevSS.perfect || overlap(ss.start, ss.stop, prevSS.start, prevSS.stop)){
							perfects_found++;
							if(perfects_found >= 2) {
								heap.clear();
								break;
							}
						}
					}
					prevSS = ss;
				}
			}
		}

		while(heap.peek().site == site){

			Triplet t2 = heap.poll();
			long row = t2.row+1;
			long col = t2.column;
			if(row < stops[col]){
				t2.row=row;
				long a = sites[row];
				long a2 = a-offsets[col];
				t2.site = a2;
				sites_tmp[col]=a2;
				heap.add(t2);
			}
			else if(heap.size() < approx_hits_cutoff){
				return current_top_score;
			}
		}
	}
}

void processRead(Config &config, long *sizes, long *sites, string &r1, vector<vector<Result> > &databaseResults, vector<Result> &resultsFinal, vector<Read> &aligned_reads, vector<Read> &unaligned_reads, string &whole_genome, int threadId, int br) {

	int split_size = r1.size() / config.SPLIT_COUNT;
	for(int i = 0; i < config.SPLIT_COUNT; i++) {
		string read = r1.substr(i * split_size, split_size);
		string read_reverse;
		reverseComplementRead(read_reverse, read);

		vector<int> offsets;
		makeOffsets(config, read, offsets);
		vector<int> read_keys;
		getReadKeys(config, read, offsets, read_keys);
		vector<int> read_keys_copy = read_keys;

		int num_hits = countHits(config, read_keys_copy, sizes, config.MAX_LEN);
		num_hits = trimHits(config, read_keys, read_keys_copy, sizes, num_hits);

		vector<vector<int> > res;

		//TODO add this if if(num_hits < read_keys.size()) {
		res = shrinkArrays(offsets, read_keys_copy, num_hits);

		offsets = res[0];
		vector<int> read_keys_final = res[1];

		vector<SiteScore> results;

		bool all_bases_covered = checkIfAllBasesCovered(config, offsets, read);
		bool pretend_all_bases_covered = calculatePretendAllBasesCovered(config, read_keys, read_keys_final, offsets, all_bases_covered, read);
		bool prescan_qscore = (num_hits >= 0);

		int quick_max_score = calcMaxQuickScore(config, offsets, read_keys_final);

		int bestScores[6] = {0};
		vector<int> precounts;
		vector<int> prescores;

		int hitsCutoff=0;
		int qscoreCutoff = (int)(config.MIN_QSCORE_MULT*quick_max_score);

		vector<int> offsets_reversed;
		reverseOffsets(offsets_reversed, offsets, config.KEYLEN, read.size());

		vector<int> keys_reversed;
		reverseComplementKeys(keys_reversed, read_keys_final, config.KEYLEN);

		if(prescan_qscore){
			vector<vector<int> > prescanResults;
			prescanAllBlocks(config, prescanResults, bestScores, read_keys_final, offsets, keys_reversed, offsets_reversed, pretend_all_bases_covered, sizes, sites);

			precounts=prescanResults[0];
			prescores=prescanResults[1];

			if(bestScores[1]<config.MIN_APPROX_HITS_TO_KEEP) {
				Read rea(br, read);
				unaligned_reads.push_back(rea);
				return;
			}
			if(bestScores[3]<quick_max_score*config.MIN_QSCORE_MULT2) {
				Read rea(br, read);
				unaligned_reads.push_back(rea);
				return;
			} //if(bestScores[3]<maxQuickScore(offsetsP, keyScoresP)*.10f){return result;}

			if(bestScores[3]>=quick_max_score && pretend_all_bases_covered){

				hitsCutoff=calcApproxHitsCutoff(config, read_keys_final.size(), bestScores[1], config.MIN_APPROX_HITS_TO_KEEP, true);
				qscoreCutoff=max(qscoreCutoff, (int)(bestScores[3]*config.DYNAMIC_QSCORE_THRESH_PERFECT));

			}else{
				hitsCutoff=calcApproxHitsCutoff(config, read_keys_final.size(), bestScores[1], config.MIN_APPROX_HITS_TO_KEEP, false);
				qscoreCutoff=max(qscoreCutoff, (int)(bestScores[3]*config.PRESCAN_QSCORE_THRESH));
			}
		}

		int max_score = calcMaxScore(config, read);

		bool fully_defined = isFullyDefined(read);

		int score = -1;
		if(precounts[0] >= hitsCutoff && prescores[0] >= qscoreCutoff) {
			score = align(config, bestScores, read_keys_final, read, offsets, sizes, sites, results, all_bases_covered, max_score, quick_max_score, fully_defined, 0, whole_genome);
		}
		if(precounts[1] >= hitsCutoff && prescores[1] >= qscoreCutoff) {
			score = align(config, bestScores, keys_reversed, read_reverse, offsets_reversed, sizes, sites, results, all_bases_covered, max_score, quick_max_score, fully_defined, 1, whole_genome);
		}

		//cout << "Number of results: " << results.size() << endl;

		if(results.size() != 0) {
			filterBadReads(config, score, results, read, read_reverse, sizes, sites, databaseResults, resultsFinal, whole_genome, threadId, br, max_score);
			config.ALIGNED_BASE_NUM += read.size();
			Read rea(br, read);
			aligned_reads.push_back(rea);
		}
		else {
			Read rea(br, read);
			unaligned_reads.push_back(rea);
		}

	}
	return;
}

void *preProcessRead(void *threadid) {
	ThreadData3 *td = (ThreadData3 *) threadid;
	for(int i = td->start; i < td->stop; i++) {

		Read read = (*(td->reads))[i];

		if(false) {
			if((i+1)%2 == 0) cout << "Read: " << read.br << " started." << endl;
		}
		processRead(*(td->config), td->sizes, td->sites, read.content, *(td->databaseResults), *(td->results), *(td->aligned_reads), *(td->unaligned_reads), *(td->whole_genome),  td->thread_id, read.br);
	}
	pthread_exit(NULL);
}

void executeMethod(Statistic global_stat, Info &info, Config &config, long **res,  map<int, FastaRead> &read_names, map<int, Result> &correct_results,
		string &whole_genome, vector<Read> &reads, bool read_index, vector<Read> &unaligned,
		vector<Read> &aligned, vector<Result> &tmp_results, vector<vector<Result> > &database, vector<ReferenceSegment> &referenceSegments) {
	long *sizes = res[1];
	long *sites = res[3];

	config.LENGTH_OF_SIZES = res[0][0];
	config.LENGTH_OF_SITES = res[2][0];

	cout << "Number of reads: "<< reads.size() << "\n";
	string build = SSTR(config.BUILD_NUMBER);
	pthread_t threads[config.THREAD_NUM];

	int rc;
	pthread_attr_t attr;
	void *status;
	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

	vector<vector<Result> > set_of_results;
	for(int i = 0; i < config.THREAD_NUM; i++) {
		vector<Result> tmp;
		set_of_results.push_back(tmp);
	}
	vector<vector<Read> > unaligned_reads;
	for(int i = 0; i < config.THREAD_NUM; i++) {
		vector<Read> tmp;
		unaligned_reads.push_back(tmp);
	}
	vector<vector<Read> > aligned_reads;
	for(int i = 0; i < config.THREAD_NUM; i++) {
		vector<Read> tmp;
		aligned_reads.push_back(tmp);
	}
	ThreadData3 datas3[config.THREAD_NUM];
	int difference = reads.size()/config.THREAD_NUM;

	vector<vector<vector<Result> > > databaseResults;
	for(int i = 0; i < config.THREAD_NUM; i++) {
		vector<vector<Result> > tmp;
		databaseResults.push_back(tmp);
	}

	for(int i = 0; i < config.THREAD_NUM; i++) {
		int start =  i*difference;
		int stop =  i*difference + difference;
		datas3[i] = ThreadData3(i, sizes, sites, &reads, start, stop, &config, &whole_genome, &set_of_results[i], &aligned_reads[i], &unaligned_reads[i], &databaseResults[i]);
	}
	datas3[config.THREAD_NUM-1].stop = reads.size();
	for(int i = 0; i < config.THREAD_NUM; i++) {
		cout << "Thread " << i << " started." << endl;
		rc = pthread_create(&threads[i], NULL, preProcessRead, (void *) &datas3[i]);
	}

	timeval t1, t2;
	gettimeofday(&t1, NULL);
	long startday = t1.tv_sec;
	long startday2 = t1.tv_usec;

	pthread_attr_destroy(&attr);
	for(int i=0; i < config.THREAD_NUM; i++ ){
		rc = pthread_join(threads[i], &status);
		if (rc){
			cout << "Error:unable to join," << rc << endl;
			exit(-1);
		}
		cout << "Main: completed thread id :" << i ;
		cout << "  exiting with status :" << status << endl;
	}

	gettimeofday(&t2, NULL);
	long endday = t2.tv_sec;
	long endday2 = t2.tv_usec;

	double timefinal = ((endday - startday) * 1000000.0 + (endday2 - startday2))/ 1000000;
	cout << "Reads aligned: " << timefinal << " seconds." << endl;
	if(!config.IS_SECOND_PHASE) {
		info.first_align_time2 = timefinal;
	}
	else if(!config.IS_THIRD_PHASE) {
		info.second_align_time2 = timefinal;
	}
	else {
		info.second_align_time2 = timefinal;
	}

	for(unsigned int i = 0; i < set_of_results.size(); i++) {
		for(unsigned int j = 0; j < set_of_results[i].size(); j++) {
			tmp_results.push_back(set_of_results[i][j]);
		}
	}
	set_of_results.clear();

	for(unsigned int i = 0; i < unaligned_reads.size(); i++) {
		for(unsigned int j = 0; j < unaligned_reads[i].size(); j++) {
			unaligned.push_back(unaligned_reads[i][j]);
		}
	}
	unaligned_reads.clear();

	writeReads(unaligned, config.OUTDIR + "//" + "unaligned.txt");

	for(unsigned int i = 0; i < aligned_reads.size(); i++) {
		for(unsigned int j = 0; j < aligned_reads[i].size(); j++) {
			aligned.push_back(aligned_reads[i][j]);
		}
	}
	aligned_reads.clear();

	for(unsigned int i = 0; i < databaseResults.size(); i++) {
		for(unsigned int j = 0; j < databaseResults[i].size(); j++) {
			database.push_back(databaseResults[i][j]);
		}
	}
	databaseResults.clear();

	//writeReads(aligned, config.OUTDIR + "//" + "aligned.txt");

	cout << "Aligned reads: " << tmp_results.size() << endl;

	string addon = "";
	if(config.IS_SECOND_PHASE) {
		addon = "x";
	}

	gettimeofday(&t1, NULL);
	startday = t1.tv_sec;
	startday2 = t1.tv_usec;
	writeSamResults(tmp_results, reads, read_names, config.OUTDIR + "//" + "results" + addon + build + ".sam", referenceSegments);
	writeResults(config, tmp_results, correct_results, config.OUTDIR + "//" + "results" + addon + build + ".txt");

	if(!config.IS_SECOND_PHASE) {
	//	writeResults(config, tmp_results, correct_results, config.OUTDIR + "//" + "results" + addon + build + ".txt");
		//writeSamResults(tmp_results, reads, read_names, config.OUTDIR + "//" + "results" + addon + build + ".sam");
	}
	else if(!config.SECOND_PHASE_MODE) {
	//	writeResults(config, tmp_results, correct_results, config.OUTDIR + "//" + "results" + addon + build + ".txt");
		//writeSamResults(tmp_results, reads, read_names, config.OUTDIR + "//" + "results" + addon + build + ".sam");
	}
	gettimeofday(&t2, NULL);
	endday = t2.tv_sec;
	endday2 = t2.tv_usec;
	timefinal = ((endday - startday) * 1000000.0 + (endday2 - startday2))/ 1000000;
	cout << "Results writen: " << timefinal << " seconds." << endl;

	gettimeofday(&t1, NULL);
	startday = t1.tv_sec;
	startday2 = t1.tv_usec;
	Statistic s = Statistic(reads.size());

	writeStatistics2(config, tmp_results, correct_results, config.OUTDIR + "//" + "statistics" + addon + build + ".txt", config.OUTDIR + "//" + "sss" + addon + build + ".txt", s);

	gettimeofday(&t2, NULL);
	endday = t2.tv_sec;
	endday2 = t2.tv_usec;
	timefinal = ((endday - startday) * 1000000.0 + (endday2 - startday2))/ 1000000;
	cout << "Statistics written: " << timefinal << " seconds. " << endl;

	int memory = getValue();
	if(!config.IS_SECOND_PHASE) {
		info.memory_after_first_align = memory;
	}
	else if(!config.IS_THIRD_PHASE) {
		info.memory_after_second_align = memory;
	}
	else {
		info.memory_after_third_align = memory;
	}

	if(read_index) {
		free(sizes);
		free(sites);
	}
	else {
		delete [] sizes;
		delete [] sites;
	}

	delete [] res[0];
	delete [] res[2];
	delete [] res;

}

int main(int argc, char *argv[]) {

	Config config;
	Info info;
	Statistic global_stat;

	timeval t3, t4;
	gettimeofday(&t3, NULL);
	long startday3 = t3.tv_sec;
	long startday4 = t3.tv_usec;

	if(readParams(argc, argv, config, info) < 0) exit(-1);

	updatePrecision2(config, config.MULTY_PRECISION);

	string genome_ref = argv[3];
	string whole_genome;
	bool read_index = createDirectories(config);
	vector<Read> tmp_unaligned_reads;
	vector<Read> tmp_aligned_reads;
	vector<ReferenceSegment> referenceSegments;

	timeval t1, t2;
	gettimeofday(&t1, NULL);
	long startday = t1.tv_sec;
	long startday2 = t1.tv_usec;

	long **res;
	if(read_index) {
		cout << "Reading index..." << endl;
		res = readIndex(whole_genome, genome_ref, config.INDEX_LOCATION, false, config, referenceSegments);
	}
	else {
		cout << "Creating index..." << endl;
		res = createIndex(true, whole_genome, false, genome_ref, config.INDEX_LOCATION, config, referenceSegments);
	}

	gettimeofday(&t2, NULL);
	long endday = t2.tv_sec;
	long endday2 = t2.tv_usec;
	double timefinal = ((endday - startday) * 1000000.0 + (endday2 - startday2))/ 1000000;

	info.index_generation_time = timefinal;
	cout << "genome size: " << whole_genome.size() << endl;
	info.genome_size = whole_genome.size();

	vector<Read> reads;
	map<int, FastaRead> read_names;
	map<int, Result> correct_results;
	readReads(config, reads, read_names, correct_results, argv[2]);
	cout << "Read names: " << read_names.size() << endl;
	vector<Result> results;

	global_stat.total_reads = reads.size() * config.SPLIT_COUNT;
	global_stat.total_bases = config.TOTAL_BASE_NUM;

	info.read_number = reads.size();
	info.first_phase_readlen = reads[0].content.size() / config.SPLIT_COUNT;

	gettimeofday(&t1, NULL);
	startday = t1.tv_sec;
	startday2 = t1.tv_usec;

	vector<vector<Result> > database;

	executeMethod(global_stat, info, config, res, read_names, correct_results, whole_genome, reads, read_index, tmp_unaligned_reads, tmp_aligned_reads, results, database, referenceSegments);

	gettimeofday(&t2, NULL);
	endday = t2.tv_sec;
	endday2 = t2.tv_usec;
	timefinal = ((endday - startday) * 1000000.0 + (endday2 - startday2))/ 1000000;
	info.first_phase_time = timefinal;

	gettimeofday(&t4, NULL);
	long endday3 = t4.tv_sec;
	long endday4 = t4.tv_usec;
	double timefinal2 = ((endday3 - startday3) * 1000000.0 + (endday4 - startday4))/ 1000000;


	info.total_execution_time = timefinal2;
	writeInfo(info, config);

	cout << "Program ended." << endl;

	return 0;
}

