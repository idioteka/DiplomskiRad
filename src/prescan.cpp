
#include "headers.h"
#include "config.h"
#include "heap.h"
#include "score.h"
#include "helpers.h"
#include "key_manager.h"

/*
 *
 * PRESCAN FUNCTIONS
 *
 */

void findMaxQscore2(Config &config, vector<long> &starts, vector<long> &stops, vector<int> &offsets,
		vector<Triplet> &triples, vector<long> &values, vector<int> &keys,
		  int prevMaxHits, bool earlyExit, bool perfectOnly, long *sizes, long *sites, vector<int> &temp){

	int numHits = offsets.size();
	Heap heap = Heap();
	heap.initHeap(offsets.size());
	heap.clear();
	for(int i = 0; i < numHits; i++){
		long start = starts[i];
		long a = sites[start];
		long a2 = a-offsets[i];

		Triplet t(i, start, a2);
		values.push_back(a2);

		heap.add(t);
	}

	int maxQuickScore_ = calcMaxQuickScore(config, offsets, keys);

	int topQscore=-999999999;

	int maxHits=0;

	int approxHitsCutoff;
	int indelCutoff;
	if(perfectOnly){
		approxHitsCutoff=numHits;
		indelCutoff=0;
	}else{
		approxHitsCutoff=max(prevMaxHits, min(config.MIN_APPROX_HITS_TO_KEEP, numHits-1)); //Faster, same accuracy
		indelCutoff=config.MAX_INDEL2;
	}

	while(!heap.isEmpty()){
		Triplet t=heap.peek();
		long site=t.site;
		int centerIndex=t.column;

		long maxNearbySite=site;

		int approxHits=0;

		long minsite=site-min(config.MAX_INDEL, indelCutoff);
		long maxsite=site+config.MAX_INDEL2;
		for(int column=0, chances=numHits-approxHitsCutoff; column<numHits && chances>=0; column++){
			long x=values[column];

			if(x>=minsite && x<=maxsite){
				if(x > maxNearbySite) {
					maxNearbySite = x;
				}
				approxHits++;
			}else{chances--;}
		}

	//	cout << "approxHitsCutoff: " << approxHitsCutoff << endl;
	//	cout << "approxHits: " << approxHits << endl;

		if(approxHits>=approxHitsCutoff){

			int qscore=quickScore(config, values, centerIndex, offsets, approxHits, numHits);

			int scoreZ = scoreZ2(config, values, centerIndex, offsets, approxHits, numHits);
			qscore+=scoreZ;

			if(qscore>topQscore){

				maxHits=max(approxHits, maxHits);
				approxHitsCutoff=max(approxHitsCutoff, approxHits-1); //Best setting for pre-scan

				topQscore=qscore;

				if(qscore>=maxQuickScore_){
					if(earlyExit){
						temp.push_back(topQscore);
						temp.push_back(maxHits);
						return;
					}
				}
			}
		}

		while(heap.peek().site==site){ //Remove all identical elements, and add subsequent elements
			Triplet t2=heap.poll();

			long row=t2.row+1;
			int col=t2.column;
			if(row<stops[col]){
				t2.row=row;
				long a=sites[row];
				long a2=a-offsets[col];

				t2.site=a2;
				values[col]=a2;
				heap.add(t2);
			}else if(earlyExit && (perfectOnly || heap.size()<approxHitsCutoff)){
				temp.push_back(topQscore);
				temp.push_back(maxHits);
				return;
			}
			if(heap.isEmpty()){break;}
		}

	}
	temp.push_back(topQscore);
	temp.push_back(maxHits);
	return;
}

void prescanAllBlocks(Config &config, vector<vector<int> > &prescanResults, int bestScores[], vector<int> &keysP, vector<int> &offsetsP,
		vector<int> &keysM, vector<int> &offsetsM, bool allBasesCovered, long *sizes, long *sites){

	int bestqscore = 0;
	int maxHits = 0;
	int minHitsToScore = config.MIN_APPROX_HITS_TO_KEEP;

	int maxQuickScore_=calcMaxQuickScore(config, offsetsP, keysP);

	vector<int> counts;
	vector<int> scores;
	prescanResults.push_back(counts);
	prescanResults.push_back(scores);

	vector<int> keys = keysP;
	vector<int> offsets = offsetsP;

	vector<long> starts;
	vector<long> stops;

	int numHits = getHits(config, keys, std::numeric_limits<int>::max(), starts, stops, sizes);

	if(numHits < minHitsToScore){
		prescanResults[0].push_back(0);
		prescanResults[1].push_back(-9999);
	}else{

		if(numHits < (int) keys.size()){
			vector<vector<long> > r = shrink(starts, stops, offsets, keys, offsets.size());
			if(r.size() != 0){
				starts = r[0];
				stops = r[1];
				for(int i = 0; i < offsets.size(); i++) {
					offsets.push_back(r[2][i]);
				}
			}
		}

		vector<Triplet> triples;
		vector<long> values;

		vector<int> temp;

		findMaxQscore2(config, starts, stops, offsets, triples, values, keys, minHitsToScore, true,
				bestqscore>=maxQuickScore_ && allBasesCovered, sizes, sites, temp);

		prescanResults[0].push_back(temp[1]);
		prescanResults[1].push_back(temp[0]);

		bestqscore = max(temp[0], bestqscore);
		maxHits = max(maxHits, temp[1]);
		if(bestqscore >= maxQuickScore_ && allBasesCovered){

			minHitsToScore=max(minHitsToScore, maxHits);
			bestScores[1]=max(bestScores[1], maxHits);
			bestScores[3]=max(bestScores[3], bestqscore);
			return;
		}
	}
	keys.clear();
	keys = keysM;
	offsets.clear();
	offsets = offsetsM;

	starts.clear();
	stops.clear();

	numHits = getHits(config, keys, std::numeric_limits<int>::max(), starts, stops, sizes);

	if(numHits < minHitsToScore){
		prescanResults[0].push_back(0);
		prescanResults[1].push_back(-9999);
	}else{

		if(numHits < (int) keys.size()){
			vector<vector<long> > r = shrink(starts, stops, offsets, keys, offsets.size());
			if(r.size() != 0){
				starts = r[0];
				stops = r[1];
				for(int i = 0; i < offsets.size(); i++) {
					offsets.push_back(r[2][i]);
				}
			}
		}

		vector<Triplet> triples;
		vector<long> values;

		vector<int> temp;
		findMaxQscore2(config, starts, stops, offsets, triples, values, keys, minHitsToScore, true,
				bestqscore>=maxQuickScore_ && allBasesCovered, sizes, sites, temp);

		prescanResults[0].push_back(temp[1]);
		prescanResults[1].push_back(temp[0]);

		bestqscore = max(temp[0], bestqscore);
		maxHits = max(maxHits, temp[1]);
		if(bestqscore >= maxQuickScore_ && allBasesCovered){

			minHitsToScore=max(minHitsToScore, maxHits);
			bestScores[1]=max(bestScores[1], maxHits);
			bestScores[3]=max(bestScores[3], bestqscore);
			return;
		}
	}


	bestScores[1]=max(bestScores[1], maxHits);
	bestScores[3]=max(bestScores[3], bestqscore);

	return;
}
