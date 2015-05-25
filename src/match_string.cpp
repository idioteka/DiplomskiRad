
#include "headers.h"
#include "config.h"
#include "score.h"
#include "const.h"

/*
 *
 * MAKE MATCH STRING
 *
 */

int scoreNoIndelsAndMakeMatchString(Config &config, string &read, string &ref, long refStart, string &match){
	int score = 0;
	int mode = -1;
	int timeInMode = 0;

	long readStart = 0;
	long readStop = read.size();
	long refStop = refStart + read.size();
	if(refStart < 0 || refStop > (long) ref.size()) {return -99999;}
	if(refStart < 0){
		readStart = 0 - refStart;
		score += config.POINTS_NOREF * readStart;
	}
	if(refStop > (long) ref.size()){
		long dif = (refStop - ref.size());
		readStop -= dif;
		score += config.POINTS_NOREF * dif;
	}

	for(long i = readStart; i < readStop; i++){
		char c = read[i];
		char r = ref[refStart+i];

		if(c == r && c != 'N'){
			if(mode == config.MODE_MS){
				timeInMode++;
				score += config.POINTS_MATCH2;
			}else{
				timeInMode = 0;
				score += config.POINTS_MATCH;
			}
			match.push_back('m');
			mode = config.MODE_MS;
		}else if(c<0 || c=='N'){
			score+=config.POINTS_NOCALL;
			match.push_back('N');
		}else if(r<0 || r=='N'){
			score+=config.POINTS_NOREF;
			match.push_back('N');
		}else{
			match.push_back('s');
			if(mode == config.MODE_SUB) {timeInMode++;}
			else{timeInMode=0;}

			score+=(POINTS_SUB_ARRAY[timeInMode+1]);

			mode=config.MODE_SUB;
		}
	}
	return score;
}

long makeGref(Config &config, string &ref, vector<long> &gaps, long refStartLoc, long refEndLoc, string &gref) {

	long g0_old = gaps[0];
	long gN_old = gaps[gaps.size()-1];
	gaps[0] = min(gaps[0], refStartLoc);
	gaps[gaps.size()-1] = max(gN_old, refEndLoc);

	long gpos = 0;

	for(unsigned int i=0; i<gaps.size(); i+=2){
		long x = gaps[i];
		long y = gaps[i+1];

		for(long r=x; r<=y; r++, gpos++){
			gref.push_back(ref[r]);
		}

		if(i+2 < gaps.size()){
			long z = gaps[i+2];
			long gap = z-y-1;

			long rem = gap%config.GAPLEN;
			long lim = y+ config.GAPBUFFER +rem;

			long div = (gap-config.GAPBUFFER2)/config.GAPLEN;

			for(long r=y+1; r<=lim; r++, gpos++){
				gref.push_back(ref[r]);
			}
			for(long g=0; g<div; g++, gpos++){
				gref.push_back(config.GAPC);
			}
			for(long r=z-config.GAPBUFFER; r<z; r++, gpos++){
				gref.push_back(ref[r]);
			}
		}
	}
	long greflimit = gpos;

	gaps[0] = g0_old;
	gaps[gaps.size()-1]=gN_old;

	//cout << gref << endl;

	return greflimit;
}

void fillUnlimited(Config &config, string &read, string &ref, long refStartLoc, long refEndLoc, long max[], vector<vector<vector<long> > > &packed) {
	int rows = read.size();
	long columns = refEndLoc-refStartLoc+1;

	long maxGain = (read.size()-1) * config.POINTSoff_MATCH2+config.POINTSoff_MATCH;
	long subfloor = 0-2*maxGain;
	long BARRIER_I2 = rows-config.BARRIER_I1;
	long BARRIER_I2b = columns-1;
	long BARRIER_D2 = rows-config.BARRIER_D1;

	for(int row=1; row<=rows; row++){
		for(int col=1; col<=columns; col++){
			char call0=(row<2 ? '?' : toupper(read[row-2]));
			char call1=toupper(read[row-1]);
			char ref0=(col<2 ? '!' : toupper(ref[refStartLoc+col-2]));
			char ref1=toupper(ref[refStartLoc+col-1]);

			bool match=(call1==ref1 && ref1!='N');
			bool prevMatch=(call0==ref0 && ref0!='N');

			bool gap=(ref1==config.GAPC);

			if(gap){
				packed[config.MODE_MS][row][col]=subfloor;
			}else{//Calculate match and sub scores

				long scoreFromDiag=packed[config.MODE_MS][row-1][col-1]&config.SCOREMASK;
				long scoreFromDel=packed[config.MODE_DEL][row-1][col-1]&config.SCOREMASK;
				long scoreFromIns=packed[config.MODE_INS][row-1][col-1]&config.SCOREMASK;
				long streak=(packed[config.MODE_MS][row-1][col-1]&config.TIMEMASK);

				{//Calculate match/sub score

					if(match){

						long scoreMS=scoreFromDiag+(prevMatch ? config.POINTSoff_MATCH2 : config.POINTSoff_MATCH);

						long scoreD=scoreFromDel+config.POINTSoff_MATCH;
						long scoreI=scoreFromIns+config.POINTSoff_MATCH;

						long score;
						long time;
						if(scoreMS>=scoreD && scoreMS>=scoreI){
							score=scoreMS;
							time=(prevMatch ? streak+1 : 1);
						}else if(scoreD>=scoreI){
							score=scoreD;
							time=1;
						}else{
							score=scoreI;
							time=1;
						}

						if(time>config.MAX_TIME){time=config.MAX_TIME-config.MASK5;}
						packed[config.MODE_MS][row][col]=(score|time);

					}else{

						long scoreMS;
						if(ref1!='N' && call1!='N'){
							scoreMS=scoreFromDiag+(prevMatch ? (streak<=1 ? config.POINTSoff_SUBR : config.POINTSoff_SUB) :
								POINTSoff_SUB_ARRAY[streak+1]);
						}else{
							scoreMS=scoreFromDiag+config.POINTSoff_NOCALL;
						}

						long scoreD=scoreFromDel+config.POINTSoff_SUB; //+2 to move it as close as possible to the deletion / insertion
						long scoreI=scoreFromIns+config.POINTSoff_SUB;

						long score;
						long time;
						if(scoreMS>=scoreD && scoreMS>=scoreI){
							score=scoreMS;
							time=(prevMatch ? 1 : streak+1);
						}else if(scoreD>=scoreI){
							score=scoreD;
							time=1;
						}else{
							score=scoreI;
							time=1;
						}

						if(time>config.MAX_TIME){time=config.MAX_TIME-config.MASK5;}

						packed[config.MODE_MS][row][col]=(score|time);

					}
				}
			}

			if(row<config.BARRIER_D1 || row>BARRIER_D2){
				packed[config.MODE_DEL][row][col]=subfloor;
			}else{//Calculate DEL score

				long streak=packed[config.MODE_DEL][row][col-1]&config.TIMEMASK;

				long scoreFromDiag=packed[config.MODE_MS][row][col-1]&config.SCOREMASK;
				long scoreFromDel=packed[config.MODE_DEL][row][col-1]&config.SCOREMASK;

				long scoreMS=scoreFromDiag+config.POINTSoff_DEL;

				long scoreD=scoreFromDel+(streak==0 ? config.POINTSoff_DEL :
					streak<config.LIMIT_FOR_COST_3 ? config.POINTSoff_DEL2 :
						streak<config.LIMIT_FOR_COST_4 ?config. POINTSoff_DEL3 :
							streak<config.LIMIT_FOR_COST_5 ? config.POINTSoff_DEL4 :
								((streak&config.MASK5)==0 ? config.POINTSoff_DEL5 : 0));

				if(ref1=='N'){
					scoreMS+=config.POINTSoff_DEL_REF_N;
					scoreD+=config.POINTSoff_DEL_REF_N;
				}else if(gap){
					scoreMS+=config.POINTSoff_GAP;
					scoreD+=config.POINTSoff_GAP;
				}

				long score;
				long time;
				if(scoreMS>=scoreD){
					score=scoreMS;
					time=1;
				}else{
					score=scoreD;
					time=streak+1;
				}

				if(time>config.MAX_TIME){time=config.MAX_TIME-config.MASK5;}
				packed[config.MODE_DEL][row][col]=(score|time);
			}

			//Calculate INS score
			if(gap || (row<config.BARRIER_I1 && col>1) || (row>BARRIER_I2 && col<BARRIER_I2b)){
				packed[config.MODE_INS][row][col]=subfloor;
			}else{//Calculate INS score

				long streak=packed[config.MODE_INS][row-1][col]&config.TIMEMASK;

				long scoreFromDiag=packed[config.MODE_MS][row-1][col]&config.SCOREMASK;
				long scoreFromIns=packed[config.MODE_INS][row-1][col]&config.SCOREMASK;

				long scoreMS=scoreFromDiag+config.POINTSoff_INS;
				long scoreI=scoreFromIns+POINTSoff_INS_ARRAY[streak+1];

				long score;
				long time;
				if(scoreMS>=scoreI){
					score=scoreMS;
					time=1;
				}else{
					score=scoreI;
					time=streak+1;
				}

				if(time>config.MAX_TIME){time=config.MAX_TIME-config.MASK5;}
				packed[config.MODE_INS][row][col]=(score|time);
			}
		}
	}

	int maxCol=-1;
	int maxState=-1;
	int maxScore = std::numeric_limits<int>::min();

	for(int state=0; state<3; state++){
		for(int col=1; col<=columns; col++){
			long x=packed[state][rows][col]&config.SCOREMASK;
			if(x>maxScore){
				maxScore=x;
				maxCol=col;
				maxState=state;
			}
		}
	}
	maxScore>>=config.SCOREOFFSET;

	max[0] = rows;
	max[1] = maxCol;
	max[2] = maxState;
	max[3] = maxScore;
}

void traceback(Config &config, string &read, string &ref, long refStartLoc, long refEndLoc, int row, int col, int state, Result &r, vector<vector<vector<long> > > &packed) {

	long outPos=0;
	long gaps=0;

	long addon = refEndLoc - col;

	string out;

	while(row>0 && col>0){

		long time=packed[state][row][col]&config.TIMEMASK;
		long prev;

		if(state==config.MODE_MS){
			if(time>1){prev=(long)state;}
			else{
				long scoreFromDiag=packed[config.MODE_MS][row-1][col-1]&config.SCOREMASK;
				long scoreFromDel=packed[config.MODE_DEL][row-1][col-1]&config.SCOREMASK;
				long scoreFromIns=packed[config.MODE_INS][row-1][col-1]&config.SCOREMASK;
				if(scoreFromDiag>=scoreFromDel && scoreFromDiag>=scoreFromIns){prev=config.MODE_MS;}
				else if(scoreFromDel>=scoreFromIns){prev=config.MODE_DEL;}
				else{prev=config.MODE_INS;}
			}

			char c=toupper(read[row-1]);
			char r=toupper(ref[refStartLoc+col-1]);
			if(c==r){
				out.push_back('m');
			}else{
				out.push_back('s');
			}
			row--;
			col--;
		}else if(state==config.MODE_DEL){
			if(time>1){prev=(long)state;}
			else{
				long scoreFromDiag=packed[config.MODE_MS][row][col-1]&config.SCOREMASK;
				long scoreFromDel=packed[config.MODE_DEL][row][col-1]&config.SCOREMASK;
				if(scoreFromDiag>=scoreFromDel){prev=config.MODE_MS;}
				else{prev=config.MODE_DEL;}
			}

			char r=toupper(ref[refStartLoc+col-1]);
			if(r==config.GAPC){
				out.push_back('-');
				gaps++;
			}else{
				out.push_back('d');
			}
			col--;
		}else{
			if(time>1){prev=(long)state;}
			else{
				long scoreFromDiag=packed[config.MODE_MS][row-1][col]&config.SCOREMASK;
				long scoreFromIns=packed[config.MODE_INS][row-1][col]&config.SCOREMASK;
				if(scoreFromDiag>=scoreFromIns){prev=config.MODE_MS;}
				else{prev=config.MODE_INS;}
			}

			if(col==0){
				out.push_back('x');
			}else{
				out.push_back('i');
			}
			row--;
		}

		state=prev;
		outPos++;
	}

	if(col!=row){
		while(row>0){
			out.push_back('x');
			outPos++;
			row--;
			col--;
		}
		if(col>0){
			//do nothing
		}
	}

	reverse(out.begin(),out.end());

//	r.precise_start = r.start-SLOW_ALIGN_PADDING + col;
//	r.precise_stop = r.stop+SLOW_ALIGN_PADDING -addon;

	if(gaps==0){
		r.matchString = out;
		return;
	}

	string out3;
	int j = 0;
	for(unsigned int i=0; i<out.size(); i++){
		char c=out[i];
		if(c!=config.GAPC){
			out3.push_back(c);
			j++;
		}else{
			int lim=j+config.GAPLEN;
			for(; j<lim; j++){
				out3.push_back('d');
			}
		}
	}
	r.matchString = out3;
	return;
}


void fillLimited(Config &config, string &read, string &ref, long refStartLoc, long refEndLoc, int score, vector<long> &gaps, long max[], Result &r) {
	string gref;
	long grefLimit = (refEndLoc - refStartLoc) + 2 * config.SLOW_ALIGN_PADDING;
	//bool gapped = true;
	if((int)gaps.size() == 0) {
		for(long i = refStartLoc-config.SLOW_ALIGN_PADDING; i < refEndLoc + config.SLOW_ALIGN_PADDING; i++) {
			gref.push_back(ref[i]);
		}
	}
	else {
		grefLimit = makeGref(config, ref, gaps, (refStartLoc - config.SLOW_ALIGN_PADDING), (refEndLoc + config.SLOW_ALIGN_PADDING), gref);
	}

	int rows = read.size();
	long columns = grefLimit+1;
	vector<vector<vector<long> > > packed;

	for(int matrix=0; matrix<3; matrix++){
		vector<long> row;
		for(int i = 0; i <= columns; i++) {
			row.push_back(0);
		}
		vector<vector<long> > mat;
		mat.push_back(row);
		for(int i=1; i<=rows; i++){
			vector<long> row;
			for(int j=0; j<columns+1; j++){
				row.push_back(config.BADoff);
			}
			mat.push_back(row);
		}
		for(int i=0; i<=rows; i++){
			int prevScore=(i<2 ? 0 : mat[i-1][0]);
			int score=prevScore+POINTSoff_INS_ARRAY[i];
			mat[i][0]=score;
		}
		packed.push_back(mat);
	}

	//cout << read << endl;
	//cout << gref << endl;

	fillUnlimited(config, read, gref, 0, grefLimit, max, packed);
	traceback(config, read, gref, 0, grefLimit, max[0], max[1], max[2], r, packed);
	//if(r.gapArray.size() == 0) r.precise_stop--;


}

void makeMatchStringForSite(Config &config, SiteScore ss, string &read, long *sizes, long *sites, Result &r, string &whole_genome, int threadId) {

	if(ss.perfect) {
		r.start = ss.start;
		r.stop = ss.stop;
	//	r.precise_start = r.start;
	//	r.precise_stop = r.stop;
		for(unsigned int i = 0; i < read.size(); i++) {
			r.matchString.push_back('m');
		}
		return;
	}

	string match;
	int score = scoreNoIndelsAndMakeMatchString(config, read, whole_genome, ss.start, match);
	int imperfectScore = maxImperfectScore(config, read);
	if(score > imperfectScore) {
		//cout << "perfect" << endl;
		//cout << match << endl;
		r.start = ss.start;
		r.stop = ss.stop;
	//	r.precise_start = r.start;
	//	r.precise_stop = r.stop;
		r.matchString = match;
		return;
	}

	long max[4];
	fillLimited(config, read, whole_genome, ss.start, ss.stop, score, ss.gapArray, max, r);
}

void makeMatchString(vector<SiteScore> &results, string &read, string &read_reverse, long *sizes, long *sites, vector<Result> &resultsFinal, string &whole_genome, int threadId, int br, int maxScore) {
	int max_score = -9999;
	SiteScore ss;

	for(unsigned int i = 0; i < results.size(); i++) {
		SiteScore s = results[i];
		if(s.score > max_score) {
			max_score = s.score;
			SiteScore ress = results[i];
			ss = results[i];
		}
	}
	Result r = Result(br, ss.start, ss.stop, ss.score, maxScore);
	r.gapArray = ss.gapArray;

/*
	if(ss.strand == 0) {
		makeMatchStringForSite(ss, read, sizes, sites, r, whole_genome, threadId);
	}
	else {
		makeMatchStringForSite(ss, read_reverse, sizes, sites, r, whole_genome, threadId);
	}
*/
	resultsFinal.push_back(r);
}

void filterBadReads(Config &config, int score, vector<SiteScore> &results, string &read, string &read_reverse, long *sizes, long *sites, vector<vector<Result> > &databaseResults, vector<Result> &resultsFinal, string &whole_genome, int threadId, int br, int maxScore) {
	int max_score = -9999;
	SiteScore ss;

	vector<Result> filtered;
	for(int i = 0; i < results.size(); i++) {
		SiteScore s = results[i];
		Result r = Result(br, s.start, s.stop, s.score, s.score);
		if(s.strand == 0) {
			makeMatchStringForSite(config, s, read, sizes, sites, r, whole_genome, threadId);
		}
		else {
			makeMatchStringForSite(config, s, read_reverse, sizes, sites, r, whole_genome, threadId);
		}
		r.gapArray = s.gapArray;
		resultsFinal.push_back(r);
	}
/*
	for(unsigned int i = 0; i < results.size(); i++) {
		SiteScore s = results[i];
		if(s.score > max_score) {
			max_score = s.score;
			SiteScore ress = results[i];
			ss = results[i];
		}
		if (s.score > score * 0.0) {
			Result r = Result(br, s.start, s.stop, s.score, s.score);
			cout << "start " << s.start << endl;
			r.gapArray = s.gapArray;
			filtered.push_back(r);
		}
	}
	Result r = Result(br, ss.start, ss.stop, ss.score, maxScore);
	r.gapArray = ss.gapArray;
	resultsFinal.push_back(r);
	*/
	databaseResults.push_back(filtered);
}
