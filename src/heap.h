

struct Triplet {
	int column;
	int row;
	int site;
	Triplet() {
		column = -1;
		row = -1;
		site = -1;
	}
	Triplet(int column_, int row_, int site_) {
		column = column_;
		row = row_;
		site = site_;
	}
};

bool isNull(Triplet &t);


void initHeap(int max_size);
void add(Triplet &t);
Triplet peek();
Triplet poll();
bool isEmpty();
void clear();
int size();
