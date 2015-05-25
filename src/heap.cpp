
#include "headers.h"
#include "heap.h"

Heap::Heap() {
	size_ = 0;
}

int Heap::compare(Triplet &t1, Triplet &t2) {
	int x = t1.site - t2.site;
	if(x != 0) return x;
	return t1.column - t2.column;
}

void Heap::percDown(int loc){

	if(loc == 1) {
		return;
	}

	int next = loc/2;
	Triplet a = array[loc];
	Triplet b = array[next];

	while(loc > 1 && compare(a, b) < 0){
		array[loc] = b;
		loc = next;
		next = next/2;
		b = array[next];
	}

	array[loc]=a;
}

void Heap::percUp(int loc){
	int next1 = loc * 2;
	int next2 = next1 + 1;
	if(next1 > size_) {
		return;
	}
	Triplet a = array[loc];
	Triplet b = array[next1];
	if(next2 > size_) {
		if(compare(a, b) > 0){
			array[next1] = a;
			array[loc] = b;
			percUp(next1);
		}
	}
	else {
		Triplet c = array[next2];
		if(compare(b, c) < 1) {
			if(compare(a, b) > 0){
				array[next1] = a;
				array[loc] = b;
				percUp(next1);
			}
		}
		else{
			if(compare(a, c)>0){
				array[next2] = a;
				array[loc] = c;
				percUp(next2);
			}
		}
	}
}

void Heap::initHeap(int max_size) {
	CAPACITY = max_size;
	Triplet t;
	array.clear();
	array.push_back(t);

}

void Heap::add(Triplet &t) {
	size_++;
	array.push_back(t);
	percDown(size_);
}

Triplet Heap::peek() {
	if(size_ == 0) {
		return array[0];
	}
	return array[1];
}
Triplet Heap::poll(){
	if(size_ == 0) {
		return array[0];
	}
	Triplet t = array[1];
	array[1] = array[size_];
	array.erase(array.begin()+size_);
	size_--;
	if(size_ > 0) {
		percUp(1);
	}
	return t;
}

bool Heap::isEmpty() {
	return size_ == 0;
}

void Heap::clear() {
	size_ = 0;
}

int Heap::size() {
	return size_;
}

