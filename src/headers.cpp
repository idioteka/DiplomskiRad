
#include <ctype.h>

int getCodeFromBase(char base) {
	char lower_base = tolower(base);
	if(lower_base == 'a') {
		return 0;
	} else if(lower_base == 'c') {
		return 1;
	}
	else if(lower_base == 'g') {
		return 2;
	}
	else if(lower_base == 't') {
		return 3;
	}
	else return -1;
}
