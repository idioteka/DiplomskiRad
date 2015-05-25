#include <assert.h>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <algorithm>
#include <vector>
#include <fstream>
#include <sstream>
#include <string>
#include <map>
#include <iostream>
#include <string>
#include <ctime>
#include <cmath>
#include <limits>
#include <cstdlib>
#include <pthread.h>
#include <sys/time.h>
#include <sys/stat.h>

using namespace std;

#define SSTR( x ) dynamic_cast< std::ostringstream & >( \
        ( std::ostringstream() << std::dec << x ) ).str()

int getCodeFromBase(char base);
