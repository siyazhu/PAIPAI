#include <string>
#include "element.h"
using namespace std;
int elementnum(string name) {
	int i;
	string nametemp;
	for (i = 0; i <= 85; i++) {
		nametemp = string(PERIODICTABLE[i]);
		if (name == nametemp)
			return (i+1);
	}
	return 0;
}


