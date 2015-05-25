
/*
 *
 * PARAMS FUNCTIONS
 *
 */

void checkParameter(Config &config, string command, string value);
int readParams(int argc, char *argv[], Config &config, Info &info);
bool createDirectories(Config &config);

void updatePrecision2(Config &config, int prec);
void updatePrecision(Config &config, bool prec);
