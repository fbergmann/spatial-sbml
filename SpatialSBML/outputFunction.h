#ifndef OUTPUTFUNCTION_H_
#define OUTPUTFUNCTION_H_

void outputTimeCource(FILE *gp, Model *model, vector<variableInfo*> &varInfoList, vector<const char*> memList, variableInfo *xInfo, variableInfo *yInfo, variableInfo *zInfo, double *sim_time, double end_time, double dt, double range_max, int dimension, int Xindex, int Yindex, int Zindex, double Xsize, double Ysize, double Zsize, int file_num, string fname);

#endif
