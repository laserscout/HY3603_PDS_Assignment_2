#ifndef _FIND_MEDIAN_H_
#define _FIND_MEDIAN_H_

float selection(float *array,double **data,int number);
float masterPart(int noProcesses,int processId,int size,int partLength,float *numberPart,double **data,MPI_Comm comm);
void slavePart(int processId,int partLength, float *numberPart,int size, double **data,MPI_Comm comm);
void validation(float median,int partLength,int size,float *numberPart,int processId,MPI_Comm comm);
void validationST(float median, int size, float *numberPart);

#endif
