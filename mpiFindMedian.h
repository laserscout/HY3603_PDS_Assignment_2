#ifndef _MPI_FIND_MEDIAN_H_
#define _MPI_FIND_MEDIAN_H_

float selection(float *array,int number);
void validationST(float median, int size, float *numberPart);
float masterPart(int noProcesses,int processId,int size,int partLength,float *numberPart,MPI_Comm comm);
void slavePart(int processId,int partLength, float *numberPart,int size, MPI_Comm comm);

#endif
