/*
  The MIT License (MIT)

  Copyright (c) 2014

  Athanassios Kintsakis
  Contact
  athanassios.kintsakis@gmail.com
  akintsakis@issel.ee.auth.gr


  Permission is hereby granted, free of charge, to any person obtaining a copy
  of this software and associated documentation files (the "Software"), to deal
  in the Software without restriction, including without limitation the rights
  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
  copies of the Software, and to permit persons to whom the Software is
  furnished to do so, subject to the following conditions:

  The above copyright notice and this permission notice shall be included in all
  copies or substantial portions of the Software.

  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
  SOFTWARE.
*/

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <assert.h>
#include "findMedian.h"

MPI_Status Stat;
void partition(float *array, int elements, float pivot, float **arraysmall, float **arraybig, int *endsmall, int *endbig, double **data);
float selection(float *array, double **data, int number);

void removeElement(int *array, int *size, int element)
{
  int i;
  int flag=0;
  for(i=0;i<*size;i++)
    {
      if(flag==1)
	array[i]=array[i+1];
      if(array[i]==element&& flag==0)
        {
	  array[i]=array[i+1];
	  flag=1;
        }
    }
  *size=*size-1;
}

void swap_values(float *array,double **data_array,int x,int y)
{
  float temp;
  temp=array[x];
  array[x]=array[y];
  array[y]=temp;

  double *data_ptr_temp;
  data_ptr_temp=data_array[x];
  data_array[x]=data_array[y];
  data_array[y]=data_ptr_temp;

}

float masterPart(int noProcesses,int processId,int size,int partLength,float *numberPart,double **data,MPI_Comm comm) 
{
  int elements,i,keepBigSet,sumSets,finalize,randomNode,k;
  float pivot, tempPivot, median;
  int endSmall=0;
  int dropoutFlag=0;
  int endBig=0;
  float *arraySmall,*arrayBig,*arrayToUse;
  int *activeNodes;
  int activeSize=noProcesses;
  int stillActive=1;
  int oldSumSets=-1;
  int checkIdentical=0;
  int useNewPivot=0;
  int *pivotArray;
  k=(int)size/2+1; 
  elements=partLength;
  activeNodes=(int *)malloc(noProcesses*sizeof(int));  
							 
  arrayToUse=numberPart;
  pivotArray=(int*)malloc(noProcesses*sizeof(int));  
						       
  for(i=0;i<activeSize;i++)
    {
      activeNodes[i]=i;
    }
  int randomCounter=0;
  int randomCounter2=0;

  for(;;)   
    {
      int counter=0;
      useNewPivot=0;
      if(stillActive==1&&checkIdentical!=0)  
        {
	  for(i=0;i<elements;i++)
            {
	      if(pivot==arrayToUse[i])
		counter++;
	      else
                {
		  useNewPivot=1;
		  tempPivot=arrayToUse[i];
		  break;
                }
            }
        }
      if(checkIdentical!=0)
        {
	  int useNewPivotMax=0;
	  MPI_Reduce(&useNewPivot,&useNewPivotMax,1,MPI_INT,MPI_MAX,0,comm); 
	  if(useNewPivotMax!=1)    
            {
	      median=pivot;
	      finalize=1;
	      MPI_Bcast(&finalize,1,MPI_INT,0,comm); 
	      free(pivotArray);
	      return median;
            }
	  else
            {
	      finalize=0;
	      int useit=0;
	      randomCounter2++;
	      MPI_Bcast(&finalize,1,MPI_INT,0,comm);
	      MPI_Gather(&useNewPivot, 1, MPI_INT, pivotArray, 1, MPI_INT, 0, comm);
	      for(i=0;i<activeSize;i++)
                {
		  if(pivotArray[i]==1)
                    {
		      if((randomCounter2>1)&&(randomNode!=activeNodes[i]))  
			{
			  useit=1;
			  randomNode=activeNodes[i];
			  randomCounter2=0;
			  break;
                        }
		      else if(randomCounter2<2)
                        {
			  useit=1;
			  randomNode=activeNodes[i];
			  break;
                        }
                    }
                }
	      if(useit!=0)
		useNewPivot=1;
	      else
		useNewPivot=0;
            }
        }
      if(useNewPivot!=0)
	MPI_Bcast(&randomNode,1,MPI_INT,0,comm);  
      if(useNewPivot==0)  
        {
	  if(randomCounter>=activeSize)
	    randomCounter=0; 
	  randomNode=activeNodes[randomCounter];
	  randomCounter++;			
	  MPI_Bcast(&randomNode,1,MPI_INT,0,comm);   
        }
      if(randomNode==processId)  
	{
	  if(useNewPivot==0)
            {
	      srand(time(NULL));
	      pivot=arrayToUse[rand() % elements];
	      MPI_Bcast(&pivot,1,MPI_FLOAT,0,comm); 
	    }
	  else
            {
	      MPI_Bcast(&tempPivot,1,MPI_FLOAT,0,comm); 
	      pivot=tempPivot;
            }
        }
      else 
	MPI_Bcast(&pivot,1,MPI_FLOAT,randomNode,comm);  
      if(stillActive==1)  
        {
	  partition(arrayToUse,elements,pivot,&arraySmall,&arrayBig,&endSmall,&endBig,data); 
        }
      else  
        {
	  endBig=0;
	  endSmall=0;
        }
      sumSets=0;
      MPI_Reduce(&endBig,&sumSets,1,MPI_INT,MPI_SUM,0,comm);  
      MPI_Bcast(&sumSets,1,MPI_INT,0,comm);
      if(oldSumSets==sumSets)
	checkIdentical=1;
      else
        {
	  oldSumSets=sumSets;
	  checkIdentical=0;
        }
      if(sumSets>k)   
	{
	  keepBigSet=1; 
	  if(endBig==0)
	    dropoutFlag=1; 
	  else
            {
	      arrayToUse=arrayBig; 
	      elements=endBig; 
            }
	}
      else if(sumSets<k) 
	{
	  keepBigSet=0;
	  k=k-sumSets;
	  if(endSmall==0)
	    dropoutFlag=1; 
	  else
	    {
	      arrayToUse=arraySmall; 
	      elements=endSmall;
	    }
	}
      else  
	{
	  median=pivot;
	  finalize=1; 
	  MPI_Bcast(&finalize,1,MPI_INT,0,comm);
	  free(pivotArray);
	  return median;
	}
      finalize=0; 
      MPI_Bcast(&finalize,1,MPI_INT,0,comm);	
      MPI_Bcast(&keepBigSet,1,MPI_INT,0,comm);    
      if(dropoutFlag==1 && stillActive==1) 
        {
	  stillActive=0;
	  removeElement(activeNodes, &activeSize, 0);
        }
      int flag;
      for(i=0;i<activeSize;i++)
        {
	  if(activeNodes[i]!=0)
            {
	      MPI_Recv(&flag,1,MPI_INT,activeNodes[i],1,comm,&Stat);  
	      if(flag==1)
		removeElement(activeNodes, &activeSize, activeNodes[i]);
            }
        }
    }
}

void slavePart(int processId,int partLength, float *numberPart,int size, double **data, MPI_Comm comm)  
{
  int dropoutflag,elements,i,sumSets,finalize,keepBigSet,randomNode;
  float pivot, tempPivot;
  int endSmall=0;
  int endBig=0;
  float *arraySmall,*arrayBig,*arrayToUse;
  arrayToUse=numberPart;
  elements=partLength;
  int stillActive=1;
  int *pivotArray;
  int oldSumSets=-1;
  int checkIdentical=0;
  int useNewPivot;
  for(;;)
    {
      finalize=0;
      int counter=0;
      useNewPivot=0;
      if(stillActive==1&&checkIdentical!=0)  
        {
	  for(i=0;i<elements;i++)
            {
	      if(pivot==arrayToUse[i])
		counter++;
	      else
                {
		  useNewPivot=1;
		  tempPivot=arrayToUse[i];
		  break;
                }
            }
        }
      if(checkIdentical!=0)
        {
	  int useNewPivotMax=0;
	  MPI_Reduce(&useNewPivot,&useNewPivotMax,1,MPI_INT,MPI_MAX,0,comm);
	  MPI_Bcast(&finalize,1,MPI_INT,0,comm);
	  if(finalize==1)
            {
	      return ;
            }
	  else
            {
	      MPI_Gather(&useNewPivot, 1, MPI_INT, pivotArray, 1, MPI_INT, 0, comm);
            }
        }
      MPI_Bcast(&randomNode,1,MPI_INT,0,comm); 
      if(randomNode!=processId) 
	MPI_Bcast(&pivot,1,MPI_FLOAT,randomNode,comm);       
      else if(randomNode==processId) 
        {
	  if(useNewPivot==0)
            {
	      srand(time(NULL));
	      pivot=arrayToUse[rand() % elements];
	      MPI_Bcast(&pivot,1,MPI_FLOAT,processId,comm); 
            }
	  else
            {
	      MPI_Bcast(&tempPivot,1,MPI_FLOAT,processId,comm); 
	      pivot=tempPivot;
            }
        }
      if(stillActive==1)   
        {
	  partition(arrayToUse,elements,pivot,&arraySmall,&arrayBig,&endSmall,&endBig,data);
        }
      else
        {
	  endBig=0;
	  endSmall=0;
        }
      sumSets=0;
      MPI_Reduce(&endBig,&sumSets,1,MPI_INT,MPI_SUM,0,comm); 
      MPI_Bcast(&sumSets,1,MPI_INT,0,comm);
      if(oldSumSets==sumSets)
	checkIdentical=1;
      else
        {
	  oldSumSets=sumSets;
	  checkIdentical=0;
        }
      MPI_Bcast(&finalize,1,MPI_INT,0,comm);
      if(finalize==1)
        {
	  return ;
        }
      MPI_Bcast(&keepBigSet,1,MPI_INT,0,comm);
      if(stillActive==1)
        {
	  if(keepBigSet==1)
            {
	      if(endBig==0)
		dropoutflag=1;
	      else
                {
		  arrayToUse=arrayBig;
		  elements=endBig;
                }
            }
	  else if(keepBigSet==0)
            {
	      if(endSmall==0)
		dropoutflag=1;
	      else
                {
		  arrayToUse=arraySmall;
		  elements=endSmall;
                }
            }
        }
      if(dropoutflag==1 && stillActive==1)
        {
	  MPI_Send(&dropoutflag,1,MPI_INT,0,1,comm); 
	  stillActive=0;
        }
      else if(stillActive==0)
        {
	  dropoutflag=-1;
	  MPI_Send(&dropoutflag,1,MPI_INT,0,1,comm); 
        }
      else
        {
	  dropoutflag=0;
	  MPI_Send(&dropoutflag,1,MPI_INT,0,1,comm); 
        }
    }
}

void partition (float *array, int elements, float pivot, float **arraysmall, float **arraybig, int *endsmall, int *endbig, double **data)
{
  int right=elements-1;
  int left=0;
  int pos;
  if(elements==1)
    {
      if(pivot>array[0])
        {
	  *endsmall=1;  
	  *endbig=0;   
	  *arraysmall=array;   
	  *arraybig=NULL;
        }
      else if(pivot<=array[0])
        {
	  *endsmall=0;    
	  *endbig=1;
	  *arraysmall=NULL;
	  *arraybig=array;
        }
    }
  else if(elements>1)
    {
      while(left<right)
        {
	  while(array[left]<pivot)
            {
	      left++;
	      if(left>=elements)
                {
		  break;
                }
            }
	  while(array[right]>=pivot)
            {
	      right--;
	      if(right<0)
                {
		  break;
                }
            }
	  if(left<right)
            {
	      swap_values(array,data,left,right);
            }
        }
      pos=right;
      if(pos<0)                   
        {                               
	  *arraysmall=NULL;           
        }                               
      else
        {
	  *arraysmall=array;
        }
      *endsmall=pos+1;
      *arraybig=&array[pos+1];
      *endbig=elements-pos-1;
    }
}

float selection(float *array,double **data,int number)
{
  float *arraybig;
  float *arraysmall;
  int endsmall=0;
  int endbig=0;
  float *arraytobeused;
  int i;
  int counter=0;
  int k;
  float pivot;
  float median;
  k=(int)number/2+1;
  arraytobeused=array;
  for(;;)
    {
      pivot=arraytobeused[rand() % number];
      partition(arraytobeused,number,pivot,&arraysmall,&arraybig,&endsmall,&endbig,data);
      if(endbig>k)
        {
	  number=endbig;
	  arraytobeused=arraybig;
	  for(i=0;i<endbig;i++)
            {
	      if(pivot==arraybig[i])
		counter++;
	      else
		break;
            }
	  if(counter==endbig)
            {
	      median=arraybig[0];
	      break;
            }
	  else
	    counter=0;
            
        }
      else if(endbig<k)
        {
	  number=endsmall;
	  arraytobeused=arraysmall;
	  k=k-endbig;
        }
      else
        {
	  median=pivot;
	  break;
        }
    }
  return median;
}
