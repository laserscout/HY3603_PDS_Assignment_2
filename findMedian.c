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
void partition2 (float *array, int elements, float pivot, float **arraysmall, float **arraybig, int *endsmall, int *endbig, double **data);
float selection2(float *array, double **data, int number);

/****Swaps two values in an array****/
void swap_values2(float *array,double **data_array,int x,int y)
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

/****Partitions the Array into larger and smaller than the pivot values****/
void partition2 (float *array, int elements, float pivot, float **arraysmall, float **arraybig, int *endsmall, int *endbig, double **data)
{
  int right=elements-1;
  int left=0;
  int pos;
    if(elements==1)
    {
        if(pivot>array[0])
        {
            *endsmall=1;  //One value in the small part
            *endbig=0;   //Zero on the big one
            *arraysmall=array;   //There is no big array therefore NULL value
            *arraybig=NULL;
        }
        else if(pivot<=array[0])
        {
            *endsmall=0;    //The exact opposite of the above actions.
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
	      swap_values2(array,data,left,right);
            }
        }
        pos=right;
        if(pos<0)                   //Arrange the arrays so that they are split into two smaller ones.
        {                               //One containing the small ones. And one the big ones.
            *arraysmall=NULL;           //However these arrays are virtual meaning that we only save the pointer values of the beging and end
        }                               //of the "real" one.
        else
        {
            *arraysmall=array;
        }
        *endsmall=pos+1;
        *arraybig=&array[pos+1];
        *endbig=elements-pos-1;
    }
}


/***==============================================***/
/***==============================================***/
/***=============SERIAL SELECTION==============***/
/***==============================================***/
/***==============================================***/

float selection2(float *array,double **data,int number)
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
        partition2(arraytobeused,number,pivot,&arraysmall,&arraybig,&endsmall,&endbig,data);
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
            //end of count equals
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
