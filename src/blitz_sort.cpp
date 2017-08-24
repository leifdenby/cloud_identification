#include "blitz/array.h"


// partitioning for the quicksort algorithm
int partition(blitz::Array<int, 2> &a, int p, int r) {
  int j, temp;
  
  int x = a(r,0);
  int i = p -1;
  
  for(j = p;j<r;j++){
    if(a(j,0)<= x)
    {
      i++;
      temp = a(i,0);
      a(i,0) = a(j,0);
      a(j,0) = temp;
      temp = a(i,1);
      a(i,1) = a(j,1);
      a(j,1) = temp;
      temp = a(i,2);
      a(i,2) = a(j,2);
      a(j,2) = temp;
      temp = a(i,3);
      a(i,3) = a(j,3);
      a(j,3) = temp;
    }
  }
  temp = a(r,0);
  a(r,0) = a(i +1,0);
  a(i+1,0) = temp;
  temp = a(r,1);
  a(r,1) = a(i +1,1);
  a(i+1,1) = temp;
  temp = a(r,2);
  a(r,2) = a(i +1,2);
  a(i+1,2) = temp;
  temp = a(r,3);
  a(r,3) = a(i +1,3);
  a(i+1,3) = temp;
    
  return i+1;
}

// actual quicksort algorithm
void quick_sort(blitz::Array<int, 2> &a, int l, int r)
{
  int j;

  if( l < r ) 
  {
    // divide and conquer
    j = partition( a, l, r);
    quick_sort( a, l, j-1);
    quick_sort( a, j+1, r);
  }
}
