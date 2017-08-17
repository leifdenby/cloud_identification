// The code to do cloud numbering
// This code is rather procedural (=faster?)
// Works on both MONC and UM data given the right parameters

// Compile on MONSOON: g++ -pthread -fno-strict-aliasing -DNDEBUG -g -fwrapv -fPIC -I/opt/python/gnu/2.7.9/lib/python2.7/site-packages/scipy/weave/blitz -I/opt/cray/netcdf-hdf5parallel/4.3.2/CRAY/83/include   -I/opt/cray/hdf5/1.8.13/CRAY/83/include -c cpp_identify_seungbu.cpp -o cpp_identify_seungbu.o -pthread -O6 -march=native -mtune=native -funroll-all-loops -fomit-frame-pointer -march=native -mtune=native -msse4 -ftree-vectorize -ftree-vectorizer-verbose=5 -ffast-math -funroll-loops -ftracer
// Link:  g++ -o cpp_identify_seungbu.exe cpp_identify_seungbu.o -L/opt/cray/netcdf-hdf5parallel/4.3.2/CRAY/83/lib -lnetcdf

// Compile on my local system: g++ -pthread -fno-strict-aliasing -DNDEBUG -g -fwrapv -fPIC -I/usr/lib/python2.7/dist-packages/scipy/weave/blitz -I/usr/include/ -c cpp_identify_seungbu.cpp -o cpp_identify_seungbu.o -pthread -O6 -march=native -mtune=native -funroll-all-loops -fomit-frame-pointer -march=native -mtune=native -msse4 -ftree-vectorize -ftree-vectorizer-verbose=5 -ffast-math -funroll-loops -ftracer
// Link: g++ -o cpp_identify_seungbu.exe cpp_identify_seungbu.o -L/usr/lib/ -lnetcdf

// This is very recent code, and there are no publications with it yet. 
// It basically works by first assigning each point that fulfils a masking criterion (e.g. each cloudy point) to a local maximum. 
// The assignment to maxima is done through a steepest gradient approach. 
// Subsequently, a merging algorithm is applied to get rid of the smaller local peaks.

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <netcdf.h>
#include <stdint.h>
#include "blitz/array.h"

#define ERRCODE 2
#define ERR(e) {printf("Error: %s\n", nc_strerror(e)); exit(ERRCODE);}

// PROPERTIES OF INPUT
static const int imax = 514;
static const int jmax = 514;
static const int kmax = 2;
static const float inv_scaling_parameter=1.0/3200.0; // data conversion parameter
                                                       // corresponding to scaling in pre-processing script
static const bool lperiodic = true; // is the domain periodic?

// PARAMETERS FOR CLUSTERING ALGORITHM
static const float mincolratio = 0.7; // fractional height used for deciding whether to merge cols
                                      // 0 corresponds to merging everything
                                      // 1 corresponds to not merging (steepest gradient association
                                      // with local maximum

// PARAMETERS WHICH DETERMINE MEMORY FOOTPRINT
// By what factor far can we make arrays smaller? (reduce memory footprint)
static const int blobfac = 32; // number of blobs wrt number of points
static const int borderfac = 8; // number of border points wrt number of points
static const int colfac = 32; // number of cols wrt number of points

typedef uint32_t indexint; //use unsigned array indices

// fill a blitz array with "short integer" netcdf data
int blitzncshort(int ncId, const char *nomVariable,
blitz::Array<short,3> & tab)
{
  int status, varId, nbDims;
  int *dimIds;
  size_t start[] = {0,0,0};
  size_t count[3];

  status = nc_inq_varid(ncId, nomVariable, &varId);
  if (status != NC_NOERR)
    std::cout << nomVariable << " is absent" << std::endl ;
  else
    {
      status = nc_inq_varndims(ncId, varId, &nbDims);

      if(nbDims != 3)
    {
      fprintf(stderr,"Error in data, the variable %s should have 3 dimensions\n", nomVariable);
      exit(-1);
    }

      dimIds = new int[nbDims*sizeof(*dimIds)];

      status = nc_inq_vardimid(ncId, varId, dimIds);

      delete [] dimIds;

      tab.resize(imax, jmax, kmax);

      count[0] = imax ; count[1] = jmax; count[2] = kmax;
      status = nc_get_vara_short(ncId, varId, start, count, tab.data());
    }
  return 0;
}

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

// MAIN PROGRAM    
int main()
  {
    int retval; // error terurn value
    int ncid,ncid2,ncid3,ncid4; //netcdf ids (decided not to reuse these)
    int varid,varid2,varid3,varid4; //variable ids (decided not to reuse these)
    int err; 
    
    int itarget,jtarget,ktarget,ltarget; //calculating back target i,j,k,l that corresponds to maxima
    indexint counter; // counts cells with certain properties (used for sanity checks)
    char tempdirection; // direction of steepest ascent temp-vars are used in inner loops
    indexint tempdata; // used in inner loop
    int tempmax,tempfield; // used in inner loop
    bool moreswaps; // used to check if all identities have been assigned, or further checking necessary

    // array that holds "cloud" mask as short (temporarily)
    blitz::Array<short,3> maskshort(imax,jmax,kmax);

    /* Open the file. NC_NOWRITE tells netCDF we want read-only access
     * to the file.*/
    if ((retval = nc_open("mask_field.nc", NC_NOWRITE, &ncid)))
      ERR(retval);

    err=blitzncshort(ncid,"maskext",maskshort);

    /* Close the file, freeing all resources. */
    if ((retval = nc_close(ncid)))
      ERR(retval);
    
    // array that holds the mask as boolean (including halo cells)
    blitz::Array<bool,3> maskext(imax,jmax,kmax);
    maskext=false;

    // transfer mask elements to boolean and count them
    counter=1;
    for (int i=1; i<imax-1; ++i) {
    for (int j=1; j<jmax-1; ++j) {
    for (int k=0; k<kmax; ++k) {
      if(maskshort(i,j,k)>0) {
        maskext(i,j,k)=true;
        counter=counter+1;
      }
    }
    }
    }
    printf("number of cloudy cells +1 %d \n",counter); 

    maskshort.free();

    // array that holds field values
    blitz::Array<short,3> fieldext(imax,jmax,kmax);
    fieldext=0;

    /* Open the file. NC_NOWRITE tells netCDF we want read-only access
     * to the file.*/
    if ((retval = nc_open("mask_field.nc", NC_NOWRITE, &ncid2)))
       ERR(retval);

    err=blitzncshort(ncid2,"fieldext",fieldext);

    /* Close the file, freeing all resources. */
    if ((retval = nc_close(ncid2)))
       ERR(retval);

    // array that hold actual numbers
    blitz::Array<indexint,3> dataext(imax,jmax,kmax);
    dataext=0;
    
    // array that holds direction of steepest ascent
    blitz::Array<char,3> direction(imax,jmax,kmax);
    direction=0;
   
    // first numbering
    for (int i=1; i<imax-1; ++i) {
    for (int j=1; j<jmax-1; ++j) {
    for (int k=0; k<kmax; ++k) {
      if(maskext(i,j,k)==true) {
        dataext(i,j,k)=i*jmax*kmax+j*kmax+k;
      }
    }
    }
    }
    
    // halo swap to get "inner domain values" into halos
    if(lperiodic==true) {
      for (int i=0; i<imax; ++i) {
      for (int k=0; k<kmax; ++k) {
    	if(maskext(i,jmax-2,k)==true) {
    	  dataext(i,0,k)=dataext(i,jmax-2,k);
    	  fieldext(i,0,k)=fieldext(i,jmax-2,k);
    	  maskext(i,0,k)=maskext(i,jmax-2,k);
    	}
    	if(maskext(i,1,k)==true) {
    	  dataext(i,jmax-1,k)=dataext(i,1,k);
    	  fieldext(i,jmax-1,k)=fieldext(i,1,k);
    	  maskext(i,jmax-1,k)=maskext(i,1,k);
    	}
      }
      }

      for (int j=0; j<jmax; ++j) {
      for (int k=0; k<kmax; ++k) {
    	if(maskext(imax-2,j,k)==true) {
    	  dataext(0,j,k)=dataext(imax-2,j,k);
    	  fieldext(0,j,k)=fieldext(imax-2,j,k);
    	  maskext(0,j,k)=maskext(imax-2,j,k);
    	}
    	if(maskext(1,j,k)==true) {
    	  dataext(imax-1,j,k)=dataext(1,j,k);
    	  fieldext(imax-1,j,k)=fieldext(1,j,k);
    	  maskext(imax-1,j,k)=maskext(1,j,k);
    	}
      }
      }
    }
   
    // the sum of directions can be used to do checks on a 
    // domain that repeats itself in the horizontal directions
    int sumdir; 
    sumdir=0;

    // detection of steepest gradient direction in inner domain
    // 0 corresponds to local maxima
    for (int i=1; i<imax-1; ++i) {
    for (int j=1; j<jmax-1; ++j) {
    for (int k=0; k<kmax; ++k) {
      if(maskext(i,j,k)==true) {
        tempdirection=0;
        tempmax=fieldext(i,j,k);
        if(maskext(i-1,j,k)==true) {
          tempfield=fieldext(i-1,j,k);   
          if(tempfield>tempmax) {
            tempdirection=1;
            tempmax=tempfield;
          }
        }
	if(maskext(i,j-1,k)==true) {
          tempfield=fieldext(i,j-1,k);   
          if(tempfield>tempmax) {
            tempdirection=2;
            tempmax=tempfield;
          }
        }
	if(k>0){
          if(maskext(i,j,k-1)==true) {
            tempfield=fieldext(i,j,k-1);   
            if(tempfield>tempmax) {
              tempdirection=3;
              tempmax=tempfield;
            }
          }
        }
        if(maskext(i+1,j,k)==true) {
          tempfield=fieldext(i+1,j,k);   
          if(tempfield>tempmax) {
            tempdirection=4;
            tempmax=tempfield;
          }
        }  
	if(maskext(i,j+1,k)==true) {
          tempfield=fieldext(i,j+1,k);   
          if(tempfield>tempmax) {
            tempdirection=5;
            tempmax=tempfield;
          }
        }
        if(k<kmax-1){   
          if(maskext(i,j,k+1)==true) {
            tempfield=fieldext(i,j,k+1);   
            if(tempfield>tempmax) {
              tempdirection=6;
              tempmax=tempfield;
            }
          }
        }                             
        direction(i,j,k)=tempdirection;
        sumdir=sumdir+tempdirection;
      }
    }
    }
    }
    printf("sum of directions %d \n",sumdir);
    fieldext.free();
    
    // assign points to local maximum   
    // order: inner domain forward sweep, halo update, inner domain reverse sweep, haloupdate

    // forward sweep
    moreswaps=true; 
    while(moreswaps==true) {
      moreswaps=false; //reset swap

      for (int i=1; i<imax-1; ++i) {
      for (int j=1; j<jmax-1; ++j) {
      for (int k=0; k<kmax; ++k) {
        if(maskext(i,j,k)==true) {
          tempdirection=direction(i,j,k);
          if(tempdirection>0 and tempdirection<4) {
            if(tempdirection==1) {
              if(dataext(i,j,k)!=dataext(i-1,j,k)) {   
                dataext(i,j,k)=dataext(i-1,j,k);
                moreswaps=true;
              }
            }
            else if(tempdirection==2) {
              if(dataext(i,j,k)!=dataext(i,j-1,k)) {   
                dataext(i,j,k)=dataext(i,j-1,k);
                moreswaps=true;
              }
            }
            else if(tempdirection==3) {
              if(dataext(i,j,k)!=dataext(i,j,k-1)) {   
                dataext(i,j,k)=dataext(i,j,k-1);
                moreswaps=true;
              }
            }        
          }        
        }
      }
      }
      }    

      // halo update
      if(lperiodic==true) {
     	for (int i=0; i<imax; ++i) {
     	for (int k=0; k<kmax; ++k) {
     	  if(maskext(i,jmax-2,k)==true) {
     	    dataext(i,0,k)=dataext(i,jmax-2,k);
     	  }
     	  if(maskext(i,1,k)==true) {
     	    dataext(i,jmax-1,k)=dataext(i,1,k);
     	  }
     	}
     	}
  
     	for (int j=0; j<jmax; ++j) {
     	for (int k=0; k<kmax; ++k) {
     	  if(maskext(imax-2,j,k)==true) {
     	    dataext(0,j,k)=dataext(imax-2,j,k);
     	  }
     	  if(maskext(1,j,k)==true) {
     	    dataext(imax-1,j,k)=dataext(1,j,k);
     	  }
     	}
     	}
      }
      
      // reverse sweep
      for (int i=imax-2; i>= 1; --i) {
      for (int j=jmax-2; j>= 1; --j) {
      for (int k=kmax-1; k>= 0; --k) {
        if(maskext(i,j,k)==true) {
          tempdirection=direction(i,j,k);
          if(tempdirection>3) {
            if(tempdirection==4) {
              if(dataext(i,j,k)!=dataext(i+1,j,k)) {   
                dataext(i,j,k)=dataext(i+1,j,k);
                moreswaps=true;
              }
            }
            else if(tempdirection==5) {
              if(dataext(i,j,k)!=dataext(i,j+1,k)) {   
                dataext(i,j,k)=dataext(i,j+1,k);
                moreswaps=true;
              }
            }
            else if(tempdirection==6) {            
              if(dataext(i,j,k)!=dataext(i,j,k+1)) {   
                dataext(i,j,k)=dataext(i,j,k+1);
                moreswaps=true;
              }
            }        
          }        
        }
      }
      }
      } 
      
      // halo update
      if(lperiodic==true) {
        for (int i=0; i<imax; ++i) {
        for (int k=0; k<kmax; ++k) {
          if(maskext(i,jmax-2,k)==true) {
            dataext(i,0,k)=dataext(i,jmax-2,k);
          }
          if(maskext(i,1,k)==true) {
            dataext(i,jmax-1,k)=dataext(i,1,k);
          }
        }
        }
  
        for (int j=0; j<jmax; ++j) {
        for (int k=0; k<kmax; ++k) {
          if(maskext(imax-2,j,k)==true) {
            dataext(0,j,k)=dataext(imax-2,j,k);
          }
          if(maskext(1,j,k)==true) {
            dataext(imax-1,j,k)=dataext(1,j,k);
          }
        }
        }
      }
    }
    
    // renumber the cells (remember halo behavior)
    // do not assign zero
    int cloudcounter;
    cloudcounter=1;
    for (int i=1; i<imax-1; ++i) {
    for (int j=1; j<jmax-1; ++j) {
    for (int k=0; k<kmax; ++k) {
      if(maskext(i,j,k)==true) {
        if(dataext(i,j,k)==i*jmax*kmax+j*kmax+k){
           dataext(i,j,k)=cloudcounter;
           cloudcounter=cloudcounter+1;
           maskext(i,j,k)=false;
        }
      }
    }
    }
    }
    
    printf("cloudcounter +1 %d \n",cloudcounter);

    // assign identity of maximum to each "cloudy" cell
    for (int i=1; i<imax-1; ++i) {
    for (int j=1; j<jmax-1; ++j) {
    for (int k=0; k<kmax; ++k) {
      if(maskext(i,j,k)==true) {
        tempdata=dataext(i,j,k);
        if(tempdata>0){
           itarget=int(tempdata/(jmax*kmax));
           jtarget=int((tempdata-itarget*jmax*kmax)/kmax);
           ktarget=tempdata-itarget*jmax*kmax-jtarget*kmax;
           dataext(i,j,k)=dataext(itarget,jtarget,ktarget);
        }
      }
    }
    }
    }

    // halo swap to get "inner domains" at halos
    if(lperiodic==true) {
      for (int i=0; i<imax; ++i) {
      for (int k=0; k<kmax; ++k) {
     	if(dataext(i,jmax-2,k)>0) {
     	  dataext(i,0,k)=dataext(i,jmax-2,k);
     	}
     	if(dataext(i,1,k)>0) {      
     	  dataext(i,jmax-1,k)=dataext(i,1,k);
     	}
      }
      }
     
      for (int j=0; j<jmax; ++j) {
      for (int k=0; k<kmax; ++k) {
     	if(dataext(imax-2,j,k)>0) {
     	  dataext(0,j,k)=dataext(imax-2,j,k);
     	}
     	if(dataext(1,j,k)>0) {
     	  dataext(imax-1,j,k)=dataext(1,j,k);
     	}
      }
      }
    }

    // free up space (masking and direction)
    direction.free();
    maskext.free();
    
    // now the more difficult part: col identification
    // a border is identified by its position and direction
    // by clustering in 4d (i,j,k, direction), the border between cloud pairs is identified
    long borderindex;
    borderindex=(long(imax)*long(jmax)*long(kmax))/long(borderfac);
    blitz::Array<char,3> nrborders(imax,jmax,kmax); // nr of borders of current grid cell
    blitz::Array<int,3> startindex(imax,jmax,kmax); // start index (si) of current grid cell
    blitz::Array<int,1> borderfield(borderindex); // identity of the border (related to i,j,k and direction) 
    blitz::Array<indexint,1> assocfield(borderindex); // identy of border from the neighbouring grid vell

    // initialise blitz arrays
    nrborders=0; 
    startindex=0; 
    borderfield=0; 
    assocfield=0;

    int bf0,bf1,bf00; // borderfield temporary variable
    int nrb0,nrb1,nrb; // nuber of borders temporary variables
    int si0,si1; // start index temporary variables
    int colindex; // col index (to be calculated)
    indexint nrfieldtemp0,nrfieldtemp1; // cloud number temporary variables
                                        // in first loop
    indexint af0,af1; // associated fields temporary variables
    indexint nrf0,nrf1; // more cloud number temporary variables

    
    si1=1;
    // first numbering of borders on the inner domain
    for (int i=1; i<imax-1; ++i) {
    for (int j=1; j<jmax-1; ++j) {
    for (int k=0; k<kmax; ++k) {
      nrfieldtemp0=dataext(i,j,k);
      if(nrfieldtemp0>0) {
        nrb=0;
        startindex(i,j,k)=si1;
        
        nrfieldtemp1=dataext(i-1,j,k);
        if(nrfieldtemp1>0){
          if(nrfieldtemp0!=nrfieldtemp1) {
            borderfield(si1+nrb)=si1+nrb;
            assocfield(si1+nrb)=nrfieldtemp1;
            nrb=nrb+1;
          }
        }
        nrfieldtemp1=dataext(i,j-1,k);
        if(nrfieldtemp1>0){
          if(nrfieldtemp0!=nrfieldtemp1) {
            borderfield(si1+nrb)=si1+nrb;
            assocfield(si1+nrb)=nrfieldtemp1;
            nrb=nrb+1;
          }
        }
        if(k>0) {
          nrfieldtemp1=dataext(i,j,k-1);
          if(nrfieldtemp1>0){
            if(nrfieldtemp0!=nrfieldtemp1) {
              borderfield(si1+nrb)=si1+nrb;
              assocfield(si1+nrb)=nrfieldtemp1;
              nrb=nrb+1;
            }
          }
        }
        nrfieldtemp1=dataext(i+1,j,k);
        if(nrfieldtemp1>0){
          if(nrfieldtemp0!=nrfieldtemp1) {
            borderfield(si1+nrb)=si1+nrb;
            assocfield(si1+nrb)=nrfieldtemp1;
            nrb=nrb+1;
          }
        }
        nrfieldtemp1=dataext(i,j+1,k);
        if(nrfieldtemp1>0){
          if(nrfieldtemp0!=nrfieldtemp1) {
            borderfield(si1+nrb)=si1+nrb;
            assocfield(si1+nrb)=nrfieldtemp1;
            nrb=nrb+1;
          }
        }
        if(k<kmax-1) {
          nrfieldtemp1=dataext(i,j,k+1);
          if(nrfieldtemp1>0){
            if(nrfieldtemp0!=nrfieldtemp1) {
              borderfield(si1+nrb)=si1+nrb;
              assocfield(si1+nrb)=nrfieldtemp1;
              nrb=nrb+1;
            }
          }
        }
        nrborders(i,j,k)=nrb;
        si1=si1+nrb;
      }
    }
    }
    }
    printf("start index of last cell %d \n",si1);
 
    // halo update for start index and number of borders
    if(lperiodic==true) {
      for (int i=0; i<imax; ++i) {
      for (int k=0; k<kmax; ++k) {
    	if(nrborders(i,jmax-2,k)>0) {
    	  nrborders(i,0,k)=nrborders(i,jmax-2,k);
    	  startindex(i,0,k)=startindex(i,jmax-2,k);
    	}
    	if(nrborders(i,1,k)>0) {
    	  nrborders(i,jmax-1,k)=nrborders(i,1,k);
    	  startindex(i,jmax-1,k)=startindex(i,1,k);
    	}
      }
      }
    
      for (int j=0; j<jmax; ++j) {
      for (int k=0; k<kmax; ++k) {
    	 if(nrborders(imax-2,j,k)>0) {
    	   nrborders(0,j,k)=nrborders(imax-2,j,k);
    	   startindex(0,j,k)=startindex(imax-2,j,k);
    	 }
    	 if(nrborders(1,j,k)>0) {
    	   nrborders(imax-1,j,k)=nrborders(1,j,k);
    	   startindex(imax-1,j,k)=startindex(1,j,k);
    	 }
      }
      }
    }
    
    // update the borderfield by sweeping back and forth
    // innerdomain sweep, innerdomain reverse sweep
    // halo update not necesary because of code above
    // this will need further comments

    moreswaps=true;
    while(moreswaps==true) {
      moreswaps=false; //reset swap
      
      // forward swap         
      for (int i=1; i<imax-1; ++i) {
      for (int j=1; j<jmax-1; ++j) {
      for (int k=0; k<kmax; ++k) {
        nrb0=nrborders(i,j,k);
        if(nrb0>0) {
          nrf0=dataext(i,j,k);
          si0=startindex(i,j,k);      
        
          nrb1=nrborders(i-1,j,k);
          if(nrb1>0) {
            nrf1=dataext(i-1,j,k);
            for (int l=0; l<nrb0; ++l) {
              bf0=borderfield(si0+l);
              bf00=bf0;
              af0=assocfield(si0+l);
              si1=startindex(i-1,j,k);      
              for (int ll=0; ll<nrb1; ++ll) {      
                bf1=borderfield(si1+ll);
                af1=assocfield(si1+ll);
                if(bf1<bf0) {
                  // check they are on the same border
                  if((af0==nrf1 && nrf0==af1) || (af0==af1 && nrf0==nrf1)) {
                    bf0=bf1;
                    moreswaps=true;
                  }
                }
              }
              if(bf00!=bf0) {
                borderfield(si0+l)=bf0;
              }
            }
          }
  
          nrb1=nrborders(i,j-1,k);
          if(nrb1>0) {
            nrf1=dataext(i,j-1,k);
            for (int l=0; l<nrb0; ++l) {      
              bf0=borderfield(si0+l);
              bf00=bf0;
              af0=assocfield(si0+l);
              for (int ll=0; ll<nrb1; ++ll) {      
                si1=startindex(i,j-1,k);      
                bf1=borderfield(si1+ll);
                af1=assocfield(si1+ll);
                if(bf1<bf0) {
                  // check they are on the same border
                  if((af0==nrf1 && nrf0==af1) || (af0==af1 && nrf0==nrf1)) {
                    bf0=bf1;
                    moreswaps=true;
                  }
                }
              }
              if(bf00!=bf0) {
                borderfield(si0+l)=bf0;
              }
            }
          }
  
          if(k>0) {
            nrb1=nrborders(i,j,k-1);
            if(nrb1>0) {
              nrf1=dataext(i,j,k-1);
              for (int l=0; l<nrb0; ++l) {      
                bf0=borderfield(si0+l);
                bf00=bf0;
                af0=assocfield(si0+l);
                si1=startindex(i,j,k-1);      
                for (int ll=0; ll<nrb1; ++ll) {      
                  bf1=borderfield(si1+ll);
                  af1=assocfield(si1+ll);
                  if(bf1<bf0) {
                    // check they are on the same border
                    if((af0==nrf1 && nrf0==af1) || (af0==af1 && nrf0==nrf1)) {
                      bf0=bf1;
                      moreswaps=true;
                    }
                  }
                }
                if(bf00!=bf0) {
                  borderfield(si0+l)=bf0;
                }
              }
            }                              
          }
        
        }
      }
      }
      }
  
      //reverse sweep
      for (int i=imax-2; i>= 1; --i) {
      for (int j=jmax-2; j>= 1; --j) {
      for (int k=kmax-1; k>= 0; --k) {
        nrb0=nrborders(i,j,k);
  
        if(nrb0>0) {
          nrf0=dataext(i,j,k);   
          si0=startindex(i,j,k);      
                         
          nrb1=nrborders(i+1,j,k);
          if(nrb1>0) {
            nrf1=dataext(i+1,j,k);
            for (int l=0; l<nrb0; ++l) {      
              bf0=borderfield(si0+l);
              bf00=bf0;
              af0=assocfield(si0+l);
              si1=startindex(i+1,j,k);      
              for (int ll=0; ll<nrb1; ++ll) {      
                bf1=borderfield(si1+ll);
                af1=assocfield(si1+ll);
                if(bf1<bf0) {
                  // check they are on the same border
                  if((af0==nrf1 && nrf0==af1) || (af0==af1 && nrf0==nrf1)) {
                    bf0=bf1;
                    moreswaps=true;
                  }
                }
              }
              if(bf00!=bf0) {
                borderfield(si0+l)=bf0;
              }
            }
          }
  
          nrb1=nrborders(i,j+1,k);
          if(nrb1>0) {
            nrf1=dataext(i,j+1,k);
            for (int l=0; l<nrb0; ++l) {      
              bf0=borderfield(si0+l);
              bf00=bf0;
              af0=assocfield(si0+l);
              si1=startindex(i,j+1,k);      
              for (int ll=0; ll<nrb1; ++ll) {      
                bf1=borderfield(si1+ll);
                af1=assocfield(si1+ll);
                if(bf1<bf0) {
                  // check they are on the same border
                  if((af0==nrf1 && nrf0==af1) || (af0==af1 && nrf0==nrf1)) {
                    bf0=bf1;
                    moreswaps=true;
                  }
                }
              }
              if(bf00!=bf0) {
                borderfield(si0+l)=bf0;
              }
            }
          }
  
          if(k<kmax-1) {
            nrb1=nrborders(i,j,k+1);
            if(nrb1>0) {
              nrf1=dataext(i,j,k+1);
              for (int l=0; l<nrb0; ++l) {      
                bf0=borderfield(si0+l);
                bf00=bf0;
                af0=assocfield(si0+l);
                si1=startindex(i,j,k+1);      
                for (int ll=0; ll<nrb1; ++ll) {      
                  bf1=borderfield(si1+ll);
                  af1=assocfield(si1+ll);
                  if(bf1<bf0) {
                    // check they are on the same border
                    if((af0==nrf1 && nrf0==af1) || (af0==af1 && nrf0==nrf1)) {
                      bf0=bf1;
                      moreswaps=true;
                    }
                  }
                }
                if(bf00!=bf0) {
                  borderfield(si0+l)=bf0;
                }
              }
            }                              
          }
        }      
      }
      }
      }
    }
    
    // renumber the borders
    counter=1; //avoid zero
    for (int i=1; i<imax-1; ++i) {
    for (int j=1; j<jmax-1; ++j) {
    for (int k=0; k<kmax; ++k) {
      nrb=nrborders(i,j,k);
      if(nrb>0) {
        si1=startindex(i,j,k);
        for (int l=0; l<nrb; ++l) {                
          if(borderfield(si1+l)==si1+l) {
            borderfield(si1+l)=-counter;
            counter=counter+1;
          }
        }
      }
    }
    }
    }

    printf("borderfield counter +1 %d \n",counter);
           
    for (int i=1; i<imax-1; ++i) {
    for (int j=1; j<jmax-1; ++j) {
    for (int k=0; k<kmax; ++k) {
      nrb=nrborders(i,j,k);
      if(nrb>0) {
        si1=startindex(i,j,k);
        for (int l=0; l<nrb; ++l) {      
          if(borderfield(si1+l)>0) {
            ltarget=borderfield(si1+l);
            borderfield(si1+l)=borderfield(ltarget);
          }
        }
      }
    }
    }
    }
    
    for (int i=1; i<imax-1; ++i) {
    for (int j=1; j<jmax-1; ++j) {
    for (int k=0; k<kmax; ++k) {
      nrb=nrborders(i,j,k);
      if(nrb>0) {
        si1=startindex(i,j,k);
        for (int l=0; l<nrb; ++l) {      
          if(borderfield(si1+l)<0){
            borderfield(si1+l)=-borderfield(si1+l);
          }
        }
      }
    }
    }
    }

    nrborders.free();
    startindex.free();

    // reopen the field data, now memory has become available again
    blitz::Array<short,3> fieldext2(imax,jmax,kmax);

    /* Open the file. NC_NOWRITE tells netCDF we want read-only access
     * to the file.*/
    if ((retval = nc_open("mask_field.nc", NC_NOWRITE, &ncid3)))
       ERR(retval);

   err=blitzncshort(ncid3,"fieldext",fieldext2);

    /* Close the file, freeing all resources. */
    if ((retval = nc_close(ncid3)))
       ERR(retval);
    
    // initialise array which holds data on cols
    // their heights and the associated "clouds"
    long colmax;
    colmax=(long(imax)*long(jmax)*long(kmax))/long(colfac);            
    blitz::Array<int,2> coldata(colmax,4);
    
    coldata=-2147483648; // lowest possible value

    // repeat calculation of si1 and nrb on the fly here to save memory
    si1=1;
    for (int i=1; i<imax-1; ++i) {
    for (int j=1; j<jmax-1; ++j) {
    for (int k=0; k<kmax; ++k) {
      nrfieldtemp0=dataext(i,j,k);
      if(nrfieldtemp0>0) {
        nrb=0;        
        nrfieldtemp1=dataext(i-1,j,k);
        if(nrfieldtemp1>0){
          if(nrfieldtemp0!=nrfieldtemp1) {
            nrb=nrb+1;
          }
        }
        nrfieldtemp1=dataext(i,j-1,k);
        if(nrfieldtemp1>0){
          if(nrfieldtemp0!=nrfieldtemp1) {
            nrb=nrb+1;
          }
        }
        if(k>0) {
          nrfieldtemp1=dataext(i,j,k-1);
          if(nrfieldtemp1>0){
            if(nrfieldtemp0!=nrfieldtemp1) {
              nrb=nrb+1;
            }
          }
        }
        nrfieldtemp1=dataext(i+1,j,k);
        if(nrfieldtemp1>0){
          if(nrfieldtemp0!=nrfieldtemp1) {
            nrb=nrb+1;
          }
        }
        nrfieldtemp1=dataext(i,j+1,k);
        if(nrfieldtemp1>0){
          if(nrfieldtemp0!=nrfieldtemp1) {
            nrb=nrb+1;
          }
        }
        if(k<kmax-1) {
          nrfieldtemp1=dataext(i,j,k+1);
          if(nrfieldtemp1>0){
            if(nrfieldtemp0!=nrfieldtemp1) {
              nrb=nrb+1;
            }
          }
        }
        if(nrb>0) {
          for (int l=0; l<nrb; ++l) {                
            colindex=borderfield(si1+l);    
            coldata(colindex,2)=dataext(i,j,k);
            coldata(colindex,3)=assocfield(si1+l);
            if(dataext(i,j,k)<assocfield(si1+l)) {
              if(fieldext2(i,j,k)>coldata(colindex,0)) {
                coldata(colindex,0)=fieldext2(i,j,k);
              }
            }
            else {
              if(fieldext2(i,j,k)>coldata(colindex,1)) {
                coldata(colindex,1)=fieldext2(i,j,k);
              }
            }
          } 
        }
        si1=si1+nrb;
      }
    }
    }
    }
    printf("start index of last cell %d \n",si1);
    
    borderfield.free();
    assocfield.free();
    
    for (int colindex=0; colindex<counter; ++colindex) {
      if(coldata(colindex,0)>coldata(colindex,1)) {
        coldata(colindex,0)=coldata(colindex,1);
      }
    }
    
    // calculate minima and maxima associated with different "clouds"
    long blobindex;
    blobindex=(long(imax)*long(jmax)*long(kmax))/long(blobfac);
    blitz::Array<short,1> blobmins(blobindex);
    blitz::Array<short,1> blobmaxs(blobindex);
    
    blobmins=32767;
    blobmaxs=-32768;
    
    indexint nrfind;
       
    for (int i=1; i<imax-1; ++i) {
    for (int j=1; j<jmax-1; ++j) {
    for (int k=0; k<kmax; ++k) {
       nrfind=dataext(i,j,k);
       if(nrfind>0) {
          if(fieldext2(i,j,k)>blobmaxs(nrfind)) {
             blobmaxs(nrfind)=fieldext2(i,j,k);
          }
          if(fieldext2(i,j,k)<blobmins(nrfind)) {
             blobmins(nrfind)=fieldext2(i,j,k);
          }
       }
    }
    }
    } 

    fieldext2.free();

    int cld1,cld2;
    
    // SORT COLDATA BY FIRST COLUMN
    quick_sort(coldata, 1, counter-2);
    
    // array which points to the "parent cloud"             
    blitz::Array<int,1> targetcld(blobindex);

    targetcld=0;
    for (long i=0; i<blobindex; ++i) {
      targetcld(i)=i;
    } 
      
    int targetcld1,targetcld2,targetcld11,targetcld22;
    short blobmintemp;  
    int m1,m2,m3;
    float m11,m22,m33,col,colratio;
    
    for (int i=counter-1; i>=1; --i) {
      cld1=coldata(i,2);
      cld2=coldata(i,3);
      targetcld1=targetcld(cld1);
      while(targetcld1!=targetcld(targetcld1)) {
        targetcld1=targetcld(targetcld1);
      }
      targetcld2=targetcld(cld2);
      while(targetcld2!=targetcld(targetcld2)) {
        targetcld2=targetcld(targetcld2);
      }

      // colratio determines if clouds get merged
      m1=std::max(blobmaxs(targetcld1),blobmaxs(targetcld2));
      m11=inv_scaling_parameter*(m1*m1*(2*(m1>0)-1)); 
      m2=std::min(blobmins(targetcld1),blobmins(targetcld2));
      m22=inv_scaling_parameter*(m2*m2*(2*(m2>0)-1));     
      m3=std::min(blobmaxs(targetcld1),blobmaxs(targetcld2));
      m33=inv_scaling_parameter*(m3*m3*(2*(m3>0)-1)); 
      col=inv_scaling_parameter*(coldata(i,0)*coldata(i,0))*(2*(coldata(i,0)>0)-1);
      colratio=(m33-col)/(m33-m22); // (lowest peak-col)/(lowest peak-lowest point)

      // first mergers determined here
      if (colratio<=mincolratio) {
        if(blobmaxs(targetcld1)<blobmaxs(targetcld2)) { // merge into higher peak
          targetcld1=targetcld(cld1);
          blobmintemp=blobmins(cld1);
          targetcld(cld1)=targetcld2;
          while(targetcld1!=targetcld2) {
             targetcld11=targetcld(targetcld1);
             blobmintemp=std::min(blobmins(targetcld1),blobmintemp);
             targetcld(targetcld1)=targetcld2;
             targetcld1=targetcld11;
          }
          blobmins(targetcld2)=std::min(blobmintemp,blobmins(targetcld2));
        }
        else{
          targetcld2=targetcld(cld2);
          blobmintemp=blobmins(cld2);
          targetcld(cld2)=targetcld1;
          while(targetcld2!=targetcld1) {
             targetcld22=targetcld(targetcld2);
             blobmintemp=std::min(blobmins(targetcld2),blobmintemp);
             targetcld(targetcld2)=targetcld1;
             targetcld2=targetcld22;
          }
          blobmins(targetcld1)=std::min(blobmintemp,blobmins(targetcld1));
        }
      }
    }

    // update all remaining targets
    // successive merges
    int tcld;
    moreswaps=true;
    while(moreswaps==true) {
      moreswaps=false;
      for (int i=cloudcounter-1; i>=1; --i) {
        tcld=targetcld(i);
        if(tcld!=targetcld(tcld)){
          targetcld(i)=targetcld(tcld);
          moreswaps=true;
        }
      }
      if(moreswaps==true) {
        moreswaps=false;
        for (int i=1; i<cloudcounter; ++i) {
          tcld=targetcld(i);
          if(tcld!=targetcld(tcld)){
            targetcld(i)=targetcld(tcld);
            moreswaps=true;
          }          
        }        
      }
    }
                      
    coldata.free();
    blobmins.free();
    blobmaxs.free();

    int refercloud;
    int encounter;
    
    // let the target cloud replace the cloud number   
    for (int i=1; i<imax-1; ++i) {
    for (int j=1; j<jmax-1; ++j) {
    for (int k=0; k<kmax; ++k) {
       if(dataext(i,j,k)>0) {
          tcld=dataext(i,j,k);
          tcld=targetcld(tcld);
          dataext(i,j,k)=tcld;
       }
    }
    }
    }
    
    targetcld.free();
    blitz::Array<int,1> encounterorder(blobindex);
    encounterorder=0;

    encounter=1;
    // renumbering procedure
    for (int i=1; i<imax-1; ++i) {
    for (int j=1; j<jmax-1; ++j) {
    for (int k=0; k<kmax; ++k) {
       if(dataext(i,j,k)>0) {
          tcld=dataext(i,j,k);
          if (encounterorder(tcld)==0) {
            encounterorder(tcld)=encounter;
            dataext(i,j,k)=encounter;
            encounter=encounter+1;
          }
          else {
            dataext(i,j,k)=encounterorder(tcld);
          }              
       }
    }
    }
    }
    printf("number of clouds after merge +1 %d \n",encounter);

    // set halo values to zero, so we can check against input easily
    for (int i=0; i<imax; ++i) {
    for (int k=0; k<kmax; ++k) {
      dataext(i,0,k)=0;
      dataext(i,jmax-1,k)=0;
    }
    }

    for (int j=0; j<jmax; ++j) {
    for (int k=0; k<kmax; ++k) {
        dataext(0,j,k)=0;
        dataext(imax-1,j,k)=0;
    }
    }

    int x_dimid, y_dimid, z_dimid;
    int dimids[3];

    encounterorder.free();

   /* Always check the return code of every netCDF function call. In
    * this example program, any retval which is not equal to NC_NOERR
    * (0) will cause the program to print an error message and exit
    * with a non-zero return code. */

   /* Create the file. The NC_CLOBBER parameter tells netCDF to
    * overwrite this file, if it already exists.*/
   if ((retval = nc_create("output.nc",NC_CLOBBER|NC_NETCDF4, &ncid4)))
      ERR(retval);

   /* Define the dimensions. NetCDF will hand back an ID for each. */
   if ((retval = nc_def_dim(ncid4, "x", imax, &x_dimid)))
      ERR(retval);
   if ((retval = nc_def_dim(ncid4, "y", jmax, &y_dimid)))
      ERR(retval);
   if ((retval = nc_def_dim(ncid4, "z", kmax, &z_dimid)))
      ERR(retval);
      
   /* The dimids array is used to pass the IDs of the dimensions of
    * the variable. */
   dimids[0] = x_dimid;
   dimids[1] = y_dimid;
   dimids[2] = z_dimid;

   /* Define the variable. The type of the variable in this case is
    * NC_INT (4-byte integer). */
   if ((retval = nc_def_var(ncid4, "data", NC_UINT, 3, 
			    dimids, &varid)))
      ERR(retval);

   /* End define mode. This tells netCDF we are done defining
    * metadata. */
   if ((retval = nc_enddef(ncid4)))
      ERR(retval);

   /* Write the pretend data to the file. Although netCDF supports
    * reading and writing subsets of data, in this case we write all
    * the data in one operation. */
   
   if ((retval = nc_put_var_uint(ncid4, varid, dataext.data())))
      ERR(retval);

   /* Close the file. This frees up any internal netCDF resources
    * associated with the file, and flushes any buffers. */
   if ((retval = nc_close(ncid4)))
      ERR(retval);

   printf("*** SUCCESS writing example file output.nc!\n");
   return 0;
   
  }
