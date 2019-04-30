#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <stdint.h>
#include <exception>
#include "blitz/array.h"

#include "common.h"
#include "blitz_sort.h"
#include "file_io.h"


/* Exception used when scaling if input is outside of valid range
 */
struct InvalidInputException : public std::exception {
  double m_smin, m_smax;

  InvalidInputException(double smin, double smax) : m_smin(smin), m_smax(smax) {
  }

  const char* what() const throw() {
    std::string s = "Input range [" + std::to_string(m_smin) + ":"
                    + std::to_string(m_smax) + "] is outside of expected scaling "
                    " range [" + std::to_string(min_scalar_value) + ":"
                    + std::to_string(max_scalar_value) + "]."
                    " Please change `max_scalar_value` in common.h and recompile";

    // need a char pointer... WTF, thank god for stack overflow (https://stackoverflow.com/a/16502000)
    return &s[0u];
  }
};

/* Apply scaling to short ints for optimisation. Checks bounds on input values
 * are within expected range so overflows don't occour
 */
void scale_field(
  blitz::Array<double,3> &field_input,
  blitz::Array<short,3> &field_scaled
) {
  double s_min = min(field_input);
  double s_max = max(field_input);

  if (s_min < min_scalar_value || s_max > max_scalar_value) {
    throw InvalidInputException(s_min, s_max);
  }
  short sign = 0;
  blitz::TinyVector<int,3> shape = field_input.shape();

  for (int i=0; i<shape[0]; i++) {
    for (int j=0; j<shape[1]; j++) {
      for (int k=0; k<shape[2]; k++) {
        if (field_input(i,j,k) > 0.0) {
          sign = 1;
        }
        else {
          sign = -1;
        }
        field_scaled(i,j,k) = sign*(short)(scaling_parameter*std::sqrt(std::abs(field_input(i,j,k))));
      }
    }
  }
}

float descale(short i) {
  return inv_scaling_parameter*((float)(i*i)*(2.0*(float)(i>0)-1));
}


/** halo swap to get "inner domain values" into halos
 */
template <typename Type>
void halo_swap(blitz::Array<Type,3> &field) {
  blitz::TinyVector<int,3> shape = field.shape();
  int imax = shape[0];
  int jmax = shape[1];
  int kmax = shape[2];

  for (int i=0; i<imax; ++i) {
    for (int k=0; k<kmax; ++k) {
      if(field(i,jmax-2,k)==true) {
        field(i,0,k)=field(i,jmax-2,k);
      }
      if(field(i,1,k)==true) {
        field(i,jmax-1,k)=field(i,1,k);
      }
    }
  }

  for (int j=0; j<jmax; ++j) {
    for (int k=0; k<kmax; ++k) {
      if(field(imax-2,j,k)==true) {
        field(0,j,k)=field(imax-2,j,k);
      }
      if(field(1,j,k)==true) {
        field(imax-1,j,k)=field(1,j,k);
      }
    }
  }
}


/** detection of steepest gradient direction in inner domain
 * 0 corresponds to local maxima
 */
void identify_steepest_ascent(
  const blitz::Array<bool, 3> &maskext,
  const blitz::Array<short,3> &fieldext,
  blitz::Array<char, 3> direction
  ) {
  blitz::TinyVector<int,3> shape = fieldext.shape();
  int imax = shape[0];
  int jmax = shape[1];
  int kmax = shape[2];

  char tempdirection; // direction of steepest ascent temp-vars are used in inner loops
  int tempmax,tempfield; // used in inner loop

  // the sum of directions can be used to do checks on a 
  // domain that repeats itself in the horizontal directions
  int sumdir; 
  sumdir=0;

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
}


/* assign points to local maximum   
 * order: inner domain forward sweep, halo update, inner domain reverse sweep, haloupdate
 *
 * XXX: `maskext` is modified so that it is set false at local maxima
 * 
 */
void assign_to_local_maxima(
  blitz::Array<bool,    3> &maskext,
  blitz::Array<indexint,3> &dataext,
  const blitz::Array<char,3> direction,
  int *cloudcounter
) {
  blitz::TinyVector<int,3> shape = maskext.shape();
  int imax = shape[0];
  int jmax = shape[1];
  int kmax = shape[2];

  int itarget,jtarget,ktarget; //calculating back target i,j,k,l that corresponds to maxima
  char tempdirection; // direction of steepest ascent temp-vars are used in inner loops

  indexint tempdata; // used in inner loop
  bool moreswaps; // used to check if all identities have been assigned, or further checking necessary

  printf("Assigning in-cloud points to local maxima\n");

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
      halo_swap<indexint>(dataext);
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
      halo_swap<indexint>(dataext);
    }
  }

  // renumber the cells (remember halo behavior)
  // do not assign zero
  (*cloudcounter)=1;
  for (int i=1; i<imax-1; ++i) {
    for (int j=1; j<jmax-1; ++j) {
      for (int k=0; k<kmax; ++k) {
        if(maskext(i,j,k)==true) {
          if(dataext(i,j,k)==i*jmax*kmax+j*kmax+k){
            dataext(i,j,k)=*cloudcounter;
            (*cloudcounter)=(*cloudcounter)+1;
            maskext(i,j,k)=false;
          }
        }
      }
    }
  }

  printf("cloudcounter +1 %d \n", *cloudcounter);

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
    halo_swap(dataext);
  }
}

void init_output(
  blitz::Array<bool,    3> &maskext,
  blitz::Array<indexint,3> &dataext
) {
  blitz::TinyVector<int,3> shape = maskext.shape();
  int imax = shape[0];
  int jmax = shape[1];
  int kmax = shape[2];

  dataext=0;

  // initial numbering
  for (int i=1; i<imax-1; ++i) {
    for (int j=1; j<jmax-1; ++j) {
      for (int k=0; k<kmax; ++k) {
        if(maskext(i,j,k)==true) {
          dataext(i,j,k)=i*jmax*kmax+j*kmax+k;
        }
      }
    }
  }
}


void identify_borders(
    // count number of borders on each cell
    // a border being a cell interface with a different cloud 
    // borders get assigned a border number (borderfield) and the associated
    // cloud number at the other side of the border (assocfield)
    const blitz::Array<indexint,3> &dataext,
    blitz::Array<int,1> borderfield, 
    blitz::Array<indexint,1> assocfield,
    indexint *borderfield_counter,
    const int cloudcounter
) {
  blitz::TinyVector<int,3> shape = dataext.shape();
  int imax = shape[0];
  int jmax = shape[1];
  int kmax = shape[2];

    int ltarget; //calculating back target i,j,k,l that corresponds to maxima
    bool moreswaps; // used to check if all identities have been assigned, or further checking necessary

    blitz::Array<char,3> nrborders(imax,jmax,kmax); // nr of borders of current grid cell
    blitz::Array<int,3> startindex(imax,jmax,kmax); // start index (si) of current grid cell

    // initialise blitz arrays
    nrborders=0; 
    startindex=0; 
    borderfield=0; 
    assocfield=0;

    int bf0,bf1,bf00; // borderfield temporary variable
    int nrb0,nrb1,nrb; // nuber of borders temporary variables
    int si0,si1; // start index temporary variables
    indexint nrfieldtemp0,nrfieldtemp1; // cloud number temporary variables
                                        // in first loop
    indexint af0,af1; // associated fields temporary variables
    indexint nrf0,nrf1; // more cloud number temporary variables

    
    si1=1;
    // first numbering of borders on the inner domain
    // si1 keeps track of the total number of border cells in the domain
    for (int i=1; i<imax-1; ++i) {
    for (int j=1; j<jmax-1; ++j) {
    for (int k=0; k<kmax; ++k) {
      nrfieldtemp0=dataext(i,j,k);
      if(nrfieldtemp0>0) {
        nrb=0;
        startindex(i,j,k)=si1;
        
        nrfieldtemp1=dataext(i-1,j,k);
        // check if neighboring cell is a cloud
        if(nrfieldtemp1>0){
          //check it is a different cloud
          if(nrfieldtemp0!=nrfieldtemp1) {
            //if so, add to the borderfield/assocfield arrays
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
        // add the number of borders of this cell to the start index
        si1=si1+nrb;

#ifdef DEBUG_CHECK_LIMITS
            if (si1+6 >= borderfield.size()) {
              printf("Error: `borderfield` array too small, decrease `borderfac` to increase size\n");
              exit(-1);
            }
#endif
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
    // halo update not necessary because of code above
    // this will need further comments
    
    // this is similar to the original identification
    // of cloud numbers
    // but now adding an index for each edge of the cell
    // (keeping in mind actual borders are sparse)
    // each cell is assigned the lowest value of bf on
    // the edges of its
    // neighbouring cells. Criteria:
    // 1) same cloud number and same associated cloud
    // 2) cloud number and associated cloud number swapped
    //    i.e. the different side of the same edge.
    
    // Possibly we should check for
    // actual adjacency
    // i.e. whether the edges actually touch
    // to identify actual borders
    
    // However, this should not to affect results 
    // when only the highest col matters
    
    // Example below
    //  1 2 2 2 1
    //  1 1 2 1 1
    //  1 - - - 1
    //  1 1 1 1 1
    //  Both sides of the 2 on the second row may be identified as same border here
    
    moreswaps=true;
    while(moreswaps==true) {
      moreswaps=false; //reset swap
      
      // forward swap         
      for (int i=1; i<imax-1; ++i) {
      for (int j=1; j<jmax-1; ++j) {
      for (int k=0; k<kmax; ++k) {
        nrb0=nrborders(i,j,k);
        if(nrb0>0) {
          nrf0=dataext(i,j,k); // cloud number of current cell
          si0=startindex(i,j,k); // first index of "borderfield" associated with
                                 // this cell       
        
          nrb1=nrborders(i-1,j,k); 
          
          // for each neighbouring cell
          // check if the two cells share neighbour pairs for
          // which condition 1 or 2 above is met.
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
              // only update in 3D array if really needed
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
    // use negative indices here
    (*borderfield_counter)=1; //avoid zero
    for (int i=1; i<imax-1; ++i) {
    for (int j=1; j<jmax-1; ++j) {
    for (int k=0; k<kmax; ++k) {
      nrb=nrborders(i,j,k);
      if(nrb>0) {
        si1=startindex(i,j,k);
        for (int l=0; l<nrb; ++l) {                
          if(borderfield(si1+l)==si1+l) {
            borderfield(si1+l)=-(*borderfield_counter);
            (*borderfield_counter)=(*borderfield_counter)+1;
          }
        }
      }
    }
    }
    }

    printf("borderfield borderfield_counter +1 %d \n", *borderfield_counter);
    
    // update locations that did not keep their index as well
    // giving everything negative indices      
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
    
    // swap sign of indices back to positive
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
}

void find_cols_on_borders(
    const blitz::Array<short,3> &fieldext,
    const blitz::Array<indexint,3> &dataext,
    const blitz::Array<int,1> &borderfield,
    const blitz::Array<indexint,1> &assocfield,
    blitz::Array<int,2> &coldata,
    indexint borderfield_counter
) {
  blitz::TinyVector<int,3> shape = dataext.shape();
  int imax = shape[0];
  int jmax = shape[1];
  int kmax = shape[2];

    int colindex; // col index (to be calculated)

    int nrb; // nuber of borders temporary variables
    int si1; // start index temporary variables
    indexint nrfieldtemp0,nrfieldtemp1; // cloud number temporary variables
                                        // in first loop

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
            // look at the  "lowest side" of the border
            if(dataext(i,j,k)<assocfield(si1+l)) {
              // find the highest of these "lowest side" points
              if(fieldext(i,j,k)>coldata(colindex,0)) {
                coldata(colindex,0)=fieldext(i,j,k);
              }
            }
            else {
              if(fieldext(i,j,k)>coldata(colindex,1)) {
                coldata(colindex,1)=fieldext(i,j,k);
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
    
    for (int colindex=0; colindex<borderfield_counter; ++colindex) {
      if(coldata(colindex,0)>coldata(colindex,1)) {
        coldata(colindex,0)=coldata(colindex,1);
      }
    }
    
    // SORT COLDATA BY FIRST COLUMN
    quick_sort(coldata, 1, borderfield_counter-2);
}

void merge_along_cols(
  const blitz::Array<short,3> &fieldext,
  blitz::Array<indexint,3> &dataext,
  const blitz::Array<int,2> &coldata,
  const indexint cloudcounter,
  const indexint borderfield_counter
) {
  blitz::TinyVector<int,3> shape = dataext.shape();
  int imax = shape[0];
  int jmax = shape[1];
  int kmax = shape[2];

    bool moreswaps; // used to check if all identities have been assigned, or further checking necessary

    int cld1,cld2;
    
    // array which points to the "parent cloud"             
    long blobindex = (long(fieldext.size()))/long(blobfac);
    blitz::Array<int,1> targetcld(blobindex);
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
            if(fieldext(i,j,k)>blobmaxs(nrfind)) {
              blobmaxs(nrfind)=fieldext(i,j,k); // highest point on a "blob"
            }
            if(fieldext(i,j,k)<blobmins(nrfind)) {
              blobmins(nrfind)=fieldext(i,j,k); // lowest point on a "blob"
            }
          }
        }
      }
    } 

    targetcld=0;
    for (long i=0; i<blobindex; ++i) {
      targetcld(i)=i;
    } 
      
    int targetcld1,targetcld2,targetcld11,targetcld22;
    short blobmintemp;  
    int m1,m2,m3;
    float m11,m22,m33,col,colratio;
    
    for (int i=borderfield_counter-1; i>=1; --i) {
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
      // note the non-linear transformation used here
      m1=std::max(blobmaxs(targetcld1),blobmaxs(targetcld2));
      //m11=inv_scaling_parameter*(m1*m1*(2*(m1>0)-1));
      m11 = descale(m1);
      m2=std::min(blobmins(targetcld1),blobmins(targetcld2));
      //m22=inv_scaling_parameter*(m2*m2*(2*(m2>0)-1));
      m22 = descale(m2);
      m3=std::min(blobmaxs(targetcld1),blobmaxs(targetcld2));
      //m33=inv_scaling_parameter*(m3*m3*(2*(m3>0)-1));
      m33 = descale(m3);
      //col=inv_scaling_parameter*(coldata(i,0)*coldata(i,0))*(2*(coldata(i,0)>0)-1);
      col = descale(coldata(i,0));
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
    // reverse order after each swap
    // the targetcld array contains the eventual target
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
                      
    blobmins.free();
    blobmaxs.free();

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

    encounterorder.free();
   
  }


void merge_smaller_peaks(
  blitz::Array<indexint,3> &dataext,
  blitz::Array<short,3> &fieldext,
  const int cloudcounter
) {

  long borderindex = long(dataext.size())/long(borderfac);

  blitz::Array<int,1> borderfield(borderindex); // identity of the border (related to i,j,k and direction) 
  blitz::Array<indexint,1> assocfield(borderindex); // identy of border from the neighbouring grid vell
  indexint borderfield_counter = 1; // counts cells with certain properties (used for sanity checks)

  identify_borders(dataext, borderfield, assocfield, &borderfield_counter, cloudcounter);

  // initialise array which holds data on cols
  // their heights and the associated "clouds"
  long colmax = (long(dataext.size()))/long(colfac);
  blitz::Array<int,2> coldata(colmax,4);
  coldata=-2147483648; // lowest possible value

#ifdef DEBUG_CHECK_LIMITS
    for (int l=0; l<borderindex; ++l) {                
      if (borderfield(l) >= colmax) {
        printf("Error: `colmax` array too small, decrease `colfac` to increase size\n");
        exit(-1);
      }
    }
#endif

  if (fieldext.size() == 0) {
    // reopen the field data, now memory has become available again
    load_field(fieldext);
  }

  find_cols_on_borders(fieldext, dataext, borderfield, assocfield, coldata, borderfield_counter);
  borderfield.free();
  assocfield.free();


  merge_along_cols(fieldext, dataext, coldata, cloudcounter, borderfield_counter);
  fieldext.free();
  coldata.free();
}

/*
 * Common interface for running from command line and calling from python, 
 * if fieldext and maskext are zero length then it will be attempted to load data from file.
 */
void find_objects(
  blitz::Array<short,3> &fieldext,
  blitz::Array<bool, 3> &maskext,
  blitz::Array<indexint,3> &dataext
) {

  bool dynamically_load_fields = (fieldext.size() == 0 || maskext.size() == 0);

  if (dynamically_load_fields) {
    load_field(fieldext);
    load_mask(maskext);
    dataext.resize(fieldext.shape());
  }
  init_output(maskext, dataext);

  if(lperiodic==true) {
    halo_swap<bool>(maskext);
    halo_swap<indexint>(dataext);
    halo_swap<short>(fieldext);
  }

  // array that holds direction of steepest ascent
  blitz::Array<char,3> direction(fieldext.shape());
  direction=0;
  identify_steepest_ascent(maskext, fieldext, direction);

  if (dynamically_load_fields) {
    fieldext.free();
  }

  int cloudcounter;
  assign_to_local_maxima(maskext, dataext, direction, &cloudcounter);

  if (dynamically_load_fields) {
    maskext.free();
  }

  // free up space (masking and direction)
  direction.free();

  merge_smaller_peaks(dataext, fieldext, cloudcounter);
}
