#include "common.h"
#include <netcdf.h>
#include "blitz/array.h"

#define ERRCODE 2
#define ERR(e) {printf("Error: %s\n", nc_strerror(e)); exit(ERRCODE);}


// fill a blitz array with "short integer" netcdf data
void blitzncshort(int ncId, const char *varName,
blitz::Array<short,3> & tab)
{
  int status, varId, nbDims;

  printf("Loading `%s`...\n", varName);

  status = nc_inq_varid(ncId, varName, &varId);
  if (status != NC_NOERR)
    std::cout << varName << " is absent" << std::endl ;
  else
  {
    int *dimIds;
    size_t start[] = {0,0,0};
    size_t dim_lengths[3];

    status = nc_inq_varndims(ncId, varId, &nbDims);

    if(nbDims != 3)
    {
      fprintf(stderr,"Error in data, the variable %s should have 3 dimensions\n", varName);
      exit(-1);
    }

    dimIds = new int[nbDims*sizeof(*dimIds)];
    status = nc_inq_vardimid(ncId, varId, dimIds);

    if (status != NC_NOERR) {
      ERR(status);
    }

    for (int n=0; n<3; n++) {
      nc_inq_dimlen(ncId, dimIds[n], &dim_lengths[n]);

      if (status != NC_NOERR) {
        ERR(status);
      }
    }

    delete [] dimIds;

    tab.resize(dim_lengths[0], dim_lengths[1], dim_lengths[2]);

    status = nc_get_vara_short(ncId, varId, start, dim_lengths, tab.data());
    if (status != NC_NOERR) {
      ERR(status);
    }
  }
}


int write_netcdf(blitz::Array<indexint,3> dataext) {
  blitz::TinyVector<int,3> shape = dataext.shape();
  int imax = shape[0];
  int jmax = shape[1];
  int kmax = shape[2];


    int retval; // error terurn value
    int ncid; //netcdf ids (decided not to reuse these)
    int varid; //variable ids (decided not to reuse these)

    int x_dimid, y_dimid, z_dimid;
    int dimids[3];

   /* Always check the return code of every netCDF function call. In
    * this example program, any retval which is not equal to NC_NOERR
    * (0) will cause the program to print an error message and exit
    * with a non-zero return code. */

   /* Create the file. The NC_CLOBBER parameter tells netCDF to
    * overwrite this file, if it already exists.*/
   if ((retval = nc_create("output.nc",NC_CLOBBER|NC_NETCDF4, &ncid)))
      ERR(retval);

   /* Define the dimensions. NetCDF will hand back an ID for each. */
   if ((retval = nc_def_dim(ncid, "x", imax, &x_dimid)))
      ERR(retval);
   if ((retval = nc_def_dim(ncid, "y", jmax, &y_dimid)))
      ERR(retval);
   if ((retval = nc_def_dim(ncid, "z", kmax, &z_dimid)))
      ERR(retval);
      
   /* The dimids array is used to pass the IDs of the dimensions of
    * the variable. */
   dimids[0] = x_dimid;
   dimids[1] = y_dimid;
   dimids[2] = z_dimid;

   /* Define the variable. The type of the variable in this case is
    * NC_INT (4-byte integer). */
   if ((retval = nc_def_var(ncid, "data", NC_UINT, 3, 
          dimids, &varid)))
      ERR(retval);

   /* End define mode. This tells netCDF we are done defining
    * metadata. */
   if ((retval = nc_enddef(ncid)))
      ERR(retval);

   /* Write the pretend data to the file. Although netCDF supports
    * reading and writing subsets of data, in this case we write all
    * the data in one operation. */
   
   if ((retval = nc_put_var_uint(ncid, varid, dataext.data())))
      ERR(retval);

   /* Close the file. This frees up any internal netCDF resources
    * associated with the file, and flushes any buffers. */
   if ((retval = nc_close(ncid)))
      ERR(retval);

   printf("*** SUCCESS writing example file output.nc!\n");
   return 0;
}


void load_mask(blitz::Array<bool,3> &maskext) {
  int retval; // error return value
  int ncid; //netcdf id

  indexint counter; // counts cells with certain properties (used for sanity checks)

  maskext=false;

  // array that holds "cloud" mask as short (temporarily)
  blitz::Array<short,3> maskshort;

  printf("Opening `mask_field.nc`\n");

  /* Open the file. NC_NOWRITE tells netCDF we want read-only access
   * to the file.*/
  if ((retval = nc_open("mask_field.nc", NC_NOWRITE, &ncid)))
    ERR(retval);

  blitzncshort(ncid,"maskext",maskshort);
  maskext.resize(maskshort.shape());

  /* Close the file, freeing all resources. */
  if ((retval = nc_close(ncid)))
    ERR(retval);

  blitz::TinyVector<int,3> shape = maskshort.shape();
  int imax = shape[0];
  int jmax = shape[1];
  int kmax = shape[2];

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
}

void load_field(blitz::Array<short,3> &fieldext) {
  int retval; // error return value
  int ncid; //netcdf id

  fieldext=0;

  printf("Opening `mask_field.nc`\n");

  /* Open the file. NC_NOWRITE tells netCDF we want read-only access
   * to the file.*/
  if ((retval = nc_open("mask_field.nc", NC_NOWRITE, &ncid)))
    ERR(retval);

  blitzncshort(ncid,"fieldext",fieldext);

  /* Close the file, freeing all resources. */
  if ((retval = nc_close(ncid)))
    ERR(retval);
}
