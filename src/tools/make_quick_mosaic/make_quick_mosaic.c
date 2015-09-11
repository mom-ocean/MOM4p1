/*

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                                   !!
!!                   GNU General Public License                      !!
!!                                                                   !!
!! This file is part of the Flexible Modeling System (FMS).          !!
!!                                                                   !!
!! FMS is free software; you can redistribute it and/or modify       !!
!! it and are expected to follow the terms of the GNU General Public !!
!! License as published by the Free Software Foundation.             !!
!!                                                                   !!
!! FMS is distributed in the hope that it will be useful,            !!
!! but WITHOUT ANY WARRANTY; without even the implied warranty of    !!
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the     !!
!! GNU General Public License for more details.                      !!
!!                                                                   !!
!! You should have received a copy of the GNU General Public License !!
!! along with FMS; if not, write to:                                 !!
!!          Free Software Foundation, Inc.                           !!
!!          59 Temple Place, Suite 330                               !!
!!          Boston, MA  02111-1307  USA                              !!
!! or see:                                                           !!
!!          http://www.gnu.org/licenses/gpl.txt                      !!
!!                                                                   !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

*/

#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include "mpp.h"
#include "mpp_io.h"
#include "constant.h"
#include "tool_util.h"
#include "read_mosaic.h"
#include "create_xgrid.h"
#define AREA_RATIO_THRESH (1.e-6)

char *usage[] = {
  "",
  "  make_quick_mosaic --input_mosaic input_mosaic.nc [--mosaic_name mosaic_name]",
  "                                                                              ",
  "make_quick_mosaic generate a complete grid a FMS coupler. It takes a coupled  ",
  "mosac as input. The atmosphere and ocean grid in output mosaic will be chosen ",
  "to be the same as the land grid of input mosaic. The land/sea mask will be    ",
  "decided by the land/sea mask of input land grid.                              ",
  "                                                                              ",
  "make_coupler_mosaic takes the following flags:                                ",
  "                                                                              ",
  "REQUIRED:                                                                     ",
  "                                                                              ",
  "--input_mosaic input_mosaic.nc specify the input mosaic file name.            ",
  "                                                                              ",
  "OPTIONAL FLAGS                                                                ",
  "                                                                              ",
  "--mosaic_name mosaic_name coupler mosaic name. The output coupler mosaic file ",
  "will be mosaic_name.nc. default value is 'quick_mosaic'                       ",
  "",
  NULL };


char grid_version[] = "0.2";
char tagname[] = "$Name: mom4p1_pubrel_dec2009_nnz $";

int main (int argc, char *argv[])
{
  int c, n, i;
  int errflg = (argc == 1);
  int    option_index = 0;
  int fid, vid, nfile_aXl, ntiles;

  int *nx, *ny;
  double **lonb, **latb, **land_area, **ocean_area, **cell_area;
  char **axl_file = NULL, **axo_file=NULL, **lxo_file=NULL;
  char **ocn_topog_file = NULL;
  char *input_mosaic = NULL;
  char mosaic_name[STRING] = "mosaic", mosaic_file[STRING];
  char griddir[STRING], lnd_mosaic[STRING], filepath[STRING];
  char history[512];
  
  static struct option long_options[] = {
    {"input_mosaic",       required_argument, NULL, 'i'},
    {"mosaic_name",        required_argument, NULL, 'm'},
    {NULL, 0, NULL, 0}
  };

  
  /*
   * process command line
   */

  while ((c = getopt_long(argc, argv, "i:", long_options, &option_index) ) != -1)
    switch (c) {
    case 'i': 
      input_mosaic = optarg;
      break;
    case 'm':
      strcpy(mosaic_name,optarg);
      break;
        case '?':
      errflg++;
    }
  if (errflg || !input_mosaic)  {
    char **u = usage;
    while (*u) { fprintf(stderr, "%s\n", *u); u++; }
    exit(2);
  }  

  strcpy(history,argv[0]);

  for(i=1;i<argc;i++) {
    strcat(history, " ");
    strcat(history, argv[i]);
  }
  
  /* First get land grid information */
  
  mpp_init(&argc, &argv);
  sprintf(mosaic_file, "%s.nc", mosaic_name);
  get_file_path(input_mosaic, griddir);  
  fid = mpp_open(input_mosaic, MPP_READ);
  vid = mpp_get_varid(fid, "lnd_mosaic_file");
  mpp_get_var_value(fid, vid, lnd_mosaic);
  nfile_aXl = mpp_get_dimlen( fid, "nfile_aXl");

  /*make sure the directory that stores the mosaic_file is not current directory */
  { 
    char cur_path[STRING];
    
    if(getcwd(cur_path, STRING) != cur_path ) mpp_error("make_quick_mosaic: The size of cur_path maybe is not big enough");
    printf("The current directory is %s\n", cur_path);
    printf("The mosaic file location is %s\n", griddir);
    if(strcmp(griddir, cur_path)==0 || strcmp( griddir, ".")==0)
      mpp_error("make_quick_mosaic: The input mosaic file location should not be current directory");
  }
  
  sprintf(filepath, "%s/%s", griddir, lnd_mosaic);
  ntiles = read_mosaic_ntiles(filepath);
  /* copy the lnd_mosaic file and grid file */
  {
    int fid2, vid2;
    char cmd[STRING], gridfile[STRING];
    size_t start[4], nread[4];
    sprintf(cmd, "cp %s %s", filepath, lnd_mosaic);

    system(cmd);
    fid2 = mpp_open(filepath, MPP_READ);
    vid2 = mpp_get_varid(fid2, "gridfiles");
    for(i=0; i<4; i++) {
      start[i] = 0; nread[i] = 1;
    }	  
    for(n=0; n<ntiles; n++) {  
      start[0] = n; nread[1] = STRING;
      mpp_get_var_value_block(fid2, vid2, start, nread, gridfile);
      sprintf(cmd, "cp %s/%s %s", griddir, gridfile, gridfile);
      printf("%s \n", cmd);
      system(cmd);
    }
    mpp_close(fid2);
  }

  /* ntiles should be either 1 or ntiles = nfile_aXl */
  if(ntiles != nfile_aXl && ntiles != 1)
    mpp_error("make_quick_mosaic: only support ntiles = 1 or ntiles = nfile_aXl, contact developer");

  nx = (int *)malloc(ntiles*sizeof(int));
  ny = (int *)malloc(ntiles*sizeof(int));
  read_mosaic_grid_sizes(filepath, nx, ny);
  lonb = (double **)malloc(ntiles*sizeof(double *));
  latb = (double **)malloc(ntiles*sizeof(double *));
  for(n=0; n<ntiles; n++) {
     lonb[n] = (double *)malloc((nx[n]+1)*(ny[n]+1)*sizeof(double));
     latb[n] = (double *)malloc((nx[n]+1)*(ny[n]+1)*sizeof(double));     
     read_mosaic_grid_data(filepath, "x", nx[n], ny[n], lonb[n], n, 0, 0); 
     read_mosaic_grid_data(filepath, "y", nx[n], ny[n], latb[n], n, 0, 0);
     for(i=0; i<(nx[n]+1)*(ny[n]+1); i++) {
       lonb[n][i] *= (M_PI/180.0);
       latb[n][i] *= (M_PI/180.0);
     }
  }

  
  /* read the exchange grid information and get the land/sea mask of land model*/
  land_area = (double **)malloc(ntiles*sizeof(double *));
  for(n=0; n<ntiles; n++) {
    land_area[n] = (double *)malloc(nx[n]*ny[n]*sizeof(double));
    for(i=0; i<nx[n]*ny[n]; i++) land_area[n][i] = 0;
  }

  vid = mpp_get_varid(fid, "aXl_file");
  for(n=0; n<nfile_aXl; n++) {
    size_t start[4], nread[4];
    int nxgrid;
    char aXl_file[STRING];
    start[0] = n;
    start[1] = 0;
    nread[0] = 1;
    nread[1] = STRING;
    mpp_get_var_value_block(fid, vid, start, nread, aXl_file);
    sprintf(filepath, "%s/%s", griddir, aXl_file);
    nxgrid = read_mosaic_xgrid_size(filepath);
    if(nxgrid>0) {
      int l;
      int *i1, *j1, *i2, *j2;
      double *area;

      i1 = (int *)malloc(nxgrid*sizeof(int));
      j1 = (int *)malloc(nxgrid*sizeof(int));
      i2 = (int *)malloc(nxgrid*sizeof(int));
      j2 = (int *)malloc(nxgrid*sizeof(int));
      area = (double *)malloc(nxgrid*sizeof(double));
      read_mosaic_xgrid_order1(filepath, i1, j1, i2, j2, area);
      if(ntiles == 1) {
	for(l=0; l<nxgrid; l++) land_area[0][j2[l]*nx[0]+i2[l]] += (area[l]*4*M_PI*RADIUS*RADIUS);
      }
      else {
	for(l=0; l<nxgrid; l++) land_area[n][j2[l]*nx[n]+i2[l]] += (area[l]*4*M_PI*RADIUS*RADIUS);
      }
      free(i1);
      free(j1);
      free(i2);
      free(j2);
      free(area);
    }
  }

  mpp_close(fid);
  
  /* calculate ocean area */
  ocean_area = (double **)malloc(ntiles*sizeof(double *));
  cell_area = (double **)malloc(ntiles*sizeof(double *));
  for(n=0; n<ntiles; n++) {
    ocean_area[n] = (double *)malloc(nx[n]*ny[n]*sizeof(double));
    cell_area[n] = (double *)malloc(nx[n]*ny[n]*sizeof(double));
    get_grid_area(&nx[n], &ny[n], lonb[n], latb[n], cell_area[n]);
    for(i=0; i<nx[n]*ny[n]; i++) {
      ocean_area[n][i] = cell_area[n][i];
      if( fabs(ocean_area[n][i]-land_area[n][i])/ocean_area[n][i] < AREA_RATIO_THRESH )
	ocean_area[n][i] = 0;
      else 
	ocean_area[n][i] -= land_area[n][i];
      if(ocean_area[n][i] < 0) {
	printf("at i = %d, ocean_area = %g, land_area = %g, cell_area=%g\n", i, ocean_area[n][i], land_area[n][i], cell_area[n][i]);
	mpp_error("make_quick_mosaic: ocean area is negative at some points");
      }
    }
  }

  /* write out land mask */
  {
    for(n=0; n<ntiles; n++) {
      int fid, id_mask, dims[2];
      char lnd_mask_file[STRING];
      double *mask;
      mask = (double *)malloc(nx[n]*ny[n]*sizeof(double));
      for(i=0; i<nx[n]*ny[n]; i++) mask[i] = land_area[n][i]/cell_area[n][i];
      
      if(ntiles > 1)
	sprintf(lnd_mask_file, "land_mask_tile%d.nc", n+1);
      else
	strcpy(lnd_mask_file, "land_mask.nc");
      fid = mpp_open(lnd_mask_file, MPP_WRITE);
      mpp_def_global_att(fid, "grid_version", grid_version);
      mpp_def_global_att(fid, "code_version", tagname);
      mpp_def_global_att(fid, "history", history);
      dims[1] = mpp_def_dim(fid, "nx", nx[n]); 
      dims[0] = mpp_def_dim(fid, "ny", ny[n]);
      id_mask = mpp_def_var(fid, "mask", MPP_DOUBLE, 2, dims,  2, "standard_name",
			    "land fraction at T-cell centers", "units", "none");
      mpp_end_def(fid);
      mpp_put_var_value(fid, id_mask, mask);
      free(mask);
      mpp_close(fid);
    }
  }
    
  /* write out ocean mask */
  {
    for(n=0; n<ntiles; n++) {
      int fid, id_mask, dims[2];
      char ocn_mask_file[STRING];
      double *mask;
      mask = (double *)malloc(nx[n]*ny[n]*sizeof(double));
      for(i=0; i<nx[n]*ny[n]; i++) mask[i] = ocean_area[n][i]/cell_area[n][i];
      
      if(ntiles > 1)
	sprintf(ocn_mask_file, "ocean_mask_tile%d.nc", n+1);
      else
	strcpy(ocn_mask_file, "ocean_mask.nc");
      fid = mpp_open(ocn_mask_file, MPP_WRITE);
      mpp_def_global_att(fid, "grid_version", grid_version);
      mpp_def_global_att(fid, "code_version", tagname);
      mpp_def_global_att(fid, "history", history);
      dims[1] = mpp_def_dim(fid, "nx", nx[n]); 
      dims[0] = mpp_def_dim(fid, "ny", ny[n]);
      id_mask = mpp_def_var(fid, "mask", MPP_DOUBLE, 2, dims,  2, "standard_name",
			    "ocean fraction at T-cell centers", "units", "none");
      mpp_end_def(fid);
      mpp_put_var_value(fid, id_mask, mask);
      free(mask);
      mpp_close(fid);
    }
  }
  
  ocn_topog_file = (char **)malloc(ntiles*sizeof(char *));
  axl_file = (char **)malloc(ntiles*sizeof(char *));
  axo_file = (char **)malloc(ntiles*sizeof(char *));
  lxo_file = (char **)malloc(ntiles*sizeof(char *));
  for(n=0; n<ntiles; n++) {
    axl_file[n] = (char *)malloc(STRING*sizeof(char));
    axo_file[n] = (char *)malloc(STRING*sizeof(char));
    lxo_file[n] = (char *)malloc(STRING*sizeof(char));
    ocn_topog_file[n] = (char *)malloc(STRING*sizeof(char));
    sprintf(ocn_topog_file[n], "ocean_topog_tile%d.nc", n+1);
    sprintf(axl_file[n], "atmos_mosaic_tile%dXland_mosaic_tile%d.nc", n+1, n+1);
    sprintf(axo_file[n], "atmos_mosaic_tile%dXocean_mosaic_tile%d.nc", n+1, n+1);
    sprintf(lxo_file[n], "land_mosaic_tile%dXocean_mosaic_tile%d.nc", n+1, n+1);
  }
  
  
  for(n=0; n<ntiles; n++) {
    int *i1, *j1, *i2, *j2;
    double *area, *di, *dj;
    int nxgrid, i, j;
    int fid, dim_string, dim_ncells, dim_two, dims[2];
    int id_contact, id_tile1_cell, id_tile2_cell;
    int id_xgrid_area, id_tile1_dist, id_tile2_dist;
    size_t start[4], nwrite[4];
    char contact[STRING];

    for(i=0; i<4; i++) {
      start[i] = 0; nwrite[i] = 1;
    }	  
    
    /* first calculate the atmXlnd exchange grid */
    i1 = (int *)malloc(nx[n]*ny[n]*sizeof(int));
    j1 = (int *)malloc(nx[n]*ny[n]*sizeof(int));
    i2 = (int *)malloc(nx[n]*ny[n]*sizeof(int));
    j2 = (int *)malloc(nx[n]*ny[n]*sizeof(int));
    area = (double *)malloc(nx[n]*ny[n]*sizeof(double));
    di   = (double *)malloc(nx[n]*ny[n]*sizeof(double));
    dj   = (double *)malloc(nx[n]*ny[n]*sizeof(double));

    /* write out the atmosXland exchange grid file, The file name will be atmos_mosaic_tile#Xland_mosaic_tile#.nc */   
    nxgrid = 0;
    for(j=0; j<ny[n]; j++) for(i=0; i<nx[n]; i++) {
      if(land_area[n][j*nx[n]+i] >0) {
	i1[nxgrid] = i+1;
	j1[nxgrid] = j+1;
	i2[nxgrid] = i+1;
	j2[nxgrid] = j+1;
	area[nxgrid] = land_area[n][j*nx[n]+i];
	di[nxgrid] = 0;
	dj[nxgrid] = 0;
	nxgrid++;
      }
    }
 
    fid = mpp_open(axl_file[n], MPP_WRITE);
    sprintf(contact, "atmos_mosaic:tile%d::land_mosaic:tile%d", n+1, n+1);
    mpp_def_global_att(fid, "grid_version", grid_version);
    mpp_def_global_att(fid, "code_version", tagname);
    mpp_def_global_att(fid, "history", history);
    dim_string = mpp_def_dim(fid, "string", STRING);
    dim_ncells = mpp_def_dim(fid, "ncells", nxgrid);
    dim_two    = mpp_def_dim(fid, "two", 2);
    id_contact = mpp_def_var(fid, "contact", MPP_CHAR, 1, &dim_string, 7, "standard_name", "grid_contact_spec",
			     "contact_type", "exchange", "parent1_cell",
			     "tile1_cell", "parent2_cell", "tile2_cell", "xgrid_area_field", "xgrid_area", 
			     "distant_to_parent1_centroid", "tile1_distance", "distant_to_parent2_centroid", "tile2_distance");
	    
    dims[0] = dim_ncells; dims[1] = dim_two;
    id_tile1_cell = mpp_def_var(fid, "tile1_cell", MPP_INT, 2, dims, 1, "standard_name", "parent_cell_indices_in_mosaic1");
    id_tile2_cell = mpp_def_var(fid, "tile2_cell", MPP_INT, 2, dims, 1, "standard_name", "parent_cell_indices_in_mosaic2");
    id_xgrid_area = mpp_def_var(fid, "xgrid_area", MPP_DOUBLE, 1, &dim_ncells, 2, "standard_name",
				"exchange_grid_area", "units", "m2");
    id_tile1_dist = mpp_def_var(fid, "tile1_distance", MPP_DOUBLE, 2, dims, 1, "standard_name", "distance_from_parent1_cell_centroid");
    id_tile2_dist = mpp_def_var(fid, "tile2_distance", MPP_DOUBLE, 2, dims, 1, "standard_name", "distance_from_parent2_cell_centroid");
    mpp_end_def(fid);
    nwrite[0] = strlen(contact);
    mpp_put_var_value_block(fid, id_contact, start, nwrite, contact);
    nwrite[0] = nxgrid;
    mpp_put_var_value(fid, id_xgrid_area, area);
    mpp_put_var_value_block(fid, id_tile1_cell, start, nwrite, i1);
    mpp_put_var_value_block(fid, id_tile2_cell, start, nwrite, i2);
    mpp_put_var_value_block(fid, id_tile1_dist, start, nwrite, di);
    mpp_put_var_value_block(fid, id_tile2_dist, start, nwrite, di);    
    start[1] = 1;
    mpp_put_var_value_block(fid, id_tile1_cell, start, nwrite, j1);
    mpp_put_var_value_block(fid, id_tile2_cell, start, nwrite, j2);
    mpp_put_var_value_block(fid, id_tile1_dist, start, nwrite, dj);
    mpp_put_var_value_block(fid, id_tile2_dist, start, nwrite, dj);   
    mpp_close(fid);
    
    /* write out the atmosXocean exchange grid file, The file name will be atmos_mosaic_tile#Xocean_mosaic_tile#.nc */
    nxgrid = 0;
    for(j=0; j<ny[n]; j++) for(i=0; i<nx[n]; i++) {
      if(ocean_area[n][j*nx[n]+i] >0) {
	i1[nxgrid] = i+1;
	j1[nxgrid] = j+1;
	i2[nxgrid] = i+1;
	j2[nxgrid] = j+1;
	area[nxgrid] = ocean_area[n][j*nx[n]+i];
	di[nxgrid] = 0;
	dj[nxgrid] = 0;
	nxgrid++;
      }
    }    
    fid = mpp_open(axo_file[n], MPP_WRITE);
    sprintf(contact, "atmos_mosaic:tile%d::ocean_mosaic:tile%d", n+1, n+1);
    mpp_def_global_att(fid, "grid_version", grid_version);
    mpp_def_global_att(fid, "code_version", tagname);
    mpp_def_global_att(fid, "history", history);
    dim_string = mpp_def_dim(fid, "string", STRING);
    dim_ncells = mpp_def_dim(fid, "ncells", nxgrid);
    dim_two    = mpp_def_dim(fid, "two", 2);
    id_contact = mpp_def_var(fid, "contact", MPP_CHAR, 1, &dim_string, 7, "standard_name", "grid_contact_spec",
			     "contact_type", "exchange", "parent1_cell",
			     "tile1_cell", "parent2_cell", "tile2_cell", "xgrid_area_field", "xgrid_area", 
			     "distant_to_parent1_centroid", "tile1_distance", "distant_to_parent2_centroid", "tile2_distance");
	    
    dims[0] = dim_ncells; dims[1] = dim_two;
    id_tile1_cell = mpp_def_var(fid, "tile1_cell", MPP_INT, 2, dims, 1, "standard_name", "parent_cell_indices_in_mosaic1");
    id_tile2_cell = mpp_def_var(fid, "tile2_cell", MPP_INT, 2, dims, 1, "standard_name", "parent_cell_indices_in_mosaic2");
    id_xgrid_area = mpp_def_var(fid, "xgrid_area", MPP_DOUBLE, 1, &dim_ncells, 2, "standard_name",
				"exchange_grid_area", "units", "m2");
    id_tile1_dist = mpp_def_var(fid, "tile1_distance", MPP_DOUBLE, 2, dims, 1, "standard_name", "distance_from_parent1_cell_centroid");
    id_tile2_dist = mpp_def_var(fid, "tile2_distance", MPP_DOUBLE, 2, dims, 1, "standard_name", "distance_from_parent2_cell_centroid");
    mpp_end_def(fid);
    start[1] = 0;
    nwrite[0] = strlen(contact);
    mpp_put_var_value_block(fid, id_contact, start, nwrite, contact);
    nwrite[0] = nxgrid;
    mpp_put_var_value(fid, id_xgrid_area, area);
    mpp_put_var_value_block(fid, id_tile1_cell, start, nwrite, i1);
    mpp_put_var_value_block(fid, id_tile2_cell, start, nwrite, i2);
    mpp_put_var_value_block(fid, id_tile1_dist, start, nwrite, di);
    mpp_put_var_value_block(fid, id_tile2_dist, start, nwrite, di);    
    start[1] = 1;
    mpp_put_var_value_block(fid, id_tile1_cell, start, nwrite, j1);
    mpp_put_var_value_block(fid, id_tile2_cell, start, nwrite, j2);
    mpp_put_var_value_block(fid, id_tile1_dist, start, nwrite, dj);
    mpp_put_var_value_block(fid, id_tile2_dist, start, nwrite, dj);   
    mpp_close(fid);
    
    /* write out landXocean exchange grid information */
    fid = mpp_open(lxo_file[n], MPP_WRITE);
    sprintf(contact, "land_mosaic:tile%d::ocean_mosaic:tile%d", n+1, n+1);
    mpp_def_global_att(fid, "grid_version", grid_version);
    mpp_def_global_att(fid, "code_version", tagname);
    mpp_def_global_att(fid, "history", history);
    dim_string = mpp_def_dim(fid, "string", STRING);
    dim_ncells = mpp_def_dim(fid, "ncells", nxgrid);
    dim_two    = mpp_def_dim(fid, "two", 2);
    id_contact = mpp_def_var(fid, "contact", MPP_CHAR, 1, &dim_string, 7, "standard_name", "grid_contact_spec",
			     "contact_type", "exchange", "parent1_cell",
			     "tile1_cell", "parent2_cell", "tile2_cell", "xgrid_area_field", "xgrid_area", 
			     "distant_to_parent1_centroid", "tile1_distance", "distant_to_parent2_centroid", "tile2_distance");
	    
    dims[0] = dim_ncells; dims[1] = dim_two;
    id_tile1_cell = mpp_def_var(fid, "tile1_cell", MPP_INT, 2, dims, 1, "standard_name", "parent_cell_indices_in_mosaic1");
    id_tile2_cell = mpp_def_var(fid, "tile2_cell", MPP_INT, 2, dims, 1, "standard_name", "parent_cell_indices_in_mosaic2");
    id_xgrid_area = mpp_def_var(fid, "xgrid_area", MPP_DOUBLE, 1, &dim_ncells, 2, "standard_name",
				"exchange_grid_area", "units", "m2");
    id_tile1_dist = mpp_def_var(fid, "tile1_distance", MPP_DOUBLE, 2, dims, 1, "standard_name", "distance_from_parent1_cell_centroid");
    id_tile2_dist = mpp_def_var(fid, "tile2_distance", MPP_DOUBLE, 2, dims, 1, "standard_name", "distance_from_parent2_cell_centroid");
    mpp_end_def(fid);
    start[1] = 0;
    nwrite[0] = strlen(contact);
    mpp_put_var_value_block(fid, id_contact, start, nwrite, contact);
    nwrite[0] = nxgrid;
    mpp_put_var_value(fid, id_xgrid_area, area);
    mpp_put_var_value_block(fid, id_tile1_cell, start, nwrite, i1);
    mpp_put_var_value_block(fid, id_tile2_cell, start, nwrite, i2);
    mpp_put_var_value_block(fid, id_tile1_dist, start, nwrite, di);
    mpp_put_var_value_block(fid, id_tile2_dist, start, nwrite, di);    
    start[1] = 1;
    mpp_put_var_value_block(fid, id_tile1_cell, start, nwrite, j1);
    mpp_put_var_value_block(fid, id_tile2_cell, start, nwrite, j2);
    mpp_put_var_value_block(fid, id_tile1_dist, start, nwrite, dj);
    mpp_put_var_value_block(fid, id_tile2_dist, start, nwrite, dj);   
    mpp_close(fid);
    
    free(i1);
    free(j1);
    free(i2);
    free(j2);
    free(area);
    free(di);
    free(dj);
  }
  
  /*Fianlly create the coupler mosaic file mosaic_name.nc */
  {
    int dim_string, dim_axo, dim_axl, dim_lxo, dims[2];
    int id_amosaic_file, id_lmosaic_file, id_omosaic_file, id_otopog_file;
    int id_axo_file, id_axl_file, id_lxo_file;
    size_t start[4], nwrite[4];

    for(i=0; i<4; i++) {
      start[i] = 0; nwrite[i] = 1;
    }	      
    printf("mosaic_file is %s\n", mosaic_file);   
    fid = mpp_open(mosaic_file, MPP_WRITE);
    mpp_def_global_att(fid, "grid_version", grid_version);
    mpp_def_global_att(fid, "code_version", tagname);
    mpp_def_global_att(fid, "history", history);
    dim_string = mpp_def_dim(fid, "string", STRING);
    dim_axo = mpp_def_dim(fid, "nfile_aXo", ntiles);
    dim_axl = mpp_def_dim(fid, "nfile_aXl", ntiles);
    dim_lxo = mpp_def_dim(fid, "nfile_lXo", ntiles);    
    id_amosaic_file = mpp_def_var(fid, "atm_mosaic_file", MPP_CHAR, 1, &dim_string,
				  1, "standard_name", "atmosphere_mosaic_file_name");
    id_lmosaic_file = mpp_def_var(fid, "lnd_mosaic_file", MPP_CHAR, 1, &dim_string,
				  1, "standard_name", "land_mosaic_file_name");
    id_omosaic_file = mpp_def_var(fid, "ocn_mosaic_file", MPP_CHAR, 1, &dim_string,
				  1, "standard_name", "ocean_mosaic_file_name");
    id_otopog_file  = mpp_def_var(fid, "ocn_topog_file", MPP_CHAR, 1, &dim_string,
				  1, "standard_name", "ocean_topog_file_name");    
    dims[0] = dim_axo; dims[1] = dim_string;
    id_axo_file = mpp_def_var(fid, "aXo_file", MPP_CHAR, 2, dims, 1, "standard_name", "atmXocn_exchange_grid_file");
    dims[0] = dim_axl; dims[1] = dim_string;
    id_axl_file = mpp_def_var(fid, "aXl_file", MPP_CHAR, 2, dims, 1, "standard_name", "atmXlnd_exchange_grid_file");
    dims[0] = dim_lxo; dims[1] = dim_string;
    id_lxo_file = mpp_def_var(fid, "lXo_file", MPP_CHAR, 2, dims, 1, "standard_name", "lndXocn_exchange_grid_file");
    mpp_end_def(fid);    
    nwrite[0] = strlen(lnd_mosaic);
    mpp_put_var_value_block(fid, id_lmosaic_file, start, nwrite, lnd_mosaic);
    mpp_put_var_value_block(fid, id_amosaic_file, start, nwrite, lnd_mosaic);
    mpp_put_var_value_block(fid, id_omosaic_file, start, nwrite, lnd_mosaic);
    for(n=0; n<ntiles; n++) {
      start[0] = n; nwrite[0] =1;
      nwrite[1] = strlen(ocn_topog_file[n]);
      mpp_put_var_value_block(fid, id_otopog_file, start, nwrite, ocn_topog_file[n]);
      nwrite[1] = strlen(axl_file[n]);
      mpp_put_var_value_block(fid, id_axl_file, start, nwrite, axl_file[n]);
      nwrite[1] = strlen(axo_file[n]);
      mpp_put_var_value_block(fid, id_axo_file, start, nwrite, axo_file[n]);
      nwrite[1] = strlen(lxo_file[n]);
      mpp_put_var_value_block(fid, id_lxo_file, start, nwrite, lxo_file[n]);
    }
    mpp_close(fid);
  }
    
  for(n=0; n<ntiles; n++) {
    free(axl_file[n]);
    free(axo_file[n]);
    free(lxo_file[n]);
  }
  free(axl_file);
  free(axo_file);
  free(lxo_file);
  printf("\n***** Congratulation! You have successfully run make_quick_mosaic\n");
  mpp_end();

  return 0;

} /* main */  
