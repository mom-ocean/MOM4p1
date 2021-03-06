#!/bin/csh -f
# Minimal runscript: CM2.1p1


set platform      = ia64           # A unique identifier for your platform
set name          = ICCMp1        # Name of the experiment you want to run
                                   # Currently one of: CM2.1p1
set npes          = 54             # number of processor
                                   # Note: If you change npes you may need to change
                                   # the layout in the corresponding namelist
set type          = ICCMp1         # type of the experiment
set root          = $cwd:h         # The directory in which you checked out src
set code_dir      = $root/src      # source code directory
set workdir       = $root/work     # where the model is run and model output is produced
                                   # This is recommended to be a link to the $WORKDIR of the platform.
set expdir        = $workdir/$name
set inputDataDir  = $expdir/INPUT   # This is path to the directory that contains the input data for this experiment.
                                     # You should have downloaded and untared this directory from MOM4p1 FTP site.
set diagtable     = $inputDataDir/diag_table  # path to diagnositics table
set datatable     = $inputDataDir/data_table  # path to the data override table.
set fieldtable    = $inputDataDir/field_table # path to the field table
set namelist      = $inputDataDir/input.nml   # path to namelist file

set executable    = $root/exec_$platform/$type/fms_$type.x      # executable created after compilation
set static        = 0              # 1 if you want static memory allocation, 0 for dynamic
if($static)    set executable = $root/exec_$platform/${type}_static/fms_$type.x 

set archive       = /archive/nnz/fms/mom4p1_pubrel_dec2009/$name #Large directory to host the input and output data.

#===========================================================================
# The user need not change any of the following
#===========================================================================

#
# Users must ensure the correct environment file exists for their platform.
#
source $root/bin/environs.$platform  # environment variables and loadable modules

set mppnccombine  = $root/bin/mppnccombine.$platform  # path to executable mppnccombine
set time_stamp    = $root/bin/time_stamp.csh          # path to cshell to generate the date

# Check if the user has extracted the input data
  if ( ! -d $inputDataDir ) then

#    if( -e $archive/$name.input.tar.gz ) then
#      cd $workdir
#      tar zxvf $archive/$name.input.tar.gz
#    else  

    echo "ERROR: the experiment directory '$inputDataDir' does not exist or does not contain input and preprocessing data directories!"
    echo "Please make sure that the variable 'name' in this script is set to one of CM2.1p1 ."
    echo "Then download and extract here the tar ball corresponding to this experiment from GFDL anonymous ftp site!"
   echo " cd  $workdir"
    echo " wget ftp.gfdl.noaa.gov:/perm/MOM4/mom4p1_pubrel_dec2009/exp/$name.input.tar.gz"
    echo " tar zxvf $name.input.tar.gz" 
    echo "Then rerun this script."
    exit 1

#    endif
  endif

set echo

# setup directory structure
  if ( ! -d $expdir )         mkdir -p $expdir
  if ( ! -d $expdir/RESTART ) mkdir -p $expdir/RESTART

# --- make sure executable is up to date ---
  set makeFile      = Make_$type
  cd $executable:h
  make -f $makeFile
  if ( $status != 0 ) then
    unset echo
    echo "ERROR: make failed"
    exit 1
  endif
#-------------------------------------------

#Change to expdir

  cd $expdir

# Create INPUT directory. Make a link instead of copy
# 
if ( ! -d $expdir/INPUT   ) mkdir -p $expdir/INPUT

  if ( ! -e $namelist ) then
    echo "ERROR: required input file does not exist $namelist "
    exit 1
  endif
  if ( ! -e $datatable ) then
    echo "ERROR: required input file does not exist $datatable "
    exit 1
  endif
  if ( ! -e $diagtable ) then
    echo "ERROR: required input file does not exist $diagtable "
    exit 1
  endif
  if ( ! -e $fieldtable ) then
    echo "ERROR: required input file does not exist $fieldtable "
    exit 1
  endif

  cp $namelist   input.nml
  cp $datatable  data_table
  cp $diagtable  diag_table
  cp $fieldtable field_table 


#   --- run the model ---

  cd $expdir

  if($npes > 1) then
     mpirun -np $npes $executable > fms.out
  else
     $executable:t > fms.out
  endif

#----------------------------------------------------------------------------------------------
# generate date for file names ---
    set begindate = `$time_stamp -bf digital`
    if ( $begindate == "" ) set begindate = tmp`date '+%j%H%M%S'`
    set enddate = `$time_stamp -ef digital`
    if ( $enddate == "" ) set enddate = tmp`date '+%j%H%M%S'`
    if ( -f time_stamp.out ) rm -f time_stamp.out
#----------------------------------------------------------------------------------------------
# get a tar restart file
  cd RESTART
  cp $expdir/input.nml .
  cp $expdir/*_table .
# combine netcdf files
  if ( $npes > 1 ) then
    set file_previous = ""
    set multires = (`ls *.nc.????`)
    foreach file ( $multires )
	if ( $file:r != $file_previous:r ) then
	    set input_files = ( `ls $file:r.????` )
              if ( $#input_files > 0 ) then
                 $mppnccombine $file:r $input_files
                 if ( $status != 0 ) then
                   echo "ERROR: in execution of mppnccombine on restarts"
                   exit 1
                 endif
                 rm $input_files
              endif
           else
              continue
           endif
           set file_previous = $file
       end
  endif

  cd $expdir
  mkdir history
  mkdir ascii
#----------------------------------------------------------------------------------------------
# rename ascii files with the date
  foreach out (`ls *.out`)
     mv $out ascii/$begindate.$out
  end

#----------------------------------------------------------------------------------------------
# combine netcdf files
  if ( $npes > 1 ) then
    set file_previous = ""
    set multires = (`ls *.nc.????`)
    foreach file ( $multires )
	if ( $file:r != $file_previous:r ) then
	    set input_files = ( `ls $file:r.????` )
              if ( $#input_files > 0 ) then
                 $mppnccombine $file:r $input_files
                 if ( $status != 0 ) then
                   echo "ERROR: in execution of mppnccombine on restarts"
                   exit 1
                 endif
                 rm $input_files
              endif
           else
              continue
           endif
           set file_previous = $file
       end
  endif

#----------------------------------------------------------------------------------------------
# rename nc files with the date
  foreach ncfile (`/bin/ls *.nc`)
     mv $ncfile history/$begindate.$ncfile
  end

  unset echo


echo end_of_run
echo "NOTE: Natural end-of-script."

#Archive the results

cd $workdir
tar cvf $name.output.tar --exclude=data_table --exclude=diag_table --exclude=field_table --exclude=fms_$type.x --exclude=input.nml --exclude=INPUT $name
gzip $name.output.tar
mv $name.output.tar.gz $archive/

exit 0
  
