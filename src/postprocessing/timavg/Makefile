# The following three directory may need to set.
NETCDFPATH = /usr/local/netcdf-4.0_hdf
HDF5PATH = /usr/local/hdf5-1.8.1
POSTPPATH  = .
FFLAGS = -I${NETCDFPATH}/include -I${HDF5PATH}/include -I/usr/include -fltconsistency -fno-alias -stack_temps -safe_cray_ptr -ftz -i-static -static -assume byterecl -g -O2 -i4 -r8 -nowarn -Wp,-w
LDFLAGS =  -L${NETCDFPATH}/lib -L${HDF5PATH}/lib -lmpi  -lnetcdff -lnetcdf -lhdf5_hl -lhdf5 -lz 
DEFFLAG = -Duse_netCDF
LNFLAGS = -v
FC      = ifort

OBJS = time_average.o

TAVG.v10.exe: $(OBJS)
	$(FC) $(LNFLAGS) -o $@ $(OBJS) $(LDFLAGS)

time_average.o: $(POSTPPATH)/time_average.f90
	$(FC) $(FFLAGS) -c $(POSTPPATH)/time_average.f90
