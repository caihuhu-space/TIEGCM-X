      module he_coefs_module

      implicit none

#ifdef SUN
#include <netcdf3.inc>
#else
#include <netcdf.inc>
#endif
      integer :: nlat_he,nlon_he,ntrigs
      integer,dimension(13) :: ifax
      real,dimension(:),allocatable :: trigs
      real,dimension(:),allocatable :: glat_he,glon_he
      real,dimension(:,:,:),allocatable :: pmn,zmn

      contains
!-----------------------------------------------------------------------
      subroutine init_he_coefs

      use input_module,only: he_coefs_ncfile
!      use netcdf_f03,only: nf_open,nf_inq_dimid,nf_inq_dimlen,
!     |  nf_inq_varid,nf_get_var_double,nf_close,NF_NOWRITE

      integer :: stat,ncid,dimid,varid,length,i
      real :: res

      stat = nf_open(he_coefs_ncfile,NF_NOWRITE,ncid)

      stat = nf_inq_dimid(ncid,'lat1',dimid)
      stat = nf_inq_dimlen(ncid,dimid,nlat_he)
      stat = nf_inq_dimid(ncid,'lat2',dimid)
      stat = nf_inq_dimlen(ncid,dimid,length)
      if (length /= nlat_he-1)
     |  call shutdown('Invalid Helium coefficient file')
      stat = nf_inq_dimid(ncid,'lat3',dimid)
      stat = nf_inq_dimlen(ncid,dimid,length)
      if (length /= nlat_he)
     |  call shutdown('Invalid Helium coefficient file')

      allocate(pmn(nlat_he,0:nlat_he-2,0:nlat_he-1))
      allocate(zmn(nlat_he,0:nlat_he-2,0:nlat_he-1))

      stat = nf_inq_varid(ncid,'pmn',varid)
      stat = nf_get_var_double(ncid,varid,pmn)
      stat = nf_inq_varid(ncid,'zmn',varid)
      stat = nf_get_var_double(ncid,varid,zmn)

      stat = nf_close(ncid)

      nlon_he = nlat_he*2
      res = 180./nlat_he

      allocate(glat_he(nlat_he))
      allocate(glon_he(nlon_he))

      glat_he = (/(-90.-res/2.+i*res,i=1,nlat_he)/)
      glon_he = (/(-180.+(i-1)*res,i=1,nlon_he)/)

      ntrigs = 3*nlon_he/2+1
      allocate(trigs(ntrigs))

      end subroutine init_he_coefs
!-----------------------------------------------------------------------
      end module he_coefs_module
