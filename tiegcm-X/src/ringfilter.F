#include <defs.h>
      subroutine ringfilter(f,lev0,lev1,order,recon,lat)

! Coded by Tong Dang, May, 2017; Modified by Haonan Wu, Feb, 2022
! This program is to build a new ring filter for the high-resolution TIE-GCM.
! The filter procedures are as follows:
! 1. Divide the longitude sessions into several chunks
! 2. Calculate the avg values of each chunk
! 3. Do the restruction based on the mean value of each chunk

      use params_module,only: nlon,nlonp2,nlonp4,nlat
      use cons_module,only: dlamda
      implicit none

      integer,intent(in) :: lev0,lev1,order,lat
      real,dimension(nlonp4,lev0:lev1),intent(inout) :: f
      character(len=*),intent(in) :: recon

      integer,dimension(*),parameter :: chunk_array =
#if (NLAT==36 && NLON==72)
     |  (/9,9,18,18,36,36/) ! 6 for 5 deg 6
#elif (NLAT==72 && NLON==144)
     |  (/9,9,18,18,36,36,72,72,72,72/) ! 10 for 2.5 deg
#elif (NLAT==144 && NLON==288)
     |  (/9,9,18,18,36,36,36,36,72,72,
     |  144,144,144,144,144,144,144,144/) ! 18 for 1.25 deg
#elif (NLAT==288 && NLON==576)
     |  (/9,9,9,9,18,18,18,18,36,36,36,36,36,36,36,36,
     |  72,72,72,72,72,72,72,72,144,144,144,144,144,144,144,144,
     |  288,288,288,288,288,288,288,288/) ! 40 for 0.625 deg
#else
      UNKNOWN NLAT,NLON
#endif
      integer,parameter :: nlat_filter = size(chunk_array)

      integer :: i,m,n,nchnk,ichnk,im2,im1,ip1,ip2
      real :: ang
      real,dimension(lev0:lev1) :: a0,left,right,fa,fb
      real,dimension(order,lev0:lev1) :: a,b
      real,dimension(nlon,lev0:lev1) :: w,wm,fx,avg

      if (lat <= nlat_filter) then
        nchnk = chunk_array(lat)
      elseif (lat <= nlat-nlat_filter) then
        return
      else
        nchnk = chunk_array(nlat-lat+1)
      endif

      if (order >= nlon) call shutdown('Filtering order too high')

      w = f(3:nlonp2,:)

! Fourier series to nth order
      a0 = sum(w,1)/nlon
      a = 0
      b = 0
      do i = 1,nlon
        ang = i*dlamda
        do m = 1,order
          a(m,:) = a(m,:)+w(i,:)*cos(m*ang)
          b(m,:) = b(m,:)+w(i,:)*sin(m*ang)
        enddo
      enddo
      a = 2*a/nlon
      b = 2*b/nlon

! reconstruct based on Fourier series (taken as mean)
      do i = 1,nlon
        ang = i*dlamda
        wm(i,:) = a0
        do m = 1,order
          wm(i,:) = wm(i,:)+a(m,:)*cos(m*ang)+b(m,:)*sin(m*ang)
        enddo
      enddo

! deviation field
      fx = w-wm

! divide lon array into chunks and calculate mean in each chunk
      n = nlon/nchnk

      do ichnk = 1,nchnk
        avg(ichnk,:) = sum(fx((ichnk-1)*n+1:ichnk*n,:),1)/n
      enddo

      do ichnk = 1,nchnk
        im2 = ichnk-2
        if (im2 <= 0) im2 = im2+nchnk
        im1 = ichnk-1
        if (im1 <= 0) im1 = im1+nchnk
        ip1 = ichnk+1
        if (ip1 >= nchnk+1) ip1 = ip1-nchnk
        ip2 = ichnk+2
        if (ip2 >= nchnk+1) ip2 = ip2-nchnk
        left = (-avg(im2,:)+7*avg(im1,:)+7*avg(ichnk,:)-avg(ip1,:))/12
        right = (-avg(im1,:)+7*avg(ichnk,:)+7*avg(ip1,:)-avg(ip2,:))/12
        fa = (left+right-2*avg(ichnk,:))/n**2
        fb = (3*avg(ichnk,:)-2*left-right)/n
        select case (trim(recon))
          case ('PLM')
            do m = 1,n
              fx((ichnk-1)*n+m,:) = (left*(n-m)+right*m)/n
            enddo
          case ('PPM')
            do m = 1,n
              fx((ichnk-1)*n+m,:) = fa*(3*m**2-3*m+1)+fb*(2*m-1)+left
            enddo
          case default
            call shutdown('Unknown reconstruction type:'//trim(recon))
        endselect
      enddo

! full field
      w = fx+wm

      f(3:nlonp2,:) = w
      f(1:2,:) = w(nlon-1:nlon,:)
      f(nlonp4-1:nlonp4,:) = w(1:2,:)

      end subroutine ringfilter
