module interp_module

  implicit none

  contains
!-----------------------------------------------------------------------
  function interp3d(z,x,y,zp0,zp1,xp0,xp1,yp0,yp1,fp,zlog) result(f)
! equidistant 3d interpolation, constant extrapolation
! 3d dimensions in the order of z,x,y in accordance with model fields
! zp0,zp1,xp0,xp1,yp0,yp1: input domain; fp: input values
! z,x,y: output locations; f: output values

    real,intent(in) :: zp0,zp1,xp0,xp1,yp0,yp1
    real,dimension(:),intent(in) :: z,x,y
    real,dimension(:,:,:),intent(in) :: fp
    logical,intent(in),optional :: zlog
    real,dimension(size(z),size(x),size(y)) :: f

    integer :: nz,nx,ny,nzp,nxp,nyp,k,i,j
    real :: dz,dx,dy
    integer,dimension(size(z)) :: k0,k1
    integer,dimension(size(x)) :: i0,i1
    integer,dimension(size(y)) :: j0,j1
    real,dimension(size(z)) :: z0,z1
    real,dimension(size(x)) :: x0,x1
    real,dimension(size(y)) :: y0,y1
    real,dimension(size(z),size(x),size(y)) :: za,zb,xa,xb,ya,yb
    real,dimension(size(z)) :: zap,zbp
    real,dimension(size(x),size(y)) :: xap,xbp,yap,ybp
    real,dimension(size(fp,1),size(x),size(y)) :: fxy

    nz = size(z)
    nx = size(x)
    ny = size(y)

    nzp = size(fp,1)
    nxp = size(fp,2)
    nyp = size(fp,3)

    dz = (zp1-zp0)/(nzp-1)
    dx = (xp1-xp0)/(nxp-1)
    dy = (yp1-yp0)/(nyp-1)

    k0 = max(int((z-zp0)/dz),0)+1
    k1 = min(k0+1,nzp)

    i0 = max(int((x-xp0)/dx),0)+1
    i1 = min(i0+1,nxp)

    j0 = max(int((y-yp0)/dy),0)+1
    j1 = min(j0+1,nyp)

    z0 = zp0+(k0-1)*dz
    z1 = z0+dz

    x0 = xp0+(i0-1)*dx
    x1 = x0+dx

    y0 = yp0+(j0-1)*dy
    y1 = y0+dy

! logarithmic interpolation in z direction if zlog is set to true
    if (present(zlog)) then
      if (zlog) then
        forall (j=1:ny) xap(:,j) = x1-x
        forall (j=1:ny) xbp(:,j) = x-x0
        forall (i=1:nx) yap(i,:) = y1-y
        forall (i=1:nx) ybp(i,:) = y-y0

        forall (k=1:nzp) fxy(k,:,:) = &
          (xap*yap*fp(k,i0,j0) + xap*ybp*fp(k,i0,j1) + &
          xbp*yap*fp(k,i1,j0) + xbp*ybp*fp(k,i1,j1)) / (dx*dy)

        zap = z1-z
        zbp = z-z0

        forall (i=1:nx,j=1:ny) f(:,i,j) = &
          exp((zap*log(fxy(k0,i,j)) + zbp*log(fxy(k1,i,j))) / dz)

        return
      endif
    endif

! otherwise (zlog not present or set to false), linear interpolation
    forall (i=1:nx,j=1:ny) za(:,i,j) = z1-z
    forall (i=1:nx,j=1:ny) zb(:,i,j) = z-z0
    forall (k=1:nz,j=1:ny) xa(k,:,j) = x1-x
    forall (k=1:nz,j=1:ny) xb(k,:,j) = x-x0
    forall (k=1:nz,i=1:nx) ya(k,i,:) = y1-y
    forall (k=1:nz,i=1:nx) yb(k,i,:) = y-y0

    f = (za*xa*ya*fp(k0,i0,j0) + za*xa*yb*fp(k0,i0,j1) + &
      za*xb*ya*fp(k0,i1,j0) + za*xb*yb*fp(k0,i1,j1) + &
      zb*xa*ya*fp(k1,i0,j0) + zb*xa*yb*fp(k1,i0,j1) + &
      zb*xb*ya*fp(k1,i1,j0) + zb*xb*yb*fp(k1,i1,j1)) / (dz*dx*dy)

  end function interp3d
!-----------------------------------------------------------------------
  function interp2d(x,y,xp0,xp1,yp0,yp1,fp) result(f)
! equidistant 2d interpolation, constant extrapolation
! xp0,xp1,yp0,yp1: input domain; fp: input values
! x,y: output locations; f: output values

    real,intent(in) :: xp0,xp1,yp0,yp1
    real,dimension(:),intent(in) :: x,y
    real,dimension(:,:),intent(in) :: fp
    real,dimension(size(x),size(y)) :: f

    integer :: nx,ny,nxp,nyp,i,j
    real :: dx,dy
    integer,dimension(size(x)) :: i0,i1
    integer,dimension(size(y)) :: j0,j1
    real,dimension(size(x)) :: x0,x1
    real,dimension(size(y)) :: y0,y1
    real,dimension(size(x),size(y)) :: xa,xb,ya,yb

    nx = size(x)
    ny = size(y)

    nxp = size(fp,1)
    nyp = size(fp,2)

    dx = (xp1-xp0)/(nxp-1)
    dy = (yp1-yp0)/(nyp-1)

    i0 = max(int((x-xp0)/dx),0)+1
    i1 = min(i0+1,nxp)

    j0 = max(int((y-yp0)/dy),0)+1
    j1 = min(j0+1,nyp)

    x0 = xp0+(i0-1)*dx
    x1 = x0+dx

    y0 = yp0+(j0-1)*dy
    y1 = y0+dy

    forall (j=1:ny) xa(:,j) = x1-x
    forall (j=1:ny) xb(:,j) = x-x0
    forall (i=1:nx) ya(i,:) = y1-y
    forall (i=1:nx) yb(i,:) = y-y0

    f = (xa*ya*fp(i0,j0) + xa*yb*fp(i0,j1) + &
      xb*ya*fp(i1,j0) + xb*yb*fp(i1,j1)) / (dx*dy)

  end function interp2d
!-----------------------------------------------------------------------
  function interp1d(x,xp0,xp1,fp) result(f)
! equidistant 1d interpolation, constant extrapolation
! xp0,xp1: input domain; fp: input values
! x: output locations; f: output values

    real,intent(in) :: xp0,xp1
    real,dimension(:),intent(in) :: x,fp
    real,dimension(size(x)) :: f

    integer :: nxp
    real :: dx
    integer,dimension(size(x)) :: i0,i1
    real,dimension(size(x)) :: x0,x1,xa,xb

    nxp = size(fp)

    dx = (xp1-xp0)/(nxp-1)

    i0 = max(int((x-xp0)/dx),0)+1
    i1 = min(i0+1,nxp)

    x0 = xp0+(i0-1)*dx
    x1 = x0+dx

    xa = x1-x
    xb = x-x0

    f = (xa*fp(i0) + xb*fp(i1)) / dx

  end function interp1d
!-----------------------------------------------------------------------
end module interp_module
