! ###############################################
! Original Code by
! Copyright 2012-2018, University of Strathclyde
! Authors: Lawrence T. Campbell
! License: BSD-3-Clause
! ###############################################
! ###############################################
! Changes by
! Authors: Lars Huebner
! License: BSD-3-Clause
! ###############################################

!!!!!!!!!!!!!!!!!!! Puffin Version 1.9.0 !!!!!!!!!!!!!!!!!!!

!> @author
!> Lars Huebner
!> DESY,
!> Hamburg, DE
!> @brief
!> Interpolation for planepole field profile from user input

module undfieldinterp

use Globals

implicit none

contains

subroutine read_planepolefield(bfile, bfield)


    implicit none

! bfile        -      file to load
! bfield contains
! n     -      number of datapoints in file
! z               -      positions array
! by              -      by field array
! c              -      coefficients for spline


    character(1024_ip), intent(in) :: bfile
    type(bfieldspline), intent(out) :: bfield

    integer(kind=ip) :: ndatapoints

    integer(kind=ip)   :: ios

    open(73,FILE=bfile, IOSTAT=ios, STATUS='OLD', ACTION='READ', POSITION ='REWIND')

    if (ios /= 0) then
        print*, 'iostat = ', ios, "in read_planepolefield"
        stop "OPEN(bfile) not performed correctly, IOSTAT /= 0"
    end if

    ! if (tProcInfo_G%qRoot) print *, 'reading number of datapoints...'
    ! read the number of datapoints
    read(73,*, IOSTAT=ios) ndatapoints ! number
    if (ios < 0) then ! end of file
        if ((tProcInfo_G%qroot) .and. (ioutInfo_G > 2)) print*, "Reached end of bfield file at line 1!"
        stop "bfield file empty!"
    end if
    ! if (tProcInfo_G%qRoot) print *, 'OK, datapoints=', ndatapoints
    bfield%n = ndatapoints

    ! if (tProcInfo_G%qRoot) print *, 'allocate arrays'
    allocate(bfield%z(ndatapoints))
    allocate(bfield%by(ndatapoints))
    !allocate(bfield%bz(ndatapoints))

    ! if (tProcInfo_G%qRoot) print *, 'reading z values...'
    read(73,*,IOSTAT=ios) bfield%z
    if (ios < 0) then ! end of file
        if ((tProcInfo_G%qroot) .and. (ioutInfo_G > 2)) print*, "Reached end of bfield file at line 2!"
        stop "bfield file has only 1 line!"
    end if
    ! if (tProcInfo_G%qRoot) print *, 'OK'

    ! if (tProcInfo_G%qRoot) print *, 'reading by values...'
    read(73,*,IOSTAT=ios) bfield%by
    if (ios < 0) then ! end of file
        if ((tProcInfo_G%qroot) .and. (ioutInfo_G > 2)) print*, "Reached end of bfield file at line 3!"
        stop "bfield file has only 2 lines!"
    end if
    ! if (tProcInfo_G%qRoot) print *, 'OK'

    !read(73,*,IOSTAT=ios) bfield%bz
    !if (ios < 0) then ! end of file
    !    if ((tProcInfo_G%qroot) .and. (ioutInfo_G > 2)) print*, "Reached end of bfield file at line 4!"
    !    stop "bfield file has only 3 lines!"
    !end if

    close(73, STATUS='KEEP')

    allocate(bfield%cy(ndatapoints))
    if (tProcInfo_G%qRoot) print *, 'interpolating field'
    call splineCoeff(ndatapoints,bfield%z,bfield%by,bfield%cy)
    if (tProcInfo_G%qRoot) print *, 'OK'
    !call splineCoeff(ndatapoints,bfield%z,bfield%bz,bfield%cz)
    return
end subroutine read_planepolefield

subroutine splineCoeff(ndatapoints,z,Bi,coeff)
    implicit none
    ! calculate a cubic spline
    ! adapted from
    ! NUMERICAL RECIPES IN FORTRAN 77: THE ART OF SCIENTIFIC COMPUTING (ISBN 0-521-43064-X)
    ! page 109
    
    ! coefficients are second derivative at that point
    integer(kind=ip), intent(in) :: ndatapoints
    real(kind=wp), intent(in) :: z(:)
    real(kind=wp), intent(in) :: Bi(:)
    real(kind=wp), intent(out) :: coeff(:)
    ! helper variables
    real(kind=wp) :: sig, p
    integer(kind=ip) :: i,k
    real(kind=wp), allocatable :: u(:)

    if (ndatapoints > 65536_ip) then! maximum interpolation points. Probably gets too slow if too large
        stop "more than 2**16 (65536) interpolation points. Aborting."
    end if 

    ! helper variables
    allocate(u(ndatapoints))

    ! "natural" boundaries
    coeff(1)=0.0_wp
    coeff(ndatapoints)=0.0_wp
    u(1)=0.0_wp

    ! if (tProcInfo_G%qRoot) print *, 'constructing spline...'
    do i=2,ndatapoints-1
      ! if (tProcInfo_G%qRoot) print *, 'step', i-2_ip, 'of', ndatapoints-2_ip
      sig = (z(i)-z(i-1))/(z(i+1)-z(i-1))
      p = sig*coeff(i-1)+2.0_wp
      coeff(i)=(sig-1.0_wp)/p
      u(i) = (6.0_wp* & 
               ( &
                 (Bi(i+1)-Bi(i)) / (z(i+1)-z(i)) &
                 -(Bi(i)-Bi(i-1)) / (z(i)-z(i-1)) &
               ) &
               / (z(i+1)-z(i-1)) &
               - sig*u(i-1) &
             ) / p
    end do
    do k=ndatapoints-1,1,-1
      coeff(k)=coeff(k)*coeff(k+1)+u(k)
    end do
    deallocate(u)
    return

end subroutine splineCoeff

subroutine evaluateSplineBfield(bspline,zin,klo,khi,fyeval,fydeval)
    implicit none
    ! return value  and first derivative of splines
    ! adapted from
    ! NUMERICAL RECIPES IN FORTRAN 77: THE ART OF SCIENTIFIC COMPUTING (ISBN 0-521-43064-X)
    ! page 110
    type(bfieldspline), intent(in) :: bspline
    real(kind=wp), intent(in) :: zin
    integer(kind=ip), intent(inout) :: klo ! interval index guess low
    integer(kind=ip), intent(inout) :: khi ! interval index guess high
    real(kind=wp), intent(out) :: fyeval,fydeval !,fzeval ,fzdeval ! return values
    ! helper variables
    integer(kind=ip) :: k ! position
    real(kind=wp) :: a,b,h,z1,z2

    !write (*,fmt="(1x,a,E16.8)",advance="NO"), "zin=",zin
    !write (*,*) ""

    if ((zin <= bspline%z(1)) .or. (zin >= bspline%z(bspline%n))) then
    ! not in field anymore. return early
        ! write (*,fmt="(1x,a,E16.8,a,E16.8,a,E16.8)",advance="NO"), "zin not in bspline%z(1),bspline%z(bspline%n): zin=",zin, " bspline%z(1)=",bspline%z(1), " bspline%z(bspline%n)=",bspline%z(bspline%n)
        ! write (*,*) ""
        fyeval = 0_wp
        fydeval = 0_wp
        return
    end if
    if ((zin < bspline%z(klo)) .or. (zin > bspline%z(khi))) then
      ! restart bisection - could be solved better, but this is a valid solution
        ! print *, "restarting bisection" 
        klo = 1_ip
        khi = bspline%n
    end if
    ! bisection
    do while ((khi-klo) > 1)
      k = (khi+klo)/2_ip
      if (bspline%z(k) > zin) then
        khi = k
      else
        klo = k
      end if
    end do ! value should now be between 2 z positions of interpolation

    z1 = bspline%z(klo)
    z2 = bspline%z(khi)
    h = z2-z1
    a = (z2-zin)/h
    b = (zin-z1)/h
    fyeval = a*bspline%by(klo)+b*bspline%by(khi) &
               + ((a**3.0_wp-a)*bspline%cy(klo)+(b**3-b)*bspline%cy(khi)) &
               * (h**2.0_wp)/2.0_wp
    fydeval = ( 2.0_wp * (bspline%by(khi)-bspline%by(klo)) &
                - bspline%cy(khi)*(a**2.0_wp + 2.0_wp * a * b - 2.0_wp * b ** 2.0_wp)*h**2.0_wp &
                + bspline%cy(klo)*(b**2.0_wp + 2.0_wp * a * b - 2.0_wp * a ** 2.0_wp)*h**2.0_wp &
              ) / (2.0_wp * h)

    ! slightly increase window again for next iteration.
    khi = MIN(khi,bspline%n-3_ip)+3_ip
    klo = MAX(klo,4_ip)-3_ip
    ! I am no expert for MPI, but found that sometimes negative indices appear here.
    if (klo < 1_ip) then
      klo = 1_ip
    end if
    if (khi > bspline%n) then
      khi = bspline%n
    end if
    return
end subroutine evaluateSplineBfield

end module undfieldinterp
