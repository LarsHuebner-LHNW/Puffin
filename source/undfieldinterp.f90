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

subroutine read_planepolefield(byfile, byfield)


    implicit none

! byfile        -      file to load
! byfield contains
! n     -      number of datapoints in file
! z               -      positions array
! by              -      by field array
! c              -      coefficients for spline


    character(1024_ip), intent(in) :: byfile
    type(byfieldspline), intent(out) :: byfield

    int(kind=ip) :: ndatapoints

    integer(kind=ip)   :: ios

    open(73,FILE=byfile, IOSTAT=ios, STATUS='OLD', ACTION='READ', POSITION ='REWIND')

    if (ios /= 0) then
        print*, 'iostat = ', ios, "in loadByField"
        stop "OPEN(byfile) not performed correctly, IOSTAT /= 0"
    end if

    ! read the number of datapoints
    read(73,*, IOSTAT=ios) ndatapoints ! number
    if (ios < 0) then ! end of file
        if ((tProcInfo_G%qroot) .and. (ioutInfo_G > 2)) print*, "Reached end of byfield file at line 1!"
        stop "byfield file empty!"
    end if
    byfield%n = ndatapoints

    allocate(byfield%z(ndatapoints))
    allocate(byfield%by(ndatapoints))

    read(73,*,IOSTAT=ios) byfield%z
    if (ios < 0) then ! end of file
        if ((tProcInfo_G%qroot) .and. (ioutInfo_G > 2)) print*, "Reached end of byfield file at line 2!"
        stop "byfield file has only 1 line!"
    end if

    read(73,*,IOSTAT=ios) byfield%by
    if (ios < 0) then ! end of file
        if ((tProcInfo_G%qroot) .and. (ioutInfo_G > 2)) print*, "Reached end of byfield file at line 3!"
        stop "byfield file has only 2 lines!"
    end if

    close(73, STATUS='KEEP')

    allocate(byfield%c(ndatapoints))
    call splineCoeff(ndatapoints,byfield%z,byfield%by,byfield%c)
    return
end subroutine read_planepolefield

subroutine splineCoeff(ndatapoints,z,by,coeff)
    ! calculate a cubic spline
    ! adapted from
    ! NUMERICAL RECIPES IN FORTRAN 77: THE ART OF SCIENTIFIC COMPUTING (ISBN 0-521-43064-X)
    ! page 109
    
    ! coefficients are second derivative at that point
    integer(kind=ip), intent(in) :: ndatapoints
    real(kind=wp), intent(in) :: z(:)
    real(kind=wp), intent(in) :: by(:)
    real(kind=wp), intent(out) :: coeff(:)

    if (ndatapoints > 65536_ip) then! maximum interpolation points. Probably gets too slow if too large
        stop "more than 2**16 (65536) interpolation points. Aborting."
    end if 

    ! helper variables
    real(kind=wp) :: u(ndatapoints)
    real(kind=wp) :: sig, p
    integer(kind=ip) :: i,k

    ! "natural" boundaries
    coeff(1)=0.0_wp
    coeff(ndatapoints)=0.0_wp
    u(1)=0.0_wp

    do i=2,ndatapoints-1
      sig = (z(i)-z(i-1))/(z(i+1)-z(i-1))
      p = sig*coeff(i-1)+2.0_wp
      coeff(i)=(sig-1.0_wp)/p
      u(i) = (6.0_wp* & 
               ( &
                 (by(i+1)-by(i)) / (z(i+1)-z(i)) &
                 -(by(i)-by(i-1)) / (z(i)-z(i-1)) &
               ) &
               / (z(i+1)-z(i-1)) &
               - sig*u(i-1) &
             ) / p
    end do
    do k=ndatapoints-1,1,-1
      coeff(k)=coeff(k)*coeff(k+1)+u(k)
    end do
    return
end subroutine splineCoeff

subroutine evaluateSplineBfield(byspline,zin,klo,khi,feval,fdeval)
    ! return value (By) and first derivative (Bz) of spline
    ! adapted from
    ! NUMERICAL RECIPES IN FORTRAN 77: THE ART OF SCIENTIFIC COMPUTING (ISBN 0-521-43064-X)
    ! page 110
    type(byfieldspline), intent(in) :: byspline
    real(kind=wp), intent(in) :: zin
    real(kind=wp), intent(inout) :: klo ! interval guess low
    real(kind=wp), intent(inout) :: khi ! interval guess high
    real(kind=wp), intent(out) :: feval,fdeval

    ! helper variables
    integer(kind=ip) :: k ! position
    real(kind=wp) :: a,b,h,z1,z2

    if ((zin < byspline%z(klo)) .or. (zin > byspline%z(khi))) then
      ! restart bisection
        klo = 1
        khi = byspline%n
    end if
    ! bisection
    do while ((khi-klo) > 1)
      k = (khi+klo)/2
      if (byspline%z(k) > zin) then
        khi = k
      else
        klo = k
      end if
    end do ! value should now be between 2 z positions of interpolation

    z1 = byspline%z(klo)
    z2 = byspline%z(khi)
    h = z2-z1
    a = (z2-zin)/h
    b = (zin-z1)/h
    feval = a*byspline%by(klo)+b*byspline%by(khi) &
               + ((a**3.0_wp-a)*byspline%c(klo)+(b**3-b)*byspline%c(khi)) &
               * (h**2.0_wp)/2.0_wp
    fdeval = (byspline%by(khi)-byspline%by(klo))*zin/h &
                + ( &
                     (h**2.0_wp - 3.0_wp * (z2-zin)**2.0_wp)*byspline%c(klo) &
                   - (h**2.0_wp - 3.0_wp * (zin-z1)**2.0_wp)*byspline%c(khi) &
                ) &
                / (2.0_wp*h) 

    ! slightly increase window again for next run. 
    khi = MIN(khi,byspline%n-5)+5
    klo = MIN(klo,6)-5
    return
end subroutine splineCoeff

end module undfieldinterp
