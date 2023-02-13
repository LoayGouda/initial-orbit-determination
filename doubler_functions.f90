module doubler_functions
    implicit none

    ! ----    Variable declaration    ----
    integer, parameter  :: ikind=SELECTED_REAL_KIND(p=8)     ! Precision declaration

    ! Constants
    real, parameter     :: R_earth  = 6378.1363              ! Earth radius [km]
    real, parameter     :: mu       = 3.986004418e5          ! Earth Gravitaional constant [km3/s2]
    real, parameter     :: flt      = 3.352810665e-3         ! Earth flattening
    real, parameter     :: pctchg   = 0.005                  ! change in each iteration
    real, parameter     :: pi       = 3.141592               ! Mathematical ratio PI
    real, parameter     :: tol      = 1.0e-10                ! tolerance

contains
!--------------------------------------------------------------------
!!
!> breif         This subroutine executes the main Double-R iteration
!! 
!--------------------------------------------------------------------
subroutine doubler_itr(cc1,cc2,magrSite1,magrSite2,magr1in,magr2in,los1,los2,los3,rsite1,rsite2,rsite3,t1,t3, &
    r2,r3,f1,f2,magr1,magr2,a,e,deltae32,rho1,rho2,rho3)

    ! ----   Constants   ----
    real, parameter     :: re  = 6378.1363                   ! Earth radius [km]
    real, parameter     :: mu  = 3.986004418e5               ! Earth Gravitaional constant [km3/s2]

    ! Inputs
    real(KIND = ikind), intent (in)               :: cc1,cc2,magrSite1,magrSite2,magr1in,magr2in,t1,t3       
    real(KIND=ikind), dimension(3), intent (in)   :: rSite1,rSite2,rSite3
    real, dimension(3), intent (in)               :: los1,los2,los3
    ! Outputs
    real(KIND = ikind), intent (out)              :: f1,f2,magr1,magr2,a,e,deltae32
    real(KIND=ikind), dimension(3), intent (out)  :: r2,r3
    ! Auxiliary variables
    real(KIND = ikind), intent (out) :: rho1,rho2,rho3
    real(KIND = ikind)               :: cosdv21,sindv21,dv21,cosdv31,sindv31,dv31,cosdv32,sindv32,dv32
    real(KIND = ikind), dimension(3) :: r1,w
    real(KIND = ikind)               :: c1,c3,p,s,n,c,deltam32,deltam12,magr3
    real(KIND = ikind)               :: ecosv1,esinv1,esinv3,ecosv2,esinv2,ecosv3
    real(KIND = ikind)               :: cosde32,sinde32,cosde21,sinde21,deltae21
    real(KIND = ikind)               :: sindh32,sindh21,deltah32,deltah21
    ! -----------------------
    rho1                = (-cc1 + SQRT(cc1**2-4*(magrSite1**2-magr1in**2))) / 2.0
    rho2                = (-cc2 + SQRT(cc2**2-4*(magrSite2**2-magr2in**2))) / 2.0
    
    r1                  = rho1 * los1 + rSite1
    r2                  = rho2 * los2 + rSite2

    magr1               = NORM2(r1)
    magr2               = NORM2(r2)

    w                   = CROSS(r1,r2)/(magr1*magr2)

    ! Check that L3 is not perpendicular to W
    if ( ABS(DOT_PRODUCT(los3,w)) == 0d0 ) then   
        print *, 'Warning!! rho3 is singular. Please choose another observation t3' 
    else
        rho3                = -DOT_PRODUCT(rSite3,w)/DOT_PRODUCT(los3,w) 
    end if 

    r3                  = rho3 * los3 + rSite3
    magr3               = NORM2(r3)

    ! When rho3 lies in the orbit plane, vector R_3 and L_3 are perpendicular to w, 
    ! and the Equation for rho3 is singular. Should such a singularity occur, 
    ! a different observation time t_3 must be used. Thus the vectors R_1, R_2, and R_3 
    ! have to be determined as functions of the estimated vector magnitudes R1 and R2.
    
    ! The difference in the true anomalies can be determined as follows:
    cosdv21             = DOT_PRODUCT(r2,r1)/(magr2*magr1)
    cosdv31             = DOT_PRODUCT(r3,r1)/(magr3*magr1)
    cosdv32             = DOT_PRODUCT(r3,r2)/(magr3*magr2)
    
    if (w(3) >= 0.0) then
        sindv21         = (r1(1)*r2(2)-r2(1)*r1(2))/ABS(r1(1)*r2(2)-r2(1)*r1(2))*SQRT(1.0-cosdv21**2)
        sindv31         = (r1(1)*r3(2)-r3(1)*r1(2))/ABS(r1(1)*r3(2)-r3(1)*r1(2))*SQRT(1.0-cosdv31**2)
        sindv32         = (r2(1)*r3(2)-r3(1)*r2(2))/ABS(r2(1)*r3(2)-r3(1)*r2(2))*SQRT(1.0-cosdv32**2)
    else if (w(3) < 0.0) then
        sindv21         = -(r1(1)*r2(2)-r2(1)*r1(2))/ABS(r1(1)*r2(2)-r2(1)*r1(2))*SQRT(1.0-cosdv21**2)
        sindv31         = -(r1(1)*r3(2)-r3(1)*r1(2))/ABS(r1(1)*r3(2)-r3(1)*r1(2))*SQRT(1.0-cosdv31**2)
        sindv32         = -(r2(1)*r3(2)-r3(1)*r2(2))/ABS(r2(1)*r3(2)-r3(1)*r2(2))*SQRT(1.0-cosdv32**2)
    end if
    
    dv21                = ATAN2(sindv21,cosdv21)
    dv31                = ATAN2(sindv31,cosdv31)
    dv32                = ATAN2(sindv32,cosdv32)

    ! In order to correct the estimated values of R1 and R2, it is necessary to compute 
    ! the resulting time intervals  between (R3, R2) and (R1, R2) to obtain residuals 
    ! as actual time differences.

    ! To avoid singularity inherent in the semilatus rectum equation due to very short
    ! observational arcs (in this case Gauss Method should be used) we differentiate 
    ! between true anomaly dv31 that is > pi and those =< pi.
    if (dv31 <= pi) then
        c1              = (magr1/magr2)*(sindv31/sindv32)
        c3              = (magr1/magr3)*(sindv21/sindv32)
        p               = (magr1+c3*magr3-c1*magr2)/(1.0+c3-c1)
    else if (dv31 > pi) then
        c1              = (magr2*sindv32)/(magr1*sindv31)
        c3              = (magr2*sindv21)/(magr3*sindv31)
        p               = (c1*magr1+c3*magr3-magr2)/(c1+c3-1.0)
    end if

    ! Conic equation for the true anomalies
    ecosv1              = p/magr1-1.0
    ecosv2              = p/magr2-1.0
    ecosv3              = p/magr3-1.0

    if (dv21 /= pi) then
        esinv1          = (ecosv1*cosdv21-ecosv2)/sindv21
        esinv2          = (-ecosv2*cosdv21+ecosv1)/sindv21
    else
        esinv2          = (ecosv2*cosdv32-ecosv3)/sindv31
        esinv3          = (-ecosv3*cosdv32+ecosv2)/sindv31    
    end if

    e                   = SQRT(ecosv2**2 + esinv2**2)        ! Eccentricity      
    a                   = p/(1-e**2)                         ! Semimajor Axis
    if (e < 1.0) then
        n               = SQRT(mu/a**3)                      ! Mean motion

        s               = magr2/p*SQRT(1-e**2)*esinv2
        c               = magr2/p*(e**2+ecosv2)              ! Should an e be added here?

        ! Differences in Eccentric anomalies
        sinde32         = (magr3/SQRT(a*p))*sindv32 - (magr3/p)*(1-cosdv32)*s
        cosde32         = 1 - magr2*magr3/(a*p)*(1-cosdv32)
        deltae32        = ATAN2(sinde32,cosde32)

        sinde21         = magr1/SQRT(a*p)*sindv21 + (magr1/p)*(1-cosdv21)*s
        cosde21         = 1 - magr2*magr1/(a*p)*(1-cosdv21)
        deltae21        = ATAN2(sinde21,cosde21)

        ! Differences in Mean anomalies using Kepler Equation
        deltam32        = deltae32 + 2*s*(SIN(deltae32/2))**2 - c*SIN(deltae32)
        deltam12        = -deltae21 + 2*s*(SIN(deltae21/2))**2 + c*SIN(deltae21)
    else
        print *, 'Warning!! Orbit is hyperbolic, e = ',e
        print *, ' '
        n               = SQRT(mu/(-a)**3)

        s               = magr2/p*SQRT(e**2-1)*esinv2
        c               = magr2/p*(e**2+ecosv2)

        sindh32         = magr3/SQRT(-a*p)*sindv32-magr3/p*(1-cosdv32)*s
        sindh21         = magr1/SQRT(-a*p)*sindv21+magr1/p*(1-cosdv21)*s

        deltah32        = LOG( sindh32 + SQRT(sindh32**2 + 1) )
        deltah21        = LOG( sindh21 + SQRT(sindh21**2 + 1) )

        deltam32        = -deltah32 + 2*s*(SINH(deltah32/2))**2 + c*SINH(deltah32)
        deltam12        = deltah21 + 2*s*(SINH(deltah21/2))**2 - c*SINH(deltah21)

        deltae32        = deltah32
    end if 

    ! F1 and F2 save the difference between the calculated time difference (deltam12/n)
    ! and the ephermis time difference (t1) corresponding to the station observations.
    ! Thus r1 and r2 will be adjusted to to obtain agreement. 
    f1                  = t1 - deltam12/n
    f2                  = t3 - deltam32/n

    ! The calculated and actual time differences will agree when both f1 and f2 equal 0

end subroutine doubler_itr

!--------------------------------------------------------------------
!!
!> breif         This subroutine calculates the classical orbital 
!!               elements from position and velocity vectors 
!!               (source: Vallado)
!!
!> inputs
!!              r2       - Position vector             km
!!              v2       - Velocity vector             km
!!
!> outputs
!!              sa       - Semi-major axis             km
!!              ecc      - Eccentricity                - 
!!              incl     - Inclination                 rad
!!              raan     - RA of Ascending Node        rad
!!              argp     - Argument of Perigee         rad
!!              nu       - True Anomaly                rad
!!
!--------------------------------------------------------------------
subroutine rv2ceo(r2,v2,sa,ecc,incl,raan,argp,nu)
    ! Inputs    
    real(KIND=ikind), dimension(3), intent (in)   :: r2,v2
    ! Outputs
    real(KIND = ikind), intent (out)              :: sa,ecc,incl,raan,argp,nu  
    ! Auxiliary variables
    real(KIND = ikind)                            :: r2mag,v2mag,magh,magn,c1,rdotv,sme
    real(KIND = ikind)                            :: hk,temp
    real(KIND=ikind), dimension(3)                :: hbar,nbar,ebar
    integer                                       :: i = 1
    character(len=2)                              :: orbittype = 'ei'

    r2mag               = NORM2(r2)
    v2mag               = NORM2(v2)

    ! Angular momentum vector
    hbar                = CROSS(r2,v2)    
    magh                = NORM2(hbar)
    if (magh > tol) then 
        ! Line of Nodes vector
        nbar(1)         = -hbar(2)
        nbar(2)         =  hbar(1)
        nbar(3)         =   0.0

        magn            = NORM2(nbar)

        c1              = v2mag*v2mag - mu /r2mag
        rdotv           = DOT_PRODUCT(r2,v2)

        do i = 1, 3
            ! Eccentricity vector
            ebar(i)     = (c1*r2(i) - rdotv*v2(i))/mu
        end do
        ! -- Eccentricity --
        ecc             = NORM2(ebar)

        sme             = ( v2mag*v2mag*0.5  ) - (mu/r2mag)
        if (ABS(sme) > tol) then
            ! -- Semi-major axis --
            sa          = -mu/(2.0*sme)
        else
            ! Semi-major axis undefined
            sa          = 0.0
            sa          = sa/sa
        end if

        ! -- Inclination --
        hk              = hbar(3)/magh
        incl            = ACOS(hk)

        ! Determine orbit type
        if (ecc < tol) then 
            if (incl < tol .or. ABS(incl-pi) < tol) then
                orbittype = 'ce'            ! Equatorial circular orbit
            else
                orbittype = 'ci'            ! Inclined circular orbit
            end if
        else
            if (incl < tol .or. ABS(incl-pi) < tol) then
                orbittype = 'ee'            ! Equatorial elliptical orbit
            end if
        end if

        ! -- RAAN -- 
        if (magn > tol) then 
            temp        = nbar(1)/magn
            if (ABS(temp) > 1.0) then
                if (temp < 0.0) then
                    temp = -1.0
                else if (temp > 0.0) then
                    temp = 1.0
                else 
                    temp = 0.0
                end if
            end if
            raan        = ACOS(temp)
            if (nbar(2) < 0.0) then
                raan    = 2.0*pi - raan
            end if
        end if

        ! -- Argument of Perigee --
        if (orbittype == 'ei') then
            argp        = angl2v(nbar,ebar)
            if (ebar(3) < 0.0) then
                argp    = 2.0*pi - argp
            end if
        else
            argp        = 0.0 
            argp        = argp/argp
        end if 

        ! -- True Anomaly -- 
        if (orbittype == 'ei' .or. orbittype == 'ee') then
            nu          = angl2v(ebar,r2)
            if (rdotv < 0.0) then
                nu      = 2.0*pi - nu
            end if
        else
            nu          = 0.0
            nu          = nu/nu
        end if
    else
        ! Undefined variables
        sa              = 0.0
        sa              = sa/sa 
        ecc             = 0.0
        ecc             = ecc/ecc 
        incl            = 0.0
        incl            = incl/incl
        raan            = 0.0
        raan            = raan/raan
        argp            = 0.0
        argp            = argp/argp
        nu              = 0.0
        nu              = nu/nu

    end if

end subroutine rv2ceo
!--------------------------------------------------------------------
!!
!> breif         This function calculates observation site 
!!               position vector from geodetic coordinats.
!!
!> inputs
!!              lst      - Local Sidereal Time         degrees
!!              lon      - East longitude              degrees
!!              alt      - Altitude over N.N           meters
!!
!> outputs
!!              lst      - Local sidereal time         degrees
!! 
!--------------------------------------------------------------------
!function site_vector(lst,alt,lat) result(rSite)
!    real(KIND=ikind),dimension(3,3) :: pm
!    real(KIND=ikind),dimension(3)   :: rSite
!    real(KIND=ikind)                :: lst,lat,rlst,rlat
!    real(KIND=ikind)                :: xp,yp
!    real(KIND=ikind)                :: cosxp,cosyp,sinxp,sinyp
!    real                            :: alt

    ! Consider Polar Motion
!    xp    = 0.137495 * pi / (180*3600.0)
!    yp    = 0.342416 * pi / (180*3600.0)
!    cosxp = cos(xp)
!    sinxp = sin(xp)
!    cosyp = cos(yp)
!    sinyp = sin(yp)

!    pm(1,1) =  cosxp
!    pm(1,2) =  0.d0
!    pm(1,3) = -sinxp
!    pm(2,1) =  sinxp * sinyp
!    pm(2,2) =  cosyp
!    pm(2,3) =  cosxp * sinyp
!    pm(3,1) =  sinxp * cosyp
!    pm(3,2) = -sinyp
!    pm(3,3) =  cosxp * cosyp

    ! Turn degrees to radiens 
!   rlst                = lst*pi/180
!    rlat                = lat*pi/180

    ! Calculate Site Vector in ECI (Curtis Eq 5.56)
!    rSite(1)            = (R_earth/SQRT(1-(2*flt - flt**2)*SIN(rlat)**2) + alt)*COS(rlat)*COS(rlst)
!    rSite(2)            = (R_earth/SQRT(1-(2*flt - flt**2)*SIN(rlat)**2) + alt)*COS(rlat)*SIN(rlst)
!    rSite(3)            = ((R_earth*(1 - flt)**2)/SQRT(1-(2*flt - flt**2)*SIN(rlat)**2) + alt)*SIN(rlat)
    
!    rSite               = MATMUL(pm,rSite)
!    print *, rSite
    
!end function site_vector

function site_vector(lst,alt,lat) result(rSite)
    real(KIND=ikind),dimension(3,3) :: pm
    real(KIND=ikind),dimension(3)   :: rSite
    real(KIND=ikind)                :: lst,lat,alt,rlst,rlat
    real(KIND=ikind)                :: xp,yp
    real(KIND=ikind)                :: cosxp,cosyp,sinxp,sinyp


    ! Turn degrees to radiens 
    rlst                = lst*pi/180.0_ikind
    rlat                = lat*pi/180.0_ikind

    ! Consider Polar Motion
    xp    = 0.137495 * pi / (180*3600.0)
    yp    = 0.342416 * pi / (180*3600.0)
    cosxp = cos(xp)
    sinxp = sin(xp)
    cosyp = cos(yp)
    sinyp = sin(yp)

    pm(1,1) =  cosxp
    pm(1,2) =  0.d0
    pm(1,3) = -sinxp
    pm(2,1) =  sinxp * sinyp
    pm(2,2) =  cosyp
    pm(2,3) =  cosxp * sinyp
    pm(3,1) =  sinxp * cosyp
    pm(3,2) = -sinyp
    pm(3,3) =  cosxp * cosyp

    ! Calculate Site Vector in ECI (Curtis Eq 5.56)
    rSite(1)            = (R_earth/SQRT(1.d0-(2.d0*flt - flt**2)*SIN(rlat)**2) + alt)*COS(rlat)*COS(rlst)
    rSite(2)            = (R_earth/SQRT(1.d0-(2.d0*flt - flt**2)*SIN(rlat)**2) + alt)*COS(rlat)*SIN(rlst)
    rSite(3)            = ((R_earth*(1.d0 - flt)**2)/SQRT(1.d0-(2.d0*flt - flt**2)*SIN(rlat)**2) + alt)*SIN(rlat)
    
    !rSite               = MATMUL(pm,rSite)
    print *, rSite
    
end function site_vector
!--------------------------------------------------------------------
!!
!> breif         This function calculates the local sidereal time.
!!
!> inputs
!!              JD       - Julian Date at 0 hr UT
!!              EL       - East longitude              degrees
!!              ut       - Universal Time              hours
!!
!> outputs
!!              lst      - Local sidereal time         degrees
!! 
!--------------------------------------------------------------------
function sidereal_time(JD,ut,lon) result(lst)
REAL(KIND = ikind) :: JD,j0,g0,gst,lst,ut,lon

j0                  = (JD - 2451545)/36525
g0                  = 100.4606184 + 36000.77004*j0 + 0.000387933*(j0**2) - 2.583e-8*(j0**3)

if (g0 .ge. 360) then
    g0 = g0 - INT(g0/360)*360
else if (g0 .lt. 0) then
    g0 = g0 - (INT(g0/360) - 1)*360
else
    g0 = g0
end if

gst                 = g0 + 360.98564724*ut/24 

lst                 = gst + lon

lst                 = lst - 360*INT(lst/360)

end function sidereal_time

!--------------------------------------------------------------------
!!
!> breif         This function calculates universal time UT 
!!               from a given julian date.
!!
!> inputs
!!              JD       - Julian Date at 0 hr UT
!!
!> outputs
!!              ut       - UT at given JD              hours
!! 
!--------------------------------------------------------------------
function universal_time(JD) result(ut)
REAL(KIND = ikind) :: JD,diff,ut

diff                = JD - INT(JD)

if (0.5 <= diff .and. diff <= 0.999) then
    ut              = (diff*24.0_ikind)-12.0_ikind
else if (0.0 <= diff .and. diff < 0.5) then
    ut              = (diff*24.0_ikind)+12.0_ikind
else if (diff .eq. 0) then
    ut              = 12.0_ikind
end if


end function universal_time

!--------------------------------------------------------------------
!!
!> breif         This function calculates the cross product 
!!               between two given vectors.
!!
!> inputs
!!              a        - First vector
!!              b        - Second vector
!!
!--------------------------------------------------------------------
function CROSS(a, b)
    real(KIND=ikind), dimension(3)                :: cross
    real(KIND=ikind), dimension(3), intent (in)   :: a, b
      
    cross(1)            = a(2) * b(3) - a(3) * b(2)
    cross(2)            = a(3) * b(1) - a(1) * b(3)
    cross(3)            = a(1) * b(2) - a(2) * b(1)

end function CROSS

!--------------------------------------------------------------------
!!
!> breif         This subroutine transforms degrees to radians and 
!!               radians to degrees
!!
!> inputs
!!              angl1    - Angle after conversion      deg/rad
!!              angl2    - Angle before conversion     deg/rad 
!!              typ      - Type to convert to             char
!!                         'rad' for radians
!!                         'deg' for degrees
!!
!--------------------------------------------------------------------
subroutine degrad(angl1,angl2,typ)
    real(KIND=ikind), intent (inout) :: angl1
    real(KIND=ikind), intent (inout) :: angl2
    character(len=*)                 :: typ
    
    if (typ == 'rad') then
        angl1               = angl2*pi/180
    else if (typ == 'deg') then
        angl1               = angl2*180/pi
    end if         

end subroutine degrad

!--------------------------------------------------------------------
!!
!> breif         This functions calculates the angle 
!!               between two 3D vectors
!!
!> inputs
!!              vec1     - First vector
!!              vec2     - Second vector 
!!
!> output
!!              theta    - Angle between both bectors [rad]
!!
!--------------------------------------------------------------------
function angl2v(vec1,vec2) result(theta)
    real(KIND = ikind),dimension(3) :: vec1,vec2
    real(KIND = ikind)              :: theta

    theta = ACOS(DOT_PRODUCT(vec1,vec2)/(SQRT(vec1(1)**2+vec1(2)**2+vec1(3)**2) & 
    *SQRT(vec2(1)**2+vec2(2)**2+vec2(3)**2)))

end function angl2v

end module doubler_functions