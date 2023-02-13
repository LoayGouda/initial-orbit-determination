!-------------------------------------------------------------------------------
!!
!> breif         This program solves the problem of orbit determination from
!!               three optical sightings using the Double-r iteration method
!!
!> inputs
!!              ra1      - right ascension at t1       rad
!!              ra2      - right ascension at t2       rad
!!              ra3      - right ascension at t3       rad
!!              dec1     - declination at t1           rad
!!              dec2     - declination at t2           rad
!!              dec3     - declination at t3           rad
!!              Mjd1     - Modified julian date of t1
!!              Mjd2     - Modified julian date of t2
!!              Mjd3     - Modified julian date of t3
!!              rSite1   - ijk site1 position vector   km
!!              rSite2   - ijk site2 position vector   km
!!              rSite3   - ijk site3 position vector   km
!!
!> outputs
!!              r        - ijk position vector at t2   km
!!              v        - ijk velocity vector at t2   km/s
!-------------------------------------------------------------------------------

program doubler

    use doubler_functions 
 
    implicit none

    ! Observation Site Location Data (From Vallado)
    real(KIND=ikind), parameter    :: lon      = -110.0      ! Longitude [deg]
    real(KIND=ikind), parameter    :: lat      = 40.0        ! Latitude  [deg]
    real(KIND=ikind), parameter    :: alt      = 2.0         ! Altitude above N.N [km]
    real(KIND=ikind), dimension(3) :: rSite1,rSite2,rSite3   ! Site Position Vector [ECI]

    real                :: ra1,ra2,ra3,dec1,dec2,dec3        ! Observation Data (RA & Dec)
    real(KIND = ikind)  :: Mjd1,Mjd2,Mjd3                    ! Modified Julian Date
    real(KIND = ikind)  :: JD1,JD2,JD3                       ! Standard Julian Date

    ! Auxiliary variables
    real(KIND = ikind)  :: magr1in,magr2in,magrSite1,magrSite2,magrSite3
    real(KIND = ikind)  :: ut1,ut2,ut3,t1,t3,sa,ecc,incl,raan,argp,nu
    real(KIND = ikind)  :: lst1,lst2,lst3,rho1,rho2,rho3    ! Local sidereal time of site at each Observation
    real, dimension(3)  :: los1,los2,los3                    ! Line-of-Sight vectors for each Obs.
    real(KIND = ikind)  :: cc1,cc2,f1_in,f2_in,f,g
    real(KIND = ikind)  :: f1delr1,f2delr1,pf1pr1,pf2pr1,f1delr2,f2delr2,pf1pr2,pf2pr2
    real(KIND = ikind)  :: f1,f2,magr1,magr2,a,e,deltae32    ! double-r iteration subroutine output
    real(KIND = ikind), dimension(3) :: r2,r3,v2
    real(KIND = ikind)  :: delta,delta1,delta2,deltar1,deltar2,magr1o,magr2o
    integer             :: ktr
    character(len=5)    :: Flag                 = 'True'
    !---------------------------------------
    
    print *, ' '
    print *, 'Initializing Double-R iteration routine ...'
    print *, ' '  

    ! ----    Input Data    ----

    ! Observation elements
    ra1                 = 6.115312122058        
    ra2                 = 1.93813490179       
    ra3                 = 2.73515399095           
    dec1                = 0.217056144548      
    dec2                = 0.355574026099       
    dec3                = -0.222179139688   

    ! Observation times
    Mjd1                = 56159.64418359_ikind
    Mjd2                = 56159.66917824_ikind
    Mjd3                = 56159.69070602_ikind

    print *, 'For Object with RA and Declination of:'
    print *, ' RA                       Dec'
    print *, ra1, '      ',dec1
    print *, ra2, '      ',dec2
    print *, ra3, '      ',dec3
    print *, ' '
    print *, 'Data acquired from Observation site at:'
    print *, ' '
    print *, 'Longitude: ',lon,' degrees' 
    print *, 'Latitude: ',lat,' degrees'
    print *, 'Altitude: ',alt,'     km above N.N'
    print *, ' '

    ! Transform Modified Julian Date to standard Julian Date
    JD1                 = Mjd1 + 2400000.5_ikind
    JD2                 = Mjd2 + 2400000.5_ikind
    JD3                 = Mjd3 + 2400000.5_ikind
    
    print *, 'Calculating UT for each observation ...'
    print *, ' '
    ut1                 = universal_time(JD1)
    ut2                 = universal_time(JD2) 
    ut3                 = universal_time(JD3)
    print *, '1st Obs. UT = ', ut1, ' h'
    print *, '2nd Obs. UT = ', ut2, ' h'
    print *, '3rd Obs. UT = ', ut3, ' h'
    print *, ' '

    print *, 'Calculating Local Sidereal Time for each observation ...'
    print *, ' '
    ! Local Sidereal Time of site at each JD
    lst1                = sidereal_time(JD1,ut1,lon)
    lst2                = sidereal_time(JD2,ut2,lon)
    lst3                = sidereal_time(JD3,ut3,lon)
    print *, '1st Obs. = ', lst1, ' degrees'
    print *, '2nd Obs. = ', lst2, ' degrees'
    print *, '3rd Obs. = ', lst3, ' degrees'
    print *, ' '

    ! Calculate the Site position vector in ECI coordinates
    print *, 'Site Position Vector is:'
    print *, ' '
    rSite1              = site_vector(lst1,alt,lat) 
    rSite2              = site_vector(lst2,alt,lat)
    rSite3              = site_vector(lst3,alt,lat)
    print *, ' '

    ! Variables preparations
    t1                  = (JD1 - JD2)*86400.0
    t3                  = (JD3 - JD2)*86400.0

    ! Form the Line-Of-Sight vectors using RA and Dec (radians)
    los1                = (/COS(dec1)*COS(ra1), COS(dec1)*SIN(ra1), SIN(dec1)/)
    los2                = (/COS(dec2)*COS(ra2), COS(dec2)*SIN(ra2), SIN(dec2)/)
    los3                = (/COS(dec3)*COS(ra3), COS(dec3)*SIN(ra3), SIN(dec3)/)
    print *, 'Line-Of-Sight Vectors are: '
    print *, ' '
    print *, 'L1 = ',los1
    print *, 'L2 = ',los2
    print *, 'L3 = ',los3
    print *, ' '

    ! Begin actual Double-R algorithm  
    magr1in             = R_earth*3.0_ikind         ! Form initial guesses
    magr2in             = R_earth*3.01_ikind        ! for r1 and r2 

    magrSite1           = NORM2(rSite1)
    magrSite2           = NORM2(rSite2)
    magrSite3           = NORM2(rSite3)

    cc1                 = 2.0*DOT_PRODUCT(los1,rSite1)
    cc2                 = 2.0*DOT_PRODUCT(los2,rSite2)
    ktr                 = 0

    print*, 'Interval = ',(ABS(t1)+ABS(t3))/60.0
    print*, ' '

    ! Main loop for the iteration
    do while  (Flag == 'True')
        ktr             = ktr + 1

        CALL doubler_itr(cc1,cc2,magrSite1,magrSite2,magr1in,magr2in,los1,los2,los3,rsite1,rsite2,rsite3,t1,t3, &
            r2,r3,f1,f2,magr1,magr2,a,e,deltae32,rho1,rho2,rho3)
    
        f1_in          = f1                         ! Without incrementation of r1 or r2
        f2_in          = f2

        ! Increment r1 by 0.005 and repeat iteration (with r1 = r1 + delta r1)
        magr1o         = magr1in                    ! Save value of r1 before incrementation
        deltar1        = pctchg*magr1in             
        magr1in        = (1.0 + pctchg)*magr1in     ! Increment r1 by 0.005
        ! Call subroutine again after first incrementation
        CALL doubler_itr(cc1,cc2,magrSite1,magrSite2,magr1in,magr2in,los1,los2,los3,rsite1,rsite2,rsite3,t1,t3, &
        r2,r3,f1,f2,magr1,magr2,a,e,deltae32,rho1,rho2,rho3)

        f1delr1        = f1                         ! F1(r1+del(r),r2)
        f2delr1        = f2                         ! F2(r1+del(r),r2)

        ! Increment r2 by 0.005 and repeat iteration (with r2 = r2 + delta r2)
        magr1in        = magr1o                     ! Reset r1 to initial value
        deltar2        = pctchg*magr2in
        magr2o         = magr2in                    ! Save value of r2 before incrementation
        magr2in        = (1.0 + pctchg)*magr2in     ! Increment r2 by 0.005
        ! Call subroutine again after second incrementation
        CALL doubler_itr(cc1,cc2,magrSite1,magrSite2,magr1in,magr2in,los1,los2,los3,rsite1,rsite2,rsite3,t1,t3, &
        r2,r3,f1,f2,magr1,magr2,a,e,deltae32,rho1,rho2,rho3)

        magr2in         = magr2o                    ! Reset r2 to initial value
        
        f1delr2         = f1                        ! F1(r1,r2+del(r))
        f2delr2         = f2                        ! F2(r1,r2+del(r))

        pf1pr1         = (f1delr1-f1_in)/deltar1
        pf2pr1         = (f2delr1-f2_in)/deltar1

        pf1pr2          = (f1delr2-f1_in)/deltar2
        pf2pr2          = (f2delr2-f2_in)/deltar2

        delta           = pf1pr1*pf2pr2 - pf2pr1*pf1pr2
        delta1          = pf1pr2*f2_in - pf2pr2*f1_in
        delta2          = pf2pr1*f1_in - pf1pr1*f2_in

        deltar1         = delta1/delta
        deltar2         = delta2/delta

        if (ABS(deltar1) < tol .and. ABS(deltar2) < tol) then
            Flag        = 'False'
        else
            magr1in     = magr1in + deltar1
            magr2in     = magr2in + deltar2
            Flag        = 'True'
        end if 

        ! Limit Iterations in case loop breaks
        if (ktr > 20) then
            Flag        = 'False'
        end if
        
    end do

    CALL doubler_itr(cc1,cc2,magrSite1,magrSite2,magr1in,magr2in,los1,los2,los3,rsite1,rsite2,rsite3,t1,t3, &
        r2,r3,f1,f2,magr1,magr2,a,e,deltae32,rho1,rho2,rho3)
        
    f                   = 1.0_ikind - a/magr2*(1.0_ikind-cos(deltae32))
    g                   = t3 - SQRT(a**3/mu)*(deltae32-sin(deltae32))
    v2                  = (r3 - f*r2)/g

    ! Print solution to the command window
    print *, 'After a total of ', ktr,' iterations performed'
    print *, ' '
    Print *, 'The final estimate of the state vectors is:'
    print *, ' '
    print *, 'r2 = ',r2
    print *, 'v2 = ',v2
    print *, ' '
    print *, 'magr2 = ',NORM2(r2),' km'
    print *, 'magv2 = ',NORM2(v2),' km/s'
    print *, ' '

    CALL rv2ceo(r2,v2,sa,ecc,incl,raan,argp,nu)

    print *, 'The classical Orbital Elements are: '
    print *, 'a    : ',sa, 'e    : ',ecc
    print *, 'i    : ',incl*180.0_ikind/pi,'RAAN : ',raan*180.0_ikind/pi
    print *, 'w    : ',argp*180.0_ikind/pi,'v    : ',nu*180.0_ikind/pi
    print *, ' '

end program