module diffusion_module
  implicit none
  ! Define a double precision kind constant to reuse, for consistency in precision.
  integer, parameter :: dp = kind(1.0d0)  
  integer, parameter :: num_levels = 100
  ! Physical constants
  real(dp), parameter :: BOLTZMANN_K = 1.3806488d-23            ! J/K
  real(dp), parameter :: BIG_G = 6.67408d-11                  ! N m^2/kg^2  (gravitational constant)
  real(dp), parameter :: M_H = 1.6726219d-27                    ! kg (mass of a proton)
  real(dp), parameter :: M_C = 1.99d-26                    ! kg (mass of a carbon atom)
  real(dp), parameter :: M_O = 2.656d-26                    ! kg (mass of an oxygen atom)
  real(dp), parameter :: M_N = 2.325d-26                    ! kg (mass of a nitrogen atom)
  real(dp), parameter :: R_gas = 8.314472d0             ! N m K^-1 mol^-1  (gas constant)
  real(dp), parameter :: M_CO2 = 0.04401d0                 ! kg/mol (molar mass of CO2)
  ! real(dp), parameter :: alpha_T = -0.25d0                 ! Thermal diffusion factor of H, H2 (from Krasnopolsky 2002)
  real(dp), parameter :: alpha_T = 0.0d0                   ! others
    ! real(dp), parameter :: D0 = 8.4d0, s = 0.597d0 ! H
  ! real(dp), parameter :: D0 = 2.23d0, s = 0.75d0 ! H2
  real(dp), parameter :: D0 = 1.0d0, s = 0.75d0 ! others
  ! real(dp), parameter :: pi = 4.0d0*DATAN(10.d0)       
  real(dp), parameter :: pi = 3.141592653589793d0 
  ! Mars parameters
  real(dp), parameter :: MARS_MASS = 0.1075d0 * 5.972d24
  real(dp), parameter :: MARS_RADIUS = 3389.0d5        

contains

  !-----------------------------------------------------------------------
  !> a piecewise function for temperature as a function of altitude,
  !> using Krasnopolsky's 2010 temperatures for altitudes
  !> >htropo=90km, fixed at Ttropo=125K between htropo and
  !> htropo-htropowidth=60km, and rising at a constant lapse rate
  !> (1.4K/km) below.
  !-----------------------------------------------------------------------
  function temperature_profile(z, Tinf) result(T)
    real(dp), intent(in) :: z, Tinf
    real(dp) :: T
    real(dp), parameter :: Ttropo = 130.0d0 !110.0d0
    real(dp), parameter :: ztropo = 108.0d5 !60.0d5
    real(dp), parameter :: zwidth = 30.0d5
    real(dp), parameter :: lapse = -1.4d-5
    real(dp), parameter :: p4 = 8.0d10

    ! if (z >= ztropo) then
    !    T = Tinf - (Tinf - Ttropo) * exp(-((z - ztropo)**2) / (8.0d10 * Tinf))
    ! else if ((z >= ztropo - zwidth) .and. (ztropo > z)) then
    !    T = Ttropo
    ! else if (ztropo - zwidth > z) then
    !    T = Ttropo - lapse * (ztropo - zwidth - z)
    ! end if

    T = Tinf - (Tinf - Ttropo) * exp(-((z - ztropo)**2) / p4)
  end function temperature_profile

  function Tcustom(z, p) result(T)
    real(dp), intent(in) :: z
    real(8), intent(in)         :: p(4)
    real(dp) :: T

    ! T = p(3) - (p(3) - p(1)) * exp(-((z - p(2))**2) / p(4))
    T = p(3) - (p(3) - p(1)) * exp(-((z - p(2))**2) / (p(4)*p(3)))
  end function Tcustom


  !-----------------------------------------------------------------------
  !> https://ntrs.nasa.gov/citations/19960042695
  !-----------------------------------------------------------------------
  function temperature_MarsGRAM(z, F107, R) result(T)
    real(dp), intent(in) :: z, F107, R   ! z in cm, R=heliocentric distance of Mars in AU
    real(dp) :: T
    real(dp) :: R0, TINF, ZF, TF, SCALE, Zkm, ZZF,RF,YSC

    R0 = 3395.428d0
    TINF = 156.3d0 + 0.9427d0*F107
    ZF = 197.94d0 - 49.058d0*R + 8 ! + 8km for Viking Lander 1
    TF = 113.7d0 + 0.5791d0*F107
    SCALE = 8.38d0 + 0.09725*F107

    Zkm = z*1.0d-5 ! cm to km

    ZZF = Zkm - ZF
    RF = R0 + ZF
    YSC  = ZZF * RF / (RF + ZZF)
    T = TINF - (TINF - TF) * exp(-YSC / SCALE)
  end function temperature_MarsGRAM

  !-----------------------------------------------------------------------
  !> https://www.grc.nasa.gov/www/k-12/airplane/atmosmrm.html
  !-----------------------------------------------------------------------
  function temperature_NASA_MAM(z) result(T)
    real(dp), intent(in) :: z
    real(dp) :: T
    T = -23.4 - 0.00222 * z*1.0d-2 + 274.15
  end function temperature_NASA_MAM

  !-----------------------------------------------------------------------
  !> Initializes the altitude grid.
  !> Parameters:
  !>   zb: Bottom altitude in cm.
  !>   zt: Top altitude in cm.
  !>   n: Number of grid levels.
  !>
  !> Returns:
  !>   altitudes: Array of altitude values (cm).
  !>   dz: Grid spacing (cm).
  !-----------------------------------------------------------------------
  subroutine initialize_grid(zb, zt, n, altitudes, dz)
    ! Input arguments
    real(dp), intent(in) :: zb, zt
    integer, intent(in) :: n
    ! Outout
    real(dp), dimension(n), intent(out) :: altitudes
    real(dp), intent(out) :: dz
    ! Local
    integer :: i

    dz = (zt - zb) / real(n - 1, dp)
    altitudes(1) = zb
    
    do i = 2, n
       altitudes(i) = altitudes(i-1) + dz
    end do

  end subroutine initialize_grid

  !-----------------------------------------------------------------------
  !> Computes the scale height at a given altitude.
  !-----------------------------------------------------------------------
  subroutine scale_height_profile(z, T, H)
    real(dp), dimension(:), intent(in) :: z, T
    real(dp), dimension(:), intent(out) :: H
    integer :: i
    real(8) :: g_local

    ! do i = 1, size(z)
    !     ! Calculate local gravity, converting radii from cm to meters
    !     g_local = BIG_G * MARS_MASS / ( ((z(i) + MARS_RADIUS)*1.0d-2)**2 )
    !     ! Calculate H in meters, then convert to cm for output
    !     H(i) = ( BOLTZMANN_K * T(i) / ( 1.0d0 * M_H * g_local ) ) * 1.0d2
    ! end do
    do i = 1, size(z)
    !   ! For H
      !  H(i) = (BOLTZMANN_K * T(i) / (1.0d0 * M_H * BIG_G * MARS_MASS) * ((z(i) + MARS_RADIUS)*1.0d-2)**2) * 1.0d2
    !   ! For H2
      ! H(i) = (BOLTZMANN_K * T(i) / ((M_H * 2.0d0) * BIG_G * MARS_MASS) * ((z(i) + MARS_RADIUS)*1.0d-2)**2) * 1.0d2
    !   ! For CO2
      !  H(i) = (BOLTZMANN_K * T(i) / ((M_C + M_O * 2.0d0) * BIG_G * MARS_MASS) * ((z(i) + MARS_RADIUS)*1.0d-2)**2) * 1.0d2
    !   ! For N2
      !  H(i) = (BOLTZMANN_K * T(i) / ((M_N * 2.0d0) * BIG_G * MARS_MASS) * ((z(i) + MARS_RADIUS)*1.0d-2)**2) * 1.0d2
    !   ! For O
      !  H(i) = (BOLTZMANN_K * T(i) / ((M_O * 1.0d0) * BIG_G * MARS_MASS) * ((z(i) + MARS_RADIUS)*1.0d-2)**2) * 1.0d2
    !   ! For O2
       H(i) = (BOLTZMANN_K * T(i) / ((M_O * 2.0d0) * BIG_G * MARS_MASS) * ((z(i) + MARS_RADIUS)*1.0d-2)**2) * 1.0d2

    !   ! For H2
      ! H(i) = (BOLTZMANN_K * T(i) / (3.34d-27 * BIG_G * MARS_MASS) * ((z(i) + MARS_RADIUS)*1.0d-2)**2) * 1.0d2
    !   ! For CO2
    !   !  H(i) = (BOLTZMANN_K * T(i) / (44.0d0 * M_H * BIG_G * MARS_MASS) * ((z(i) + MARS_RADIUS)*1.0d-2)**2) * 1.0d2
    !   ! For N2
    !   !  H(i) = (BOLTZMANN_K * T(i) / (2.0d0 * 2.325d-26 * BIG_G * MARS_MASS) * ((z(i) + MARS_RADIUS)*1.0d-2)**2) * 1.0d2
    !   ! For O
      !  H(i) = (BOLTZMANN_K * T(i) / (2.656d-26 * BIG_G * MARS_MASS) * ((z(i) + MARS_RADIUS)*1.0d-2)**2) * 1.0d2
    !   ! For O2
      !  H(i) = (BOLTZMANN_K * T(i) / (2.0d0 * 2.656d-26 * BIG_G * MARS_MASS) * ((z(i) + MARS_RADIUS)*1.0d-2)**2) * 1.0d2
    end do
  end subroutine scale_height_profile

  !-----------------------------------------------------------------------
  !> Forward difference to compute n_CO2 profile (approximated total density in atmosphere)
  !> using input n_CO2 (cm^-3) at top of GCM
  !-----------------------------------------------------------------------
  subroutine total_density_CO2(altitudes, dz, n_CO2_bottom, T_profile, n_CO2)
    ! Input arguments
    real(dp), dimension(:), intent(in) :: altitudes, T_profile
    real(dp), intent(in) :: dz, n_CO2_bottom
    ! Output
    real(dp), dimension(size(altitudes)), intent(out) :: n_CO2
    ! Local variables
    integer :: i
    real(dp) :: g, T_i, dT, R_spec
    ! Specific gas constant
    R_spec = R_gas / M_CO2
    ! Set initial value
    n_CO2(1) = n_CO2_bottom
    ! Compute profile using forward difference
    do i = 2, size(altitudes)
      g = BIG_G * MARS_MASS / ((altitudes(i) + MARS_RADIUS) * 1.0d-2)**2
      T_i = T_profile(i-1)
      dT = T_profile(i) - T_i
      n_CO2(i) = n_CO2(i-1) * (1.d0 - dT / T_i - g * dz * 1.0d-2 / (R_spec * T_i))
    end do
  end subroutine total_density_CO2

  !-----------------------------------------------------------------------
  !> Computes the molecular diffusion coefficient D(z) for Hydrogen at each altitude.
  !>
  !> Parameters:
  !>     T: Temperature profile (K).
  !>     n_tot: Total number density (cm^-3).
  !>
  !> Returns:
  !>     D: Molecular diffusion coefficient D(z) (cm^2/s).
  !-----------------------------------------------------------------------
  subroutine compute_diffusion_coeff(T, n_tot, D)
    real(dp), dimension(:), intent(in) :: T, n_tot
    real(dp), dimension(:), intent(out) :: D
    integer :: i
    do i = 1, size(T)
       D(i) = (D0 * 1d17 * T(i)**s) / n_tot(i)
    end do
  end subroutine compute_diffusion_coeff

  !-----------------------------------------------------------------------
  !> Compute the effusion velocity (cm/s) for a given temperature and molecular mass.
  !-----------------------------------------------------------------------
  function effusion_velocity(T, molar_mass, z_top) result(v)
    real(dp), intent(in) :: T, molar_mass, z_top
    real(dp) :: lambda_H, vth, v
    lambda_H = (molar_mass * M_H * BIG_G * MARS_MASS) / (BOLTZMANN_K * T * 1.0d-2 * (MARS_RADIUS + z_top))
    vth = sqrt(2.0d0 * BOLTZMANN_K * T / (molar_mass * M_H))
    v = 1.0d2 * exp(-lambda_H) * vth * (lambda_H + 1.0d0) / (2.0d0 * sqrt(pi))
  end function effusion_velocity

  !-----------------------------------------------------------------------
  !> Initializes the hydrogen number density profile (can be anything, here I use 
  !> the case flux Փ [cm^-2 s^-1] = 0).
  !>
  !> Parameters:
  !>     z: Array of altitude values (cm).
  !>     n_bot: Hydrogen number density at the bottom boundary (cm^-3).
  !>     T: Temperature profile (K).
  !>     H: Scale height (cm).                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            
  !>
  !> Returns:
  !>     n_H: Hydrogen number density profile.
  !-----------------------------------------------------------------------
  subroutine initialize_H_steady(z, n_bot, T, H, n_H)
    real(dp), dimension(:), intent(in) :: z, T, H
    real(dp), intent(in) :: n_bot
    real(dp), dimension(:), intent(out) :: n_H
    integer :: i
    real(dp) :: dz, dTdz, dn_dz

    n_H(1) = n_bot
    do i = 2, size(z)
       dz = z(i) - z(i - 1)
       dTdz = (T(i) - T(i - 1)) / dz
       dn_dz = -n_H(i - 1) / H(i - 1) - n_H(i - 1) * ((1.0d0 + alpha_T) / T(i - 1)) * dTdz
       n_H(i) = n_H(i - 1) + dn_dz * dz
    end do
    ! do i = 2, size(z)
    !   n_H(i) = 0
    ! end do
  end subroutine initialize_H_steady

  !-----------------------------------------------------------------------
  !> Initializes the exponential hydrogen number density profile
  !>
  !> Parameters:
  !>     altitudes: Array of altitude values (cm).
  !>     n_bottom: Hydrogen number density at the bottom boundary (cm^-3).
  !>     H: Scale height (cm).                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            
  !>
  !> Returns:
  !>     n_H: Hydrogen number density profile (cm^-3).
  !-----------------------------------------------------------------------
  subroutine initialize_H_expo(altitudes, n_bottom, H, n_H)
    ! Inputs
    real(dp), dimension(:), intent(in) :: altitudes      ! in cm
    real(dp), intent(in) :: n_bottom                     ! in cm^-3
    real(dp), dimension(:), intent(in) :: H              ! in cm
    ! Output
    real(dp), dimension(size(altitudes)), intent(out) :: n_H          ! in cm^-3
    ! Locals
    integer :: i
    real(dp) :: alt0

    ! Compute the base altitude
    alt0 = altitudes(1)
    ! Compute the hydrogen density profile
    do i = 1, size(altitudes)
       n_H(i) = n_bottom * exp( - (altitudes(i) - alt0) / (1.0_dp * H(i)) )
    end do
  end subroutine initialize_H_expo

  !=======================================================================
  !  Subroutine: Thomas
  !  Solves a tridiagonal system Ax = d using Thomas algorithm
  !
  !  a:   (Input) Real(dp) array(N), the sub-diagonal. a(1) is ignored.
  !  b:   (Input) Real(dp) array(N), the main diagonal.
  !  c:   (Input) Real(dp) array(N), the super-diagonal. c(N) is ignored.
  !  rhs: (Input) Real(dp) array(N), the right-hand side vector.
  !  n_H_new:   (Output) Real(dp) array(N), the solution vector.
  !=======================================================================
  subroutine Thomas(a,b,c,rhs,n_H_new)
    real(dp), dimension(num_levels),intent(in) :: a, b, c, rhs
    real(dp), dimension(num_levels) :: C_prime, D_prime
    real(dp), dimension(:), intent(out) :: n_H_new
    real(dp) :: denom
    integer :: i, N

    N = size(a)
   ! Forward sweep
    C_prime(1) = c(1)/b(1)
    D_prime(1) = rhs(1)/b(1)
 
    do i = 2, N
     denom = b(i) - a(i)*C_prime(i-1)
     if (i<N) then
       C_prime(i) = c(i)/denom
     else
       C_prime(i) = 0
     end if
     D_prime(i) = (rhs(i) - a(i) * D_prime(i-1))/denom
    end do
   
    ! Backward substitution
   !  n_H_new = 0.0_dp ! Initialize to zero
    n_H_new(N) = D_prime(N)
    do i = N-1, 1, -1
     n_H_new(i) = D_prime(i) - C_prime(i) * n_H_new(i+1)
    end do
 
     ! Enforce non-negativity constraint on the solution
     do i = 1, N
       if (n_H_new(i) .lt. 0.0d0) then
         n_H_new(i) = 0.0d0
       end if
     end do

  end subroutine Thomas


  subroutine CN_matrix(&
    n_H     , D, dz, T, H,                     & 
    a,b,c,rhs)
    implicit none
    integer :: N, i
    real(dp), intent(in) :: dz
    real(dp), dimension(:), intent(in) :: n_H, D, T, H
    real(dp) :: D_i, D_ip1, D_im1, T_i, T_ip1, T_im1, H_i, H_ip1, H_im1, &
    dDdz, dTdz, dT2dz2, dHdz, phi_coeff
    real(dp), dimension(size(n_H)), intent(out) :: a, b, c, rhs

    N=size(D)

    do i = 2, N-1
      D_i   = D(i)
      D_ip1 = D(i+1)
      D_im1 = D(i-1)

      T_i   = T(i)
      T_ip1 = T(i+1)
      T_im1 = T(i-1)

      H_i   = H(i)
      H_ip1 = H(i+1)
      H_im1 = H(i-1)

      ! Centered difference approximations for derivatives
      dDdz   = (D_ip1 - D_im1) / (2.0d0 * dz)
      dTdz   = (T_ip1 - T_im1) / (2.0d0 * dz)
      dT2dz2 = (T_ip1 - 2.0d0*T_i + T_im1) / (dz**2)
      dHdz   = (H_ip1 - H_im1) / (2.0d0 * dz)

      ! φ = 1/H_i + ((1 + α_T)/T_i) * dTdz
      phi_coeff = (1.0d0 / H_i) + ((1.0d0 + alpha_T) / T_i) * dTdz

      a(i) = - D_i / dz**2 + ( dDdz + D_i*phi_coeff ) / &
                                            (2.0d0 * dz)
      b(i) = - ( -2.0d0*D_i / dz**2 + dDdz*phi_coeff + D_i*(&
              -1.0d0/(H_i**2)*dHdz + (1.0d0+alpha_T)*( &
                -1.0d0/(T_i**2)*dTdz**2 + 1.0d0/T_i*dT2dz2 ) ) )
      c(i) = - D_i / dz**2 - ( dDdz + D_i*phi_coeff ) / (2.0 * dz)

      rhs(i) = - ( a(i)*n_H(i-1) + b(i)*n_H(i) + c(i)*n_H(i+1) )

      ! if (ABS(b(i)) < (ABS(a(i)) + ABS(c(i)))) then
      !   print*, 'ERROR: Matrix is not diagonally dominant at i = ', i
      !   ! stop
      ! end if
    end do
  
  end subroutine CN_matrix


  subroutine interfaces(&
    n_H     , D, dz, T, H,                     & 
    a,b,c,rhs)
    implicit none
    integer :: N, i
    real(dp), intent(in) :: dz
    real(dp), dimension(:), intent(in) :: n_H, D, T, H
    real(dp) :: dTdz_ip1_2, dTdz_im1_2, phi_ip1_2, phi_im1_2
    real(dp) :: D_ip1_2, D_im1_2, T_ip1_2, T_im1_2, H_ip1_2, H_im1_2
    real(dp), dimension(size(n_H)), intent(out) :: a, b, c, rhs

    N=size(D)

    do i = 2, N-1
      ! --- Define properties at cell interfaces (i-1/2 and i+1/2) ---
      D_ip1_2 = (D(i + 1) + D(i)) / 2.0d0
      D_im1_2 = (D(i) + D(i - 1)) / 2.0d0

      T_ip1_2 = (T(i + 1) + T(i)) / 2.0d0
      T_im1_2 = (T(i) + T(i - 1)) / 2.0d0
      
      H_ip1_2 = (H(i + 1) + H(i)) / 2.0d0
      H_im1_2 = (H(i) + H(i - 1)) / 2.0d0

      ! --- Calculate the drift term coefficient 'phi' at the interfaces ---
      ! phi = 1/H + (1+alpha_T)/T * dT/dz
      dTdz_ip1_2 = (T(i + 1) - T(i)) / dz
      dTdz_im1_2 = (T(i) - T(i - 1)) / dz

      phi_ip1_2 = (1.0d0 / H_ip1_2) + ((1.0d0 + alpha_T) / T_ip1_2) * dTdz_ip1_2
      phi_im1_2 = (1.0d0 / H_im1_2) + ((1.0d0 + alpha_T) / T_im1_2) * dTdz_im1_2

      ! --- Calculate the tridiagonal coefficients a(i), b(i), c(i) ---
      ! Based on the flux balance equation: Flux(i+1/2) - Flux(i-1/2) = 0
      ! where Flux = -D * (dn/dz + n*phi)
      
            ! Coefficient for n(i-1)
            ! a(i) = D_im1_2 / dz - D_im1_2 * phi_im1_2 / 2.0d0
      a(i) = -D_im1_2 / dz + D_im1_2 * phi_im1_2
      
      ! Coefficient for n(i+1)
      ! c(i) = D_ip1_2 / dz + D_ip1_2 * phi_ip1_2 / 2.0d0
      c(i) = -D_ip1_2 / dz

      ! Coefficient for n(i)
!             b(i) = - (D_ip1_2+D_im1_2)/dz - (D_ip1_2*phi_ip1_2/2.0d0)
!      &                                    + (D_im1_2*phi_im1_2/2.0d0)
      b(i) = (D_ip1_2 / dz) + (D_im1_2 / dz) - D_ip1_2 * phi_ip1_2

      ! --- Set the Right-Hand Side ---
      ! For a steady-state problem, the RHS is zero for interior points.
      rhs(i) = 0.0d0

      ! if (ABS(b(i)) < (ABS(a(i)) + ABS(c(i)))) then
      !   print*, 'ERROR: Matrix is not diagonally dominant at i = ', i
      !   ! stop
      ! end if
    end do
  
  end subroutine interfaces

!--------------------------------------------------------------------------
!> @brief Constructs the tridiagonal matrix using a robust finite volume
!>        upwind scheme. This method is designed to be stable and avoid
!>        the unphysical oscillations that require scheme-switching.
!> @param D(:)      Diffusion coefficient profile [cm²/s].
!> @param dz        Grid spacing [cm].
!> @param T(:)      Temperature profile [K].
!> @param H(:)      Scale height profile [cm].
!> @param[out] a(:) The sub-diagonal.
!> @param[out] b(:) The main diagonal.
!> @param[out] c(:) The super-diagonal.
!--------------------------------------------------------------------------
  subroutine build_matrix_robust(D, dz, T, H, a, b, c,rhs)
        implicit none
        real(8), intent(in) :: dz
        real(8), dimension(:), intent(in)  :: D, T, H
        real(8), dimension(size(D)), intent(out) :: a, b, c,rhs
    
        integer :: i, N
        real(8) :: D_up, D_dn      ! Diffusion coeff at interfaces i+1/2 and i-1/2
        real(8) :: phi_up, phi_dn  ! Advection coeff 'phi' at interfaces
        real(8) :: T_up, T_dn, H_up, H_dn
        real(8) :: dTdz_up, dTdz_dn
        real(8) :: Pe_up, Pe_dn    ! Cell Peclet numbers (?) at interfaces
    
        N = size(D)
    
        ! Zero out the arrays
        a = 0.0d0
        b = 0.0d0
        c = 0.0d0
        rhs = 0.0d0
    
        do i = 2, N - 1
            ! --- Properties at the 'downstream' interface (i+1/2) ---
            D_up    = 0.5d0 * (D(i) + D(i+1))
            T_up    = 0.5d0 * (T(i) + T(i+1))
            H_up    = 0.5d0 * (H(i) + H(i+1))
            dTdz_up = (T(i+1) - T(i)) / dz
            phi_up  = (1.0d0/H_up)+((1.0d0+alpha_T)/T_up)*dTdz_up
            ! Drift term at i+1/2
            Pe_up   = D_up * phi_up
    
            ! --- Properties at the 'upstream' interface (i-1/2) ---
            D_dn    = 0.5d0 * (D(i) + D(i-1))
            T_dn    = 0.5d0 * (T(i) + T(i-1))
            H_dn    = 0.5d0 * (H(i) + H(i-1))
            dTdz_dn = (T(i) - T(i-1)) / dz
            phi_dn  = (1.0d0/H_dn)+((1.0d0+alpha_T)/T_dn)*dTdz_dn
            ! Drift term at i-1/2
            Pe_dn   = D_dn * phi_dn
    
            ! --- Calculate tridiagonal coefficients a(i), b(i), c(i) ---
            ! Based on the flux balance: F_up - F_dn = 0
            ! F = Diffusive_Flux + Advective_Flux
            ! Diffusive_Flux = -D/dz * (n_i+1 - n_i)
            ! Advective_Flux (upwinded) = -Pe * n_i   (if Pe > 0)
            !                             -Pe * n_i+1 (if Pe < 0)
    
            ! Coefficient for n(i-1)
            ! a(i) = -D_dn / dz**2 - max(Pe_dn, 0.0d0) / dz
            a(i) = -D_dn / dz**2 + min(Pe_dn, 0.0d0) / dz
    
            ! Coefficient for n(i+1)
            ! c(i) = -D_up / dz**2 + min(Pe_up, 0.0d0) / dz
            c(i) = -D_up / dz**2 - max(Pe_up, 0.0d0) / dz
    
            ! Coefficient for n(i)
            ! b(i) = (D_up + D_dn) / dz**2 + max(-Pe_up, 0.0d0)/dz - min(Pe_dn, 0.0d0) / dz
            b(i) = (D_up + D_dn) / dz**2 - min(Pe_up, 0.0d0)/dz + max(Pe_dn, 0.0d0) / dz
    
            ! The right-hand-side for the steady state is 0 for interior points.
            ! if (ABS(b(i)) < (ABS(a(i)) + ABS(c(i)))) then
            !   print*, 'ERROR: Matrix is not diagonally dominant at i = ', i
            !   ! stop
            ! end if
        end do
    
    end subroutine build_matrix_robust


      !--------------------------------------------------------------------------
      !> @brief  Assemble the tridiagonal coefficient matrices (a, b, c) and
      !>         RHS vector for the steady-state 1D advection–diffusion equation.
      !>         This is a finite-difference steady-state formulation.
      !>         Diffusion term: Centered difference (second order).
      !>         Advection term: Upwind scheme (first order).
      !> @param D(:)      Diffusion coefficient profile [cm²/s].
      !> @param dz        Grid spacing [cm].
      !> @param T(:)      Temperature profile [K].
      !> @param H(:)      Scale height profile [cm].
      !> @param[out] a(:) The sub-diagonal.
      !> @param[out] b(:) The main diagonal.
      !> @param[out] c(:) The super-diagonal.
      !--------------------------------------------------------------------------
      subroutine finite_difference(D, dz, T, H, a, b, c, rhs)
            implicit none
            real(8), intent(in) :: dz
            real(8), dimension(:), intent(in)  :: D, T, H
            real(8), dimension(size(D)), intent(out) :: a, b, c, rhs
        
            integer :: i, N
            real(8) :: D_up, D_dn      ! Diffusion coeff at interfaces i+1/2 and i-1/2
            real(8) :: phi_up, phi_dn  ! Advection coeff 'phi' at interfaces
            real(8) :: T_up, T_dn, H_up, H_dn
            real(8) :: dTdz_up, dTdz_dn
            real(8) :: v_up, v_dn    ! drift terms at interfaces (cm/s)
        
            N = size(D)
        
            ! Zero out the arrays
            a = 0.0d0
            b = 0.0d0
            c = 0.0d0
            rhs = 0.0d0
        
            do i = 2, N - 1
                ! --- Properties at the 'downstream' interface (i+1/2) ---
                D_up    = 0.5d0 * (D(i) + D(i+1))
                T_up    = 0.5d0 * (T(i) + T(i+1))
                H_up    = 0.5d0 * (H(i) + H(i+1))
                dTdz_up = (T(i+1) - T(i)) / dz
                phi_up  = (1.0d0/H_up)+((1.0d0+alpha_T)/T_up)*dTdz_up
                ! Drift term at i+1/2
                v_up   = D_up * phi_up
        
                ! --- Properties at the 'upstream' interface (i-1/2) ---
                D_dn    = 0.5d0 * (D(i) + D(i-1))
                T_dn    = 0.5d0 * (T(i) + T(i-1))
                H_dn    = 0.5d0 * (H(i) + H(i-1))
                dTdz_dn = (T(i) - T(i-1)) / dz
                phi_dn  = (1.0d0/H_dn)+((1.0d0+alpha_T)/T_dn)*dTdz_dn
                ! Drift term at i-1/2
                v_dn   = D_dn * phi_dn
        
                ! --- Calculate tridiagonal coefficients a(i), b(i), c(i) ---
                ! Based on the flux balance: F_up - F_dn = 0
                ! F = Diffusive_Flux + Advective_Flux
                ! Diffusive_Flux = -D/dz * (n_i+1 - n_i)
                ! Advective_Flux (upwinded) = -Pe * n_i   (if Pe > 0)
                !                             -Pe * n_i+1 (if Pe < 0)
        
                ! Coefficient for n(i-1)
            !     a(i) = -D_dn / dz**2 - max(v_dn, 0.0d0) / dz
                a(i) = -D_dn / dz**2 + min(v_dn, 0.0d0) / dz
        
                ! Coefficient for n(i+1)
            !     c(i) = -D_up / dz**2 + min(v_up, 0.0d0) / dz
                c(i) = -D_up / dz**2 - max(v_up, 0.0d0) / dz
        
                ! Coefficient for n(i)
            !     b(i) = (D_up + D_dn) / dz**2 + max(-v_up, 0.0d0)/dz
!      &                 - min(v_dn, 0.0d0) / dz
                b(i) = (D_up + D_dn) / dz**2 - min(v_up, 0.0d0)/dz + max(v_dn, 0.0d0) / dz
        
                ! The right-hand-side for the steady state is 0 for interior points.
            end do
        
        end subroutine finite_difference

  !-----------------------------------------------------------------------
  !> Solves continuity equation dn/dt = -dՓ/dz using Crank-Nicolson scheme
  !> and use the Thomas algorithm for solving the tridiagonal system.
  !> Parameters:
  !>     n_H: Input number density profile (cm^-3).
  !>     D: Molecular diffusion coefficient D(z) (cm^2/s).
  !>     dz: Grid spacing (cm). 
  !>     T: Temperature profile (K).                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         
  !>     H: Scale height (cm).
  !>     v_eff: Effusion velocity (cm/s).
  !>
  !> Returns:
  !>     n_H_new: New H number density profile (cm^-3).
  !-----------------------------------------------------------------------
  ! subroutine crank_nicolson(n_H, D, dz, T, H, v_eff, n_H_new)
  subroutine crank_nicolson ( &
        n_H     , D, dz, T, H,                     & 
        num_levels, bc_top_type , bc_top_val ,              &
        alpha_T , z_top,                              &  ! species‑specific thermal‐diff coeff
        n_H_new)

   implicit none
   integer, intent(in) :: num_levels
   real(dp), intent(in) :: dz, z_top, alpha_T, bc_top_val !, v_eff
   real(dp), dimension(:), intent(in) :: n_H, D, T, H
   real(dp), dimension(:), intent(out) :: n_H_new
   character(len=*), intent(in) :: bc_top_type
   real(dp) :: v_eff
 
   integer :: N, i
   real(dp), dimension(size(n_H)) :: a, b, c, rhs
   real(dp), dimension(size(n_H)) :: C_prime, D_prime
   real(dp) :: D_i, D_ip1, D_im1, T_i, T_ip1, T_im1, H_i, H_ip1, H_im1, &
              dDdz, dTdz, dT2dz2, dHdz, phi_coeff, D_top, H_top, T_top, &
              coeff, beta, gamma, denom
   real(dp) :: dTdz_ip1_2, dTdz_im1_2, phi_ip1_2, phi_im1_2
   real(dp) :: D_ip1_2, D_im1_2, T_ip1_2, T_im1_2, H_ip1_2, H_im1_2

   N = size(D)

   ! Set up tridiagonal matrix
  !  call CN_matrix(n_H, D, dz, T, H, a, b, c, rhs)
  !  call build_matrix_robust(D, dz, T, H, a, b, c, rhs)
  !  call interfaces(n_H, D, dz, T, H, a, b, c, rhs)
    call finite_difference(D, dz, T, H, a, b, c,rhs)

   ! Bottom boundary conditions
   b(1) = 1.0d0
   c(1) = 0.0d0
   a(1) = 0.0d0
   rhs(1) = n_H(1)

   ! Top Robin boundary
   D_top = (D(N) + D(N-1)) / 2.0d0
   H_top = (H(N) + H(N-1)) / 2.0d0
   T_top = (T(N) + T(N-1)) / 2.0d0
   dTdz = (T(N) - T(N-1)) / dz

   coeff = 1.0d0 / H_top + (1.0d0 + alpha_T) / T_top * dTdz
  !  beta = v_eff * dz + D_top + D_top * coeff * dz
  !  gamma = -D_top

  !  a(N) = gamma
  !  b(N) = beta
  !  c(N) = 0.0d0
  !  rhs(N) = 0.0d0

   ! ─── top boundary ( i = N ) ─────────────────────────────────────────────
    select case (bc_top_type)

    case ('n')          ! n = constant
      print*, 'select case n'
      a(N)   = -1.0d0
      b(N)   = 1.0d0
      c(N)   = 0.0d0
      rhs(N) = 0.0d0

    case ('f')               ! F = constant  (positive = upward escape)
      !  +D ∂n/∂z = F   (sign flipped because upward is +z)
      print*, 'select case f'
      a(N)   = -D_top
      b(N)   = D_top + D_top * coeff * dz
      c(N)   = 0.0d0
      rhs(N) = bc_top_val * dz

    case ('v')           ! v_escape · n  (Robin condition)
      print*, 'select case v'
      !  −D ∂n/∂z = v n   ➋
      !  FD:  (n(N) − n(N-1)) / dz = -v/D * n(N)
      v_eff = effusion_velocity(T(num_levels), bc_top_val, z_top)
      print*, 'v_eff = ', v_eff
      beta = v_eff * dz + D_top + D_top * coeff * dz
      gamma = -D_top

      ! a(N)    = -D_top/dz - D_top*coeff
      ! b(N)    = D_top/dz - v_eff

      a(N) = gamma
      b(N) = beta
      c(N) = 0.0d0
      rhs(N) = 0.0d0

    case default
      stop "Unknown top BC"

    end select


   !-------------------------------
   ! Thomas algorithm
   !-------------------------------

    call Thomas(a, b, c, rhs, n_H_new)
    ! print*, 'n_H bottom = ', n_H_new(1), &
          ! ' n_H top = ', n_H_new(num_levels)

    if (n_H_new(1) < n_H_new(2)) then
      print*,'CN failed: n_H bottom = ', n_H_new(1), &
            ' n_H next = ', n_H_new(2), &
            ' n_H top = ', n_H_new(num_levels)
      print*, 'choosing interfaces'
      call interfaces(n_H, D, dz, T, H, a, b, c, rhs)
      ! Reapply the boundary conditions
      ! Bottom boundary conditions
      b(1) = 1.0d0
      c(1) = 0.0d0
      a(1) = 0.0d0
      rhs(1) = n_H(1)

      ! Top Robin boundary
      D_top = (D(N) + D(N-1)) / 2.0d0
      H_top = (H(N) + H(N-1)) / 2.0d0
      T_top = (T(N) + T(N-1)) / 2.0d0
      dTdz = (T(N) - T(N-1)) / dz

      coeff = 1.0d0 / H_top + (1.0d0 + alpha_T) / T_top * dTdz
      !  beta = v_eff * dz + D_top + D_top * coeff * dz
      !  gamma = -D_top

      !  a(N) = gamma
      !  b(N) = beta
      !  c(N) = 0.0d0
      !  rhs(N) = 0.0d0

      ! ─── top boundary ( i = N ) ─────────────────────────────────────────────
        select case (bc_top_type)

        case ('n')          ! n = constant
          print*, 'select case n'
          a(N)   = -1.0d0
          b(N)   = 1.0d0
          c(N)   = 0.0d0
          rhs(N) = 0.0d0

        case ('f')               ! F = constant  (positive = upward escape)
          !  +D ∂n/∂z = F   (sign flipped because upward is +z)
          print*, 'select case f'
          a(N)   = -D_top
          b(N)   = D_top + D_top * coeff * dz
          c(N)   = 0.0d0
          rhs(N) = bc_top_val * dz

        case ('v')           ! v_escape · n  (Robin condition)
          print*, 'select case v'
          !  −D ∂n/∂z = v n   ➋
          !  FD:  (n(N) − n(N-1)) / dz = -v/D * n(N)
          v_eff = effusion_velocity(T(num_levels), bc_top_val, z_top)
          print*, 'v_eff = ', v_eff
          beta = v_eff * dz + D_top + D_top * coeff * dz
          gamma = -D_top
      
          a(N) = gamma
          b(N) = beta
          c(N) = 0.0d0
          rhs(N) = 0.0d0

        case default
          stop "Unknown top BC"

        end select
        call Thomas(a, b, c, rhs, n_H_new)
    end if

  end subroutine crank_nicolson

  !-----------------------------------------------------------------------
  ! Running the CN scheme with inputs from the GCM.
  !> Parameters:
  !>     num_levels: Number of altitude levels for grid spacing.
  !>     z_bottom: Bottom altitude (cm).
  !>     z_top: Top altitude (cm).
  !>     n_bottom: Input number density of H from top of GCM (cm^-3).
  !>     n_CO2_bottom: Input number density of CO2 from top of GCM (cm^-3).
  !>     tol: Tolerance for error.
  !>     max_iter: Max number of iteration for convergence.
  !>     Tinf: Temperature at exobase (K), typically 240K?
  !>
  !> Returns:
  !>     altitudes: Array of altitude values (cm).
  !>     n_H: New H number density profile (cm^-3).
  !-----------------------------------------------------------------------
  ! subroutine diffusion(num_levels,z_bottom,z_top,n_bottom,n_CO2_bottom,Tinf, &
  !         altitudes, n_H, n_H_ini, H_effusion_velocity, total_densities, T_profile)
  subroutine diffusion(num_levels,z_bottom,z_top,n_bottom,n_CO2_bottom,params, &
                      bc_top_type, bc_top_val, &
                      altitudes, n_H, n_H_ini, total_densities, T_profile, &
                      H_effusion_velocity, H_profile,D &
                      )
    implicit none
    ! Input arguments
    integer, intent(in) :: num_levels
    real(dp), intent(in) :: z_bottom, z_top
    real(dp), intent(in) :: n_bottom
    real(dp), intent(in) :: n_CO2_bottom
    ! real(dp), intent(in) :: Tinf
    real(dp), intent(in) :: params(4) ! Custom temperature parameters
    real(dp), intent(in) :: bc_top_val
    character(len=*), intent(in) :: bc_top_type
    ! Output
    real(dp), dimension(num_levels), intent(out) :: altitudes, n_H, n_H_ini, total_densities, T_profile
    real(dp), dimension(num_levels), intent(out) :: H_profile, D
    real(dp), intent(out) :: H_effusion_velocity
    ! Local variables
    real(dp) :: dz
    ! real(dp), dimension(num_levels) ::  H_profile,D
    real(dp), dimension(num_levels) :: n_H_new
    ! real(dp) :: H_effusion_velocity
    integer :: i
  
    ! Initialize grid
    call initialize_grid(z_bottom, z_top, num_levels, altitudes, dz)

    altitudes = altitudes*1.0d-5 ! cm to km
  
    ! Temperature and scale height
    do i = 1, num_levels
        ! T_profile(i) = temperature_profile(altitudes(i), Tinf)
        T_profile(i) = Tcustom(altitudes(i),params)
    end do

    altitudes = altitudes*1.0d5 ! km to cm

    ! If n_bottom is non-physical or zero, no H is present.
    if (n_bottom .gt. 0.0d0) then
        ! --- NORMAL CALCULATION ---
        ! Calculate physical profiles based on the new T(z) grid
        call scale_height_profile(altitudes, T_profile, H_profile)
      
        ! Initialize hydrogen profile
        call initialize_H_steady(altitudes, n_bottom, T_profile, H_profile, n_H_ini)
        ! call initialize_H_expo(altitudes, n_bottom, H_profile, n_H_ini)
      
        ! Total densities
        call total_density_CO2(altitudes,dz,n_CO2_bottom,T_profile,total_densities)
      
        ! Diffusion coefficients
        call compute_diffusion_coeff(T_profile, total_densities, D)
      
        ! Effusion velocity
        H_effusion_velocity = effusion_velocity(T_profile(num_levels), 1.0d0, z_top)
      
        n_H = n_H_ini

        ! call crank_nicolson(n_H, D, dz, T_profile, H_profile, H_effusion_velocity, n_H_new)
        call crank_nicolson(n_H, D, dz, T_profile, H_profile, num_levels,bc_top_type, bc_top_val, alpha_T, z_top, n_H_new)
          
        n_H = n_H_new
      else
        ! --- NO HYDROGEN AT THE BOTTOM: SET ALL H-RELATED OUTPUTS TO ZERO ---
        H_effusion_velocity = 0.0d0
        n_H(:) = 0.0d0
        n_H_ini(:) = 0.0d0
        ! Still calculate the background atmosphere properties for diagnostics
        call scale_height_profile(altitudes, T_profile, H_profile)
        call total_density_CO2(altitudes, dz, n_CO2_bottom,T_profile, total_densities)
        call compute_diffusion_coeff(T_profile, total_densities,D)
    end if


  end subroutine diffusion

end module diffusion_module

