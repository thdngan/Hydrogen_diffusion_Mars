module diffusion_module
  implicit none
  ! Define a double precision kind constant to reuse, for consistency in precision.
  integer, parameter :: dp = kind(1.0d0)  
  ! Physical constants
  real(dp), parameter :: BOLTZMANN_K = 1.38d-23            ! J/K
  real(dp), parameter :: BIG_G = 6.67d-11                  ! N m^2/kg^2  (gravitational constant)
  real(dp), parameter :: M_H = 1.67d-27                    ! kg (mass of a proton)
  real(dp), parameter :: R_gas = 8.314462618d0             ! N m K^-1 mol^-1  (gas constant)
  real(dp), parameter :: M_CO2 = 0.04401d0                 ! kg/mol (molar mass of CO2)
  real(dp), parameter :: alpha_T = -0.25d0                 ! Thermal diffusion factor of H (from Krasnopolsky 2002)
  ! real(dp), parameter :: pi = 4.0d0*DATAN(10.d0)       
  real(dp), parameter :: pi = 3.14159265358979d0 
  ! Mars parameters
  real(dp), parameter :: MARS_MASS = 0.1075d0 * 5.972d24
  real(dp), parameter :: MARS_RADIUS = 3396.0d5    

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
    real(dp), parameter :: Ttropo = 125.0d0
    real(dp), parameter :: ztropo = 90.0d5
    real(dp), parameter :: zwidth = 30.0d5
    real(dp), parameter :: lapse = -1.4d-5

    if (z >= ztropo) then
       T = Tinf - (Tinf - Ttropo) * exp(-((z - ztropo)**2) / (8.d10 * Tinf))
    else if ((z >= ztropo - zwidth) .and. (ztropo > z)) then
       T = Ttropo
    else if (ztropo - zwidth > z) then
       T = Ttropo - lapse * (ztropo - zwidth - z)
    end if
  end function temperature_profile

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
    do i = 1, size(z)
       H(i) = (BOLTZMANN_K * T(i) / (1.0d0 * M_H * BIG_G * MARS_MASS) * ((z(i) + MARS_RADIUS)*1.0d-2)**2) * 1.0d2
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
    real(dp), parameter :: D0 = 8.4d0, s = 0.597d0
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
      ! n_H(i) = 0.0d0
    end do
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
  subroutine crank_nicolson(n_H, D, dt, dz, T, H, v_eff, n_bottom, n_H_new)
   implicit none
   real(dp), intent(in) :: dz, v_eff, n_bottom, dt
   real(dp), dimension(:), intent(in) :: n_H, D, T, H
   real(dp), dimension(:), intent(out) :: n_H_new
 
   integer :: N, i
   real(dp), dimension(size(n_H)) :: a, b, c, rhs
   real(dp), dimension(size(n_H)) :: C_prime, D_prime
   real(dp) :: D_i, D_ip1, D_im1, T_i, T_ip1, T_im1, H_i, H_ip1, H_im1, &
              dDdz, dTdz, dT2dz2, dHdz, phi_coeff, D_top, H_top, T_top, &
              coeff, beta, gamma, denom

   N = size(D)

   ! Set up tridiagonal matrix
   do i = 2, N-1
      D_i = D(i)
      D_ip1 = D(i+1)
      D_im1 = D(i-1)

      T_i = T(i)
      T_ip1 = T(i+1)
      T_im1 = T(i-1)

      H_i = H(i)
      H_ip1 = H(i+1)
      H_im1 = H(i-1)

      dDdz = (D_ip1 - D_im1) / (2.0d0 * dz)
      dTdz = (T_ip1 - T_im1) / (2.0d0 * dz)
      dT2dz2 = (T_ip1 - 2.0d0 * T_i + T_im1) / dz**2.0d0
      dHdz = (H_ip1 - H_im1) / (2.0d0 * dz)

      phi_coeff = (1.0d0 / H_i) + ((1.0d0 + alpha_T) / T_i) * dTdz

      a(i) = (dt / 2.0d0) * (-D_i / dz**2.0d0 + (dDdz + D_i * phi_coeff) / (2.0d0 * dz))
      b(i) = 1 - (dt / 2.0d0) * (-2.0d0 * D_i / dz**2.0d0 + dDdz * phi_coeff + D_i * (&
            -1.0d0 / H_i**2.0d0 * dHdz + (1.0d0 + alpha_T) * (&
            -1.0d0 / T_i**2.0d0 * dTdz**2.0d0 + 1.0d0 / T_i * dT2dz2)))
      c(i) = (dt / 2.0d0) * (-D_i / dz**2.0d0 - (dDdz + D_i * phi_coeff) / (2.0d0 * dz))

      rhs(i) = -(n_H(i-1)*a(i) + n_H(i)*(b(i)-2.0d0) + n_H(i+1)*c(i))
   end do

   ! Bottom boundary conditions
   b(1) = 1.0d0
   c(1) = 0.0d0
   a(1) = 0.0d0
   rhs(1) = n_bottom

   ! Top boundary
   D_top = (D(N) + D(N-1)) / 2.0d0
   H_top = (H(N) + H(N-1)) / 2.0d0
   T_top = (T(N) + T(N-1)) / 2.0d0
   dTdz = (T(N) - T(N-1)) / dz

   coeff = 1.0d0 / H_top + (1.0d0 + alpha_T) / T_top * dTdz

   beta = v_eff * dz + D_top + D_top * coeff * dz
   gamma = -D_top

   a(N) = gamma
   b(N) = beta
   c(N) = 0.0d0
   rhs(N) = 0.0d0

   !-------------------------------
   ! Thomas algorithm
   !-------------------------------
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

  end subroutine crank_nicolson

  !-----------------------------------------------------------------------
  ! Running the CN scheme with inputs from the GCM.
  !> Parameters:
  !>     num_levels: Number of altitude levels for grid spacing.
  !>     z_bottom: Bottom altitude (cm).
  !>     z_top: Top altitude (cm).
  !>     n_bottom: Input number density of H from top of GCM (cm^-3).
  !>     n_CO2_bottom: Input number density of CO2 from top of GCM (cm^-3).
  !>     Tinf: Temperature at exobase (K), typically 240K?
  !>     dt: Time step (seconds)
  !>     n_H_prev: Number density of H from previous timestep (if
  !>               it is the very first timestep, we have to generate
  !>               an initial profile first and insert it here)
  !>
  !> Returns:
  !>     altitudes: Array of altitude values (cm).
  !>     n_H: New H number density profile (cm^-3).
  !>     H_effusion_velocity: Effusion velocity (cm/s) at the exobase
  !-----------------------------------------------------------------------
  subroutine diffusion(num_levels,z_bottom,z_top,n_bottom,n_CO2_bottom,Tinf,dt,n_H_prev, &
                      altitudes, n_H_new, H_effusion_velocity)
    implicit none
    ! Input arguments
    integer, intent(in) :: num_levels
    real(dp), intent(in) :: z_bottom, z_top
    real(dp), intent(in) :: n_bottom
    real(dp), intent(in) :: n_CO2_bottom
    real(dp), intent(in) :: Tinf
    real(dp), intent(in) :: dt
    real(dp), dimension(num_levels), intent(in) :: n_H_prev
    ! Output
    real(dp), dimension(num_levels), intent(out) :: altitudes, n_H_new
    real(dp), intent(out) :: H_effusion_velocity
    ! Local variables
    real(dp) :: dz
    real(dp), dimension(num_levels) ::  T_profile, H_profile, D, total_densities
    ! real(dp) :: H_effusion_velocity
    integer :: i
  
    ! Initialize grid
    call initialize_grid(z_bottom, z_top, num_levels, altitudes, dz)
  
    ! Temperature and scale height
    do i = 1, num_levels
        T_profile(i) = temperature_profile(altitudes(i), Tinf)
    end do
  
    call scale_height_profile(altitudes, T_profile, H_profile)
  
    ! Total densities
    call total_density_CO2(altitudes,dz,n_CO2_bottom,T_profile,total_densities)
  
    ! Diffusion coefficients
    call compute_diffusion_coeff(T_profile, total_densities, D)
  
    ! Effusion velocity
    H_effusion_velocity = effusion_velocity(T_profile(num_levels), 1.0d0, z_top)


    call crank_nicolson(n_H_prev, D,dt, dz, T_profile, H_profile, H_effusion_velocity, n_bottom, n_H_new)

  end subroutine diffusion

  ! subroutine diffusion_i(num_levels,z_bottom,z_top,n_bottom,n_CO2_bottom,Tinf,dt, &
  !                       altitudes, n_H_new, H_effusion_velocity)
  !   implicit none
  !   ! Input arguments
  !   integer, intent(in) :: num_levels
  !   real(dp), intent(in) :: z_bottom, z_top
  !   real(dp), intent(in) :: n_bottom
  !   real(dp), intent(in) :: n_CO2_bottom
  !   real(dp), intent(in) :: Tinf
  !   real(dp), intent(in) :: dt
  !   ! Output
  !   real(dp), dimension(num_levels), intent(out) :: altitudes, n_H_new
  !   real(dp), intent(out) :: H_effusion_velocity
  !   ! Local variables
  !   real(dp) :: dz
  !   real(dp), dimension(num_levels) ::  T_profile, H_profile, D, n_H_ini, total_densities
  !   ! real(dp) :: H_effusion_velocity
  !   integer :: i

  !   ! Initialize grid
  !   call initialize_grid(z_bottom, z_top, num_levels, altitudes, dz)

  !   ! Temperature and scale height
  !   do i = 1, num_levels
  !   T_profile(i) = temperature_profile(altitudes(i), Tinf)
  !   end do

  !   call scale_height_profile(altitudes, T_profile, H_profile)

  !   call initialize_H_steady(altitudes, n_bottom, T_profile, H_profile, n_H_ini)

  !   ! Total densities
  !   call total_density_CO2(altitudes,dz,n_CO2_bottom,T_profile,total_densities)

  !   ! Diffusion coefficients
  !   call compute_diffusion_coeff(T_profile, total_densities, D)

  !   ! Effusion velocity
  !   H_effusion_velocity = effusion_velocity(T_profile(num_levels), 1.0d0, z_top)


  !   call crank_nicolson(n_H_ini, D,dt, dz, T_profile, H_profile, H_effusion_velocity, n_bottom, n_H_new)

  ! end subroutine diffusion_i

end module diffusion_module

