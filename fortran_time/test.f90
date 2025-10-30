program test
  use diffusion_module
  implicit none
  integer, parameter :: num_levels = 100                       ! Number of altitude levels for grid spacing
  real(dp), parameter :: z_bottom = 150.0d5, z_top = 220.0d5   ! Bottom & top altitudes (cm)
  real(dp), parameter :: n_bottom_ini = 332645.34485837d0      ! Input number density of H from top of GCM (cm^-3)
  real(dp), parameter :: n_CO2_bottom_ini = 1.0d10             ! Input number density of CO2 from top of GCM (cm^-3)
  real(dp), parameter :: Tinf = 240.0d0                        ! Temperature at exobase (K), typically 240K?
  real(dp), parameter :: dt = 10.0d0                           ! Timestep of GCM, 10 seconds?

  real(dp), dimension(num_levels) :: altitudes, n_H, n_H_prev, n_H_ini
  real(dp), dimension(num_levels) :: T_profile, H_profile ,rel_errors
  real(dp) :: dz, v_eff
  integer :: i, j

  ! Below is an example time loop, about one Martian year here

  !################################################################################
  !##################### PARAMETERS FOR TESTING TIME LOOP #########################
  !################################################################################

  ! All of the parameters are just for this test program, we're not gonna use them in the model
  integer, parameter :: num_steps = 5935490
  real(dp), allocatable :: n_bottom_list(:), n_CO2_bottom_list(:)
  real(dp), parameter :: n_bottom_step = -1.0d-2
  real(dp) :: n_bottom, n_CO2_bottom, d
  real(dp), parameter :: tol = 1.0d-3
  
  ! These parameters are for printing the output, depends on whether we want to or not, might be used in the model
  integer, parameter :: output_interval = 100  ! Write output every 10 steps (100 000 for 11 days)
  character(len=100) :: filename  ! Variable to hold the dynamic filename
  character(len=100) :: folder ! Create output folder name
  folder = "output_0-3sol/"

  ! Create output directory (system-dependent: works on Unix/Linux/Mac)
  ! call system("mkdir -p output_0-3sol")
  ! call system("if not exist output mkdir output_100")

  !################################################################################
  !######################## VARYING DENSITY AT BOTTOM #############################
  !################################################################################
  ! Changing nH and nCO2 at 150km: here I'm just making them decrease gradually, it would
  ! definitely be different in the GCM though, like due to seasonal variations, etc.
  ! So basically this step is not needed when we actually insert the module into the GCM.
  allocate(n_bottom_list(num_steps), n_CO2_bottom_list(num_steps))
  do i = 1, num_steps
    n_bottom_list(i) = n_bottom_ini !+ (i-1)*n_bottom_step
  end do

  do i = 1, num_steps
    n_CO2_bottom_list(i) = n_CO2_bottom_ini !+ (i-1)*n_bottom_step
  end do

  !################################################################################
  !############################## INITIAL TIME STEP ###############################
  !################################################################################
  ! For the first time step, we have to initialize a n_H profile from top of GCM to exobase
  ! because we have no prior information, then the later timesteps can just use the computed
  ! n_H profile from previous timestep as initial condition

  ! Initialize grid
  call initialize_grid(z_bottom, z_top, num_levels, altitudes, dz)
  ! Temperature and scale height
  do i = 1, num_levels
      T_profile(i) = temperature_profile(altitudes(i), Tinf)
  end do
  call scale_height_profile(altitudes, T_profile, H_profile)
  call initialize_H_steady(altitudes, n_bottom_ini, T_profile, H_profile, n_H_ini)
  ! call initialize_H_expo(altitudes, n_bottom_ini, H_profile, n_H_ini)

  print*, 'Finished initializing'

  n_H_prev = n_H_ini

  !################################################################################
  !######################### TIME LOOP AND PRINT OUTPUT ###########################
  !################################################################################
  ! Put the content inside this loop into the GCM, because the GCM is a time loop itself (right?)
  do i = 1, num_steps
    ! print*, 'Doing step ', num_steps
    n_bottom = n_bottom_list(i)               ! In the GCM, just use n_H from the top altitude, basically n_bottom_ini
    n_CO2_bottom = n_CO2_bottom_list(i)       ! Same for this, use n_CO2 from the top altitude, basically n_CO2_bottom_ini
    
    call diffusion(num_levels,z_bottom,z_top,n_bottom,n_CO2_bottom,Tinf,dt,n_H_prev, &
                    altitudes, n_H, v_eff)
    rel_errors = NORM2(n_H - n_H_prev)/NORM2(n_H)

    d = n_H_prev(num_levels) - n_H(num_levels)
    n_H_prev = n_H
    ! if ((i<=10000) .and. (mod(i, output_interval)== 0)) then
    !   print *,d
    ! end if
    
    ! Only write output every output_interval steps
    ! if (mod(i, output_interval) == 0) then
    ! .and. (d>=tol))
    if ((i<=10000) .and. (mod(i, output_interval)== 0)) then 
      ! Construct filename inside 'output/' folder
      write(filename, '(A,I0,A)') trim(folder)//"nH_output_", i, ".dat"

      open(unit=10, file=filename, status='replace')
      write(10, *) v_eff
      do j = 1, num_levels
        write(10,*) altitudes(j)*1.0d-5, n_H(j), rel_errors(j)
      end do
      close(10)
      print *, 'Output written to ', trim(filename)
    end if
  end do

end program test