program test
  use diffusion_module
  implicit none
  ! real(dp), parameter :: z_bottom = 139.485d5, z_top = 239.485d5   ! Bottom & top altitudes (cm)
  ! real(dp), parameter :: n_CO2_bottom = 1.7961d9    ! Input number density of CO2 from top of GCM (cm^-3)
  ! real(dp), dimension(4), parameter :: params = [134.802d0, 104.556d0,188.695d0,2281.62d0]
  ! real(dp), parameter :: H_bottom = 12.3327d0          ! Input number density of H from top of GCM (cm^-3)

  real(dp), parameter :: z_bottom = 150.0d5, z_top = 220.0d5   ! Bottom & top altitudes (cm)
  real(dp), parameter :: n_CO2_bottom = 10000000000.0d0 ! 2381112947.409259d0    ! Input number density of CO2 from top of GCM (cm^-3)
  real(dp), dimension(4), parameter :: params = [125.0d0, 90.0d0,240.0d0,8.0d0]
  real(dp), parameter :: H_bottom = 332645.34485837002d0 !332645.34485836583d0          ! Input number density of H from top of GCM (cm^-3)

  real(dp), parameter :: H2_bottom = 962025.8509056892d0        ! H2
  real(dp), parameter :: CO2_bottom = 2059572654.8110359d0        ! CO2
  real(dp), parameter :: N2_bottom = 128911578.8371078d0       ! N2
  real(dp), parameter :: O_bottom = 84962944.0424445d0        ! O
  real(dp), parameter :: O2_bottom = 7983453.272630622d0        ! O2

  character(len=100), parameter :: output_file = "nH_test.dat"     ! Output file name

  real(dp), dimension(num_levels) :: altitudes, n_H, n_H_ini,total_densities, T_profile
  real(dp), dimension(num_levels) :: H_profile, D
  real(dp) :: v_eff
  integer :: i

  ! const speciesbclist=Dict(
  !              :CO2=>["n" 2.1e17; "f" 0.],
  !              :N2=>["n" 1.9e-2*2.1e17; "f" 0.],
  !              :O=>["f" 0.; "f" 1.2e8],
  !              :O2=>["f" 0.; "f" 0.],
  !              :H2=>["f" 0.; "v" H2_effusion_velocity],
  !              :H=>["f" 0.; "v" H_effusion_velocity],
  !              );

  ! call diffusion(num_levels,z_bottom,z_top,n_bottom,n_CO2_bottom,Tinf, &
                  ! altitudes, n_H, n_H_ini, v_eff, total_densities, T_profile)

  ! other species
  call diffusion(num_levels,z_bottom,z_top,H_bottom,n_CO2_bottom,params, &
                  'v', 1.0d0, &
                  altitudes, n_H, n_H_ini, total_densities, T_profile,&
                  v_eff, H_profile, D)

  


  ! Output result to file
  open(unit=10, file=output_file)
  ! write(10, *) v_eff
  ! write(10, *) '  altitude (km),            n_H (cm^-3),              T (K),                    H_profile (cm),'&
              !  ,'           n_H_ini (cm^-3),          total_densities (cm^-3),  D (cm^2/s)'
  do i = 1, num_levels
    write(10,*) altitudes(i)*1.0d-5, n_H(i), n_H_ini(i) !, T_profile(i), H_profile(i), n_H_ini(i), total_densities(i), D(i)
    ! write(10,*) altitudes(i), total_densities(i)
  end do
  close(10)
  print*, n_H(num_levels)
  print *, 'Output written to ',output_file

end program test

