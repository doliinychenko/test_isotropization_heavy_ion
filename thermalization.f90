program thermalization
 ! For command line options
 use kinds
 use cla

 implicit none

 double precision, parameter :: pi = dacos(-1.d0)
 integer, parameter :: phist_bins = 30
 double precision, parameter :: phist_max = 1.5d0  ! GeV

 double precision dt, dx, dz
 integer nt, nx, nz, ix_plot, iz_plot
 integer max_sort
 integer total_ev, total_files
 integer i_file
 character*80 particle_file, datafile_alias, saveload_file
 double precision gs_sigma, many_sigma_sqr, gauss_denom, gauss_norm

 double precision, dimension(:,:,:,:,:,:), allocatable :: Tmn, TmnL ! mu, nu, sort, t,x,z
 double precision, dimension(:,:,:,:,:), allocatable :: phist  ! sort, t,x,z, p
 double precision, dimension(:,:,:,:), allocatable :: temp     ! sort, t,x,z
 double precision, dimension(:,:,:,:), allocatable :: temp_err ! sort, t,x,z
 double precision, dimension(:,:,:,:,:), allocatable :: jmu ! mu, sort, t,x,z
 double precision, dimension(:,:,:,:), allocatable :: jBmu, jSmu ! mu, t,x,z
 double precision, dimension(:,:,:,:,:), allocatable :: umu ! mu, sort, t,x,z
 double precision, dimension(:,:), allocatable :: total_p ! mu, t
 integer, dimension(:,:,:,:), allocatable :: part_num_arr ! sort, t,x,z; no gaussian smearing
 integer, dimension(:), allocatable :: total_B, total_S !t
 logical load_from_saved

 call cla_init
 call cla_register('-load_from_saved',  'load Tmn from previously saved file',  cla_logical  ,'T')
 call cla_register('-urqmd_input',  'urqmd file alias', cla_char, &
                   '/scratch/hyihp/oliiny/therm_project/urqmd_2588090/urqmd-3.4/test.f14')
 call cla_register('-save_load_file',  'file to save/load Tmn, jmu, etc',  cla_char  ,'test')
 call cla_register('-nx', 'grid x dimension (-nx:nx) in dx', cla_int, '10')
 call cla_register('-nz', 'grid z dimension (-nz:nz) in dz', cla_int, '10')
 call cla_register('-nt', 'number of points in time to consider', cla_int, '20')
 call cla_register('-dx', 'grid x step [fm]', cla_float, '1.d0')
 call cla_register('-dz', 'grid z step [fm]', cla_float, '1.d0')
 call cla_validate
 call cla_get('-load_from_saved', load_from_saved)
 call cla_get('-urqmd_input', datafile_alias)
 call cla_get('-save_load_file', saveload_file)

 if (.not. load_from_saved) then
   total_ev = 0
   max_sort = 7
   call cla_get('-nx', nx)
   call cla_get('-nz', nz)
   call cla_get('-nt', nt)
   call cla_get('-dx', dx)
   call cla_get('-dz', dz)

   gs_sigma = 1.d0
   many_sigma_sqr = 3 * 3 * gs_sigma * gs_sigma
   gauss_denom = 2 * gs_sigma * gs_sigma
   gauss_norm = (2 * pi * gs_sigma * gs_sigma)**(-3./2.)

   call init_arrays()
   call system('ls '//trim(datafile_alias)//'|wc -l > file_list.txt')
   call system('ls '//trim(datafile_alias)//' >> file_list.txt')

   open(unit = 20, file = 'file_list.txt')
   read(20,*)total_files
   print *,"Total files: ", total_files
   do i_file = 1, total_files
     read(20,'(A)') particle_file
     call Tmn_from_f14(trim(adjustl(particle_file)))
   end do
   close(20)
   call normalize_to_event_number()
   call save_Tmn(trim(adjustl(saveload_file)))
 else
   call read_Tmn(trim(adjustl(saveload_file)))
 endif

 call print_conserved('conserved_quantities.txt')

 ix_plot = 6
 iz_plot = 0
 call print_Tmn('Tmn.txt', .False., ix_plot, iz_plot) ! FALSE - comp. frame

 call get_Landau_Tmn()

 call print_Tmn('TmnL.txt', .True., ix_plot, iz_plot) ! TRUE - Landau RF
 call print_collective_velocities('v_collective.txt', ix_plot, iz_plot)

 call print_vtk_map('vtk/edens.vtk', "energy_density")
 call print_vtk_map('vtk/dens.vtk',  "density")
 call print_vtk_map('vtk/p.vtk',     "average_pressure")
 call print_vtk_map('vtk/x.vtk',     "pressure_asymetry_x")
 call print_vtk_map('vtk/y.vtk',     "off_diagonality_measure_y")
 call print_vtk_map('vtk/invRe.vtk', "invRe")
 call print_vtk_map('vtk/Npart.vtk', "particle_number")

 call print_var_versus_t('plots/vz_all.dat', 'vz', 0, ix_plot, iz_plot)
 call print_var_versus_t('plots/vz_pi.dat',  'vz', 1, ix_plot, iz_plot)
 call print_var_versus_t('plots/vz_K.dat',   'vz', 2, ix_plot, iz_plot)
 call print_var_versus_t('plots/vz_rho.dat', 'vz', 3, ix_plot, iz_plot)
 call print_var_versus_t('plots/vz_N.dat',   'vz', 4, ix_plot, iz_plot)
 call print_var_versus_t('plots/vz_De.dat',  'vz', 5, ix_plot, iz_plot)
 call print_var_versus_t('plots/e_tot.dat', 'energy_density',   0, ix_plot, iz_plot)
 call print_var_versus_t('plots/p_tot.dat', 'average_pressure', 0, ix_plot, iz_plot)
 call print_var_versus_t('plots/x_tot.dat', 'pressure_asymetry_x', 0, ix_plot, iz_plot)
 call print_var_versus_t('plots/x_pi.dat',  'pressure_asymetry_x', 1, ix_plot, iz_plot)
 call print_var_versus_t('plots/y_tot.dat', 'off_diagonality_measure_y', 0, ix_plot,iz_plot)
 call print_var_versus_t('plots/invRe.dat',  'invRe', 0, ix_plot, iz_plot)

 call print_percentage_of_var_in_range_vs_t('plots/x_area.dat', 'pressure_asymetry_x', 1.d-3, 0.3d0)
 call print_percentage_of_var_in_range_vs_t('plots/y_area.dat', 'off_diagonality_measure_y', 1.d-3, 0.3d0)

 call system('ls '//trim(datafile_alias)//'|wc -l > file_list.txt')
 call system('ls '//trim(datafile_alias)//' >> file_list.txt')

 print *,"Reading UrQMD files again to build dN/dp histograms"
 open(unit = 20, file = 'file_list.txt')
 read(20,*)total_files
 print *,"Total files: ", total_files
 do i_file = 1, total_files
   read(20,'(A)') particle_file
   call get_phist(trim(adjustl(particle_file)))
 end do
 close(20)
 call print_phist('plots/test_phist.dat', 20, 0, 0, 1)  ! pions at (0,0) at 20*dt

 call delete_arrays_from_memory()

contains

subroutine init_arrays()
  allocate(Tmn(0:3, 0:3, 0:max_sort,  1:nt, -nx:nx, -nz:nz))
  allocate(TmnL(0:3, 0:3, 0:max_sort,  1:nt, -nx:nx, -nz:nz))
  allocate(phist(0:max_sort,  1:nt, -nx:nx, -nz:nz, 0:phist_bins))
  allocate(temp(0:max_sort,  1:nt, -nx:nx, -nz:nz))
  allocate(temp_err(0:max_sort,  1:nt, -nx:nx, -nz:nz))
  allocate(jmu(0:3, 0:max_sort,  1:nt, -nx:nx, -nz:nz))
  allocate(jBmu(0:3,  1:nt, -nx:nx, -nz:nz))
  allocate(jSmu(0:3,  1:nt, -nx:nx, -nz:nz))
  allocate(umu(0:3, 0:max_sort,  1:nt, -nx:nx, -nz:nz))
  allocate(total_p(0:3, 1:nt))
  allocate(total_B(1:nt))
  allocate(total_S(1:nt))
  allocate(part_num_arr(0:max_sort, 1:nt, -nx:nx, -nz:nz))
  Tmn     = 0.d0
  TmnL    = 0.d0
  jBmu    = 0.d0
  jSmu    = 0.d0
  jmu     = 0.d0
  phist   = 0.d0
  part_num_arr = 0
  total_p(0:3,:) = 0.d0
  total_B(:) = 0
  total_S(:) = 0
  umu(0,:,:,:,:) = 1.d0
  umu(1:3,:,:,:,:) = 0.d0
end subroutine

subroutine delete_arrays_from_memory()
  deallocate(Tmn)
  deallocate(TmnL)
  deallocate(phist)
  deallocate(temp)
  deallocate(temp_err)
  deallocate(jmu)
  deallocate(jBmu)
  deallocate(jSmu)
  deallocate(umu)
  deallocate(total_p)
  deallocate(total_B)
  deallocate(total_S)
  deallocate(part_num_arr)
end subroutine

subroutine save_Tmn(fname)
  character(len=*), intent(in) :: fname
  print *, "Saving Tmn data to file ", fname
  open(unit = 9, file = fname, form='unformatted')
    write(9)total_ev, max_sort, nt, nx, nz, dt, dx, dz, gs_sigma
    write(9)Tmn
    write(9)jmu
    write(9)jBmu
    write(9)jSmu
    write(9)part_num_arr
  close(9)
end subroutine

subroutine read_Tmn(fname)
  character(len=*), intent(in) :: fname
  print *, "Reading Tmn data to file ", fname
  open(unit = 9, file = fname, form='unformatted')
    read(9)total_ev, max_sort, nt, nx, nz, dt, dx, dz, gs_sigma
    many_sigma_sqr = 3 * 3 * gs_sigma * gs_sigma
    gauss_denom = 2 * gs_sigma * gs_sigma
    gauss_norm = (2 * pi * gs_sigma * gs_sigma)**(-3./2.)

    print *,"Total events: ", total_ev
    print *,"nt, nx, nz: ", nt, nx, nz
    print *,"dt, dx, dz: ", dt, dx, dz
    call init_arrays()
    read(9)Tmn
    read(9)jmu
    read(9)jBmu
    read(9)jSmu
    read(9)part_num_arr
  close(9)
end subroutine

subroutine get_Landau_Tmn()
 use Land_Eck, only: FindLandau
 implicit none

 integer it, sort, ix, iz
 do it = 1, nt
   do sort = 0, max_sort
     do ix = -nx, nx; do iz = -nz, nz
       call FindLandau( Tmn(0:3, 0:3, sort, it, ix, iz),&
                       TmnL(0:3, 0:3, sort, it, ix, iz),&
                       umu(0:3, sort, it, ix, iz))
     end do; end do
   end do
 end do
end subroutine

subroutine Tmn_from_f14(fname)
 character(len=*), intent(in) :: fname
 double precision Elab, r(0:3), p(0:3), m, dr(1:3), sf, upart(0:3)
 integer tsteps, ev, Npart, i, nu, ityp, i3, sort, Bpart, Spart, io
 integer ch, last_col_part
 integer it, ix, iz

 print *, "Reading from ", fname
 open(unit = 14, file = fname)
 do ! event cycle
   if (.not. read_f14_event_header(14, Elab, ev, tsteps, dt)) then
     exit
   endif
   if (mod(ev, 20) == 0) then
     print *, "Event number: ", ev, " Elab[AGeV]: ", Elab, "tsteps: ", tsteps
   endif
   do it = 1, tsteps
     read(14,*) Npart
     read(14,*)
     do i = 1, Npart

       read(14,*, iostat = io) r(0:3), p(0:3), m, ityp, i3, ch, last_col_part
       if (io .ne. 0) then
         print *, "error reading file, io = ", io, " ev = ", ev, " i = ", i
       endif
       if (it > nt) then; cycle; endif
       Bpart = BfromItyp(ityp)
       Spart = SfromItyp(ityp)
       total_p(0:3, it) = total_p(0:3, it) + p(0:3)
       total_B(it) = total_B(it) + Bpart
       total_S(it) = total_S(it) + Spart
       sort = get_sort(ityp)

       ! Here all particles are counted, including spectators
       do ix = -nx, nx; do iz = -nz, nz
         dr(1) = r(1) - ix * dx
         dr(2) = r(2)
         dr(3) = r(3) - iz * dz
         if (too_far(dr(1:3))) then; cycle; endif

         part_num_arr(0, it, ix, iz) = part_num_arr(0, it, ix, iz) + 1
         part_num_arr(sort, it, ix, iz) = part_num_arr(sort, it, ix, iz) + 1
       end do; end do

       ! Do not add spectators to anything else except total particle count
       if (last_col_part == 0) then; cycle; endif

       upart(0) = 1.d0
       upart(1:3) = p(1:3)/p(0)

       do ix = -nx, nx; do iz = -nz, nz ! loop over space grid
         ! dr - comp. frame vector from grid point to particle
         dr(1) = r(1) - ix * dx
         dr(2) = r(2)
         dr(3) = r(3) - iz * dz
         if (too_far(dr(1:3))) then; cycle; endif
         sf = smearing_factor(dr(1:3), p(0:3))

         jBmu(0:3, it, ix, iz) = jBmu(0:3, it, ix, iz) +&
                                       upart(0:3) * sf * Bpart
         jSmu(0:3, it, ix, iz) = jSmu(0:3, it, ix, iz) +&
                                       upart(0:3) * sf * Spart
         jmu(0:3, 0, it, ix, iz) = jmu(0:3, 0, it, ix, iz) +&
                                       upart(0:3) * sf
         part_num_arr(0, it, ix, iz) = part_num_arr(0, it, ix, iz) + 1
         do nu = 0,3
           Tmn(0:3, nu, 0, it, ix, iz) = Tmn(0:3, nu, 0, it, ix, iz) +&
                                         upart(0:3) * p(nu) * sf
         end do

         ! Sort specific analysis
         if (sort > 0) then
           jmu(0:3, sort, it, ix, iz) = jmu(0:3, sort, it, ix, iz) +&
                                            upart(0:3) * sf
           do nu = 0,3
             Tmn(0:3, nu, sort, it, ix, iz) = Tmn(0:3, nu, sort, it, ix, iz) +&
                                              upart(0:3) * p(nu) * sf
           end do
         endif ! end sort specific analysis
       end do; end do ! end loop over space grid
     end do
   end do
 end do ! end event cycle
 close(14)

 total_ev = total_ev + ev
end subroutine

subroutine normalize_to_event_number()
 print *,"normalazing to total number of events: ", total_ev
 Tmn = Tmn / total_ev
 jBmu = jBmu / total_ev
 jSmu = jSmu / total_ev
 jmu  = jmu  / total_ev
 total_p = total_p / total_ev
 total_B = total_B / total_ev
 total_S = total_S / total_ev
end subroutine

double precision function smearing_factor(dr, p)
  double precision, intent(in) :: dr(1:3), p(0:3)
  double precision gam_inv, bet(1:3), tmp, dr_RF(1:3), dr_RF_sqr

  bet(1:3) = p(1:3)/p(0)
  gam_inv = 1.d0 - bet(1)*bet(1) - bet(2)*bet(2) - bet(3)*bet(3)
  if (gam_inv > 0.d0) then
    gam_inv = sqrt(gam_inv)
  else
    gam_inv = huge(gam_inv)/2.d0 ! Very large number
    print *,"Warning: got weird particle - p = ", p
  endif
  tmp = (bet(1)*dr(1) + bet(2)*dr(2) + bet(3)*dr(3))/gam_inv/(1.d0 + gam_inv)
  dr_RF(1:3) = dr(1:3) + tmp * bet(1:3)
  dr_RF_sqr = dr_RF(1)*dr_RF(1) + dr_RF(2)*dr_RF(2) + dr_RF(3)*dr_RF(3)
  smearing_factor = gauss_norm * exp(-dr_RF_sqr/gauss_denom)
end function smearing_factor

subroutine print_conserved(fname)
  character(len=*), intent(in) :: fname
  integer it

  open(unit = 8, file = fname)
  write(8,*)"# t[fm/c] E[GeV] px[GeV/c] py[GeV/c] pz[GeV/c]    B     S"
  do it = 1, nt
    write(8,*) it*dt, total_p(0:3, it), total_B(it), total_S(it)
  end do
  close(8)
end subroutine

subroutine print_vtk_map(fname, variable_to_print)
  character(len=*), intent(in) :: fname, variable_to_print
  character(len=4) s_hlp
  integer it, ix, iz
  double precision var

  do it = 1, nt
    write(s_hlp,'(i4)')it
    open(unit = 8, file = fname//"."//trim(adjustl(s_hlp)))
    write(8,'(A)')"# vtk DataFile Version 2.0"
    write(8,'(A)')variable_to_print
    write(8,'(A)')"ASCII"
    write(8,'(A)')"DATASET STRUCTURED_POINTS"
    write(8,'(A11,3I5)')"DIMENSIONS ", 2*nx+1, 1, 2*nz+1
    write(8,'(A8,3I5)')"SPACING ", 1, 1, 1
    write(8,'(A7,3I5)')"ORIGIN ", nx, 0, nz
    write(8,'(A,I8)')"POINT_DATA", (2*nx+1) * (2*nz+1)
    write(8,'(A,A,A)')"SCALARS ", variable_to_print, " float 1"
    write(8,'(A)')"LOOKUP_TABLE default"

    do iz = -nz, nz
      do ix = -nx, nx
        var = select_var(variable_to_print, 0,it,ix,iz)
        write(8,'(f16.4)', advance = 'no') var
      end do
      write(8,*)
    end do
    close(8)
  end do

end subroutine

subroutine print_var_versus_t(fname, variable_to_print, sort, ix, iz)
  character(len=*), intent(in) :: fname, variable_to_print
  integer, intent(in) :: sort, ix, iz
  integer it
  character(len=12)sname

  call get_sort_name(sort, sname)

  open(unit = 8, file = fname)

  write(8,'(A,A)')"# ",sname
  write(8,'(A,A,2(A,F5.1),A)')"# ", variable_to_print, &
                              " @ x = ", ix*dx, " fm, z = ", iz*dz, " fm"
  write(8,'(A,A)')"# time[fm/c] ", variable_to_print
  do it = 1, nt
    write(8,'(f5.1,f10.4)') it*dt, select_var(variable_to_print, sort,it,ix,iz)
  end do
  close(8)

end subroutine

subroutine print_percentage_of_var_in_range_vs_t(fname, varname, r_down, r_up)
  ! prints #cells(r_down < var < r_up) / #cells(r_down < var) versus time
  character(len=*), intent(in) :: fname, varname
  double precision, intent(in) :: r_down, r_up
  double precision n1, n2, var
  integer ix, iz, it

  open(unit = 8, file = fname)
  write(8,'(3A,F5.2)')"# time[fm/c] % of ", varname, " below ", r_up
  do it = 1, nt
    n1 = 0.d0; n2=0.d0
    do iz = -nz, nz; do ix = -nx, nx
      var = select_var(varname, 0, it, ix, iz)
      if (var > r_down) then
        n1 = n1 + 1.d0
        if (var < r_up) then
          n2 = n2 + 1.d0
        endif
      endif
    end do; end do
    if (n1 > 0.d0) then
      write(8,'(f5.1,f10.4)') it*dt, n2/n1*1.d2 ! in %
    else
      write(8,'(f5.1,f10.4)') it*dt, 0.d0
    endif
  end do
  close(8)
end subroutine

double precision function select_var(varname, sort,it,ix,iz) result (var)
  use Land_Eck, only: EuclidProduct
  integer, intent(in) :: sort,it,ix,iz
  character(len=*), intent(in) :: varname

  select case(varname)
    case ("energy_density")
      var = get_land_e0(sort,it,ix,iz)
    case ("average_pressure")
      var = get_land_p0(sort,it,ix,iz)
    case ("density")
      var = EuclidProduct(jmu(0:3,sort,it,ix,iz), umu(0:3,sort,it,ix,iz))
    case ("particle_number")
      var = part_num_arr(sort,it,ix,iz)
    case ("pressure_asymetry_x")
      var = get_asym_x(sort,it,ix,iz)
    case ("off_diagonality_measure_y")
      var = get_offdiag_y(sort,it,ix,iz)
    case ("invRe")
      var = get_inv_Re(sort,it,ix,iz)
    ! for all velocities: minus sign, because u_mu has lower index
    case ("vx")
      var = -umu(1,sort,it,ix,iz)/umu(0,sort,it,ix,iz)
    case ("vy")
      var = -umu(2,sort,it,ix,iz)/umu(0,sort,it,ix,iz)
    case ("vz")
      var = -umu(3,sort,it,ix,iz)/umu(0,sort,it,ix,iz)
    case default
      print *,"Wrong variable name: ", varname
  end select

end function select_var

double precision function get_land_e0(sort,it,ix,iz)
  integer, intent(in) :: sort,it,ix,iz
  get_land_e0 = TmnL(0,0,sort,it,ix,iz)
end function get_land_e0

double precision function get_land_p0(sort,it,ix,iz)
  integer, intent(in) :: sort,it,ix,iz
  double precision px, py, pz
  px = TmnL(1,1,sort,it,ix,iz)
  py = TmnL(2,2,sort,it,ix,iz)
  pz = TmnL(3,3,sort,it,ix,iz)
  get_land_p0 = (px + py + pz) / 3.d0
end function get_land_p0

double precision function get_asym_x(sort,it,ix,iz)
  integer, intent(in) :: sort,it,ix,iz
  double precision px, py, pz
  px = TmnL(1,1,sort,it,ix,iz)
  py = TmnL(2,2,sort,it,ix,iz)
  pz = TmnL(3,3,sort,it,ix,iz)
  if (px + py + pz > 1.d-4) then
    get_asym_x = (abs(px-py) + abs(py-pz) + abs(pz-px))/(px + py + pz)
  else
    get_asym_x = 0.d0
  endif
end function get_asym_x

double precision function get_offdiag_y(sort,it,ix,iz)
  integer, intent(in) :: sort,it,ix,iz
  double precision px, py, pz, Txy, Txz, Tyz
  px = TmnL(1,1,sort,it,ix,iz)
  py = TmnL(2,2,sort,it,ix,iz)
  pz = TmnL(3,3,sort,it,ix,iz)
  Txy = TmnL(1,2,sort,it,ix,iz)
  Txz = TmnL(1,3,sort,it,ix,iz)
  Tyz = TmnL(2,3,sort,it,ix,iz)
  if (px + py + pz > 1.d-4) then
    get_offdiag_y = (abs(Txy) + abs(Tyz) + abs(Txz)) / (px + py + pz) * 3.d0
  else
    get_offdiag_y = 0.d0
  endif
end function get_offdiag_y

double precision function get_inv_Re(sort,it,ix,iz)
  integer, intent(in) :: sort,it,ix,iz
  double precision px, py, pz, Txy, Txz, Tyz, p0
  px = TmnL(1,1,sort,it,ix,iz)
  py = TmnL(2,2,sort,it,ix,iz)
  pz = TmnL(3,3,sort,it,ix,iz)
  Txy = TmnL(1,2,sort,it,ix,iz)
  Txz = TmnL(1,3,sort,it,ix,iz)
  Tyz = TmnL(2,3,sort,it,ix,iz)
  p0 = (px + py + pz) / 3.d0
  if (p0 > 1.d-9) then
    get_inv_Re = sqrt(Txy*Txy + Tyz*Tyz + Txz*Txz + (px-p0)*(px-p0) +&
                 (py-p0)*(py-p0) + (pz-p0)*(pz-p0)) / p0
  else
    ! Very big Re^-1 for vacuum: definitely kinetics, not hydro there
    get_inv_Re = 1.d3
  endif
end function get_inv_Re



subroutine print_collective_velocities(fname, ix, iz)
  character(len=*), intent(in) :: fname
  integer, intent(in) :: ix, iz
  integer it, sort
  character(len=12)sname
  double precision u_hlp(0:3)

  open(unit = 8, file = fname)
  write(8,'(2(A,f5.1))')"# Collective velocities of particles at x = ",&
                                                 ix*dx, ", z = ", iz*dz
  do it = 1, nt
    write(8,'(A,f5.1,A)')"t = ", it*dt, " fm/c"
    do sort = 0, max_sort
      call get_sort_name(sort, sname)
      u_hlp(0:3) = umu(0:3, sort, it, ix, iz)
      ! u_mu has lower index, so vector v = -u(1:3)/u(0)
      write(8,'(3f9.3,A4,A)') -u_hlp(1:3)/u_hlp(0),"  # ",trim(sname)
    end do
  end do
  close(8)
end subroutine

subroutine print_Tmn(fname, landau, ix, iz)
  character(len=*), intent(in) :: fname
  logical, intent(in) :: landau
  integer, intent(in) :: ix, iz
  integer it, sort, mu, nu
  double precision T_hlp(0:3,0:3)
  character(len=12)sname

  open(unit = 8, file = fname)
  if (landau) then
    write(8,'(2(A,f5.1))')"# Tmn in Landau rest frame at x = ",&
                                          ix*dx, ", z = ", iz*dz
  else
    write(8,'(2(A,f5.1))')"# Tmn in computational frame at x = ",&
                                            ix*dx, ", z = ", iz*dz
  endif

  do sort = 0, max_sort
    call get_sort_name(sort, sname)
    write(8,'(A2,A)')"# ",trim(sname)
    do it = 1, nt
      write(8,'(A,f5.1,A)')"t = ", it*dt, " fm/c"
      if (landau) then
        T_hlp(0:3,0:3) = TmnL(0:3,0:3, sort, it, ix, iz)
      else
        T_hlp(0:3,0:3) = Tmn(0:3,0:3, sort, it, ix, iz)
      endif
      do mu = 0, 3
        write(8,'(4f10.4)') (T_hlp(mu,nu), nu = 0,3)
      end do
      write(8,*)
    end do
  end do
  close(8)
end subroutine

integer function get_sort(ityp)
  integer, intent(in) :: ityp
  select case(ityp)
    case (101); get_sort = 1 ! pions
    case (106); get_sort = 2 ! kaons
    case (104); get_sort = 3 ! rhos
    case (1);   get_sort = 4 ! nucleons
    case (17);  get_sort = 5 ! Deltas
    case (27);  get_sort = 6 ! Lambdas
    case (102); get_sort = 7 ! etas
    case default; get_sort = -1
  end select
end function get_sort

subroutine get_sort_name(sort, sname)
  integer, intent(in) :: sort
  character(len=12), intent(out) :: sname

  select case(sort)
    case (0); sname = "total"
    case (1); sname = "pions"
    case (2); sname = "kaons"
    case (3); sname = "rhos"
    case (4); sname = "nucleons"
    case (5); sname = "Deltas"
    case (6); sname = "Lambdas"
    case (7); sname = "etas"
    case default; print *,"error: unknown sort ", sort
  end select
end subroutine

double precision function get_sort_mass(sort) result (mass)
  integer, intent(in) :: sort
  select case(sort)
    case (1); mass = 0.138d0
    case (2); mass = 0.495d0
    case (3); mass = 0.776d0
    case (4); mass = 0.938d0
    case (5); mass = 1.232d0
    case (6); mass = 1.115d0
    case (7); mass = 0.548d0
    case default; print *, "error: unknown sort ", sort
  end select
end function get_sort_mass

logical function too_far(dr)
  double precision, intent(in) :: dr(1:3)
  if (dr(1)*dr(1) + dr(2)*dr(2) + dr(3)*dr(3) > many_sigma_sqr) then
    too_far = .TRUE.
  else
    too_far = .FALSE.
  endif
end function too_far

logical function read_f14_event_header(uread, Elab, ev, tsteps, dt)
  ! Reads UrQMD file 14 header, returns true at success
  ! uread - number of unit to be read
  ! Elab - collision energy
  ! ev - event number
  ! tsteps - number of output moments per event in f14
  ! dt - time difference between output moments
  integer, intent(in) :: uread
  double precision, intent(out) :: Elab, dt
  integer, intent(out) :: tsteps, ev

  character testchar
  character*100 tmp
  integer i, io

  read(uread,*,iostat = io) testchar;
  if ((io .ne. 0) .or. (testchar .ne. 'U')) then
    read_f14_event_header = .FALSE.
    return
  endif

  do i=1,3; read(uread,*); end do
  read(uread,'(A39, e11.4)') tmp, Elab
  read(uread,'(a7,i9,a54,i7,a20,f11.3)') tmp, ev, tmp, tsteps, tmp, dt
  tsteps = nint(tsteps/dt)
  do i=1,11; read(uread,*); end do
  read_f14_event_header = .TRUE.

end function read_f14_event_header

integer function BfromItyp(ityp) result (Bch)
 implicit none
 integer, intent(in) :: ityp

    if (abs(ityp) .ge. 100) then
      Bch = 0
    else if (ityp > 0) then
      Bch = 1
    else
      Bch = -1
    endif
    return

end function BfromItyp

integer function SfromItyp(ityp) result (Sch)
 implicit none
 integer, intent(in) :: ityp
 integer ia

 integer :: strres(1:55), strmes(100:137)
 data strres/ 26*0,13*1,9*1,6*2,3/
 data strmes/ 6*0,-1,0,-1,0,-1,0,0,-1,0,0,0,-1,0,0,0,&
                  -1,0,0,0,-1,0,0,0,-1,0,0,0,0,0,0,0,0/

  !return zero strangeness in case of PYTHIA PDG particles
  if (abs(ityp).gt.1000)then
      Sch=0
      return
  endif

  !standart UrQMD particles
  if(ityp.eq.0) then;  Sch=0; return;  endif
  ia=iabs(ityp)
  if(ia .ge. 100)then;  Sch = strmes(ia);  else;   Sch = strres(ia);  endif
  if (ityp < 0) then; Sch = - Sch; endif
  return

end function SfromItyp

integer function histindex_from_p(p) result (histindex)
 implicit none
 double precision, intent(in) :: p
 histindex = floor(p/phist_max*phist_bins)
end function histindex_from_p

double precision function p_from_histindex(ind) result (p)
 implicit none
 integer, intent(in) :: ind
 p = (ind + 0.5d0) * phist_max / phist_bins
end function p_from_histindex

subroutine print_phist(fname, it, ix, iz, sort)
 implicit none
 character(len=*), intent(in) :: fname
 integer, intent(in) :: it, ix, iz, sort
 integer i
 double precision p, dN_dp_over_p2, m, E, bin_width

 print *,"Writing momentum histograms to file ", fname
 open(unit = 8, file = fname)
 m = get_sort_mass(sort)
 bin_width = phist_max / phist_bins
 do i = 0, phist_bins
   p = p_from_histindex(i)
   E = sqrt(p*p + m*m)
   dN_dp_over_p2 = phist(sort, it, ix, iz, i) / bin_width / p / p
   write(8,'(2f12.5)') E, dN_dp_over_p2
 end do
 close(8)
end subroutine

subroutine get_phist(fname)
 use Land_Eck, only: GetBoostMatrix
 implicit none
 character(len=*), intent(in) :: fname
 double precision Elab, r(0:3), p(0:3), m, dr(1:3), sf, mom_abs
 double precision M_boost(0:3,0:3), p_rest(0:3)
 integer tsteps, ev, Npart, i, ityp, i3, sort, io
 integer ch, last_col_part
 integer it, ix, iz, hist_index

 print *, "Momentum histograms building: reading from ", fname
 open(unit = 14, file = fname)
 do ! event cycle
   if (.not. read_f14_event_header(14, Elab, ev, tsteps, dt)) then
     exit
   endif
   if (mod(ev, 20) == 0) then
     print *, "Event number: ", ev, " Elab[AGeV]: ", Elab, "tsteps: ", tsteps
   endif
   do it = 1, tsteps
     read(14,*) Npart
     read(14,*)
     do i = 1, Npart

       read(14,*, iostat = io) r(0:3), p(0:3), m, ityp, i3, ch, last_col_part
       if (io .ne. 0) then
         print *, "error reading file, io = ", io, " ev = ", ev, " i = ", i
       endif
       if (it > nt) then; cycle; endif
       sort = get_sort(ityp)
       if (sort < 0) then; cycle; endif

       ! Do not add spectators to histogram
       if (last_col_part == 0) then; cycle; endif

       do ix = -nx, nx; do iz = -nz, nz ! loop over space grid
         ! dr - comp. frame vector from grid point to particle
         dr(1) = r(1) - ix * dx
         dr(2) = r(2)
         dr(3) = r(3) - iz * dz
         if (too_far(dr(1:3))) then; cycle; endif
         sf = smearing_factor(dr(1:3), p(0:3))

         ! Add particle to momentum histogram
         call GetBoostMatrix(umu(0:3, sort,  it, ix, iz), M_boost)
         p_rest = matmul(M_boost, p(0:3))  ! boost to local Landau RF
         mom_abs = sqrt(p_rest(1)*p_rest(1) + &
                        p_rest(2)*p_rest(2) + &
                        p_rest(3)*p_rest(3))
         hist_index = histindex_from_p(mom_abs)
         if (hist_index .le. phist_bins) then
           ! print *, "Adding ", ityp, " with momentum ", p_rest(0:3),&
           !         " to cell (ix,iz) = (", ix, ",", iz, ") at time it = ", it
           phist(sort, it, ix, iz, hist_index) = &
            phist(sort, it, ix, iz, hist_index) + sf
         endif
       end do; end do ! end loop over space grid
     end do
   end do
 end do ! end event cycle
 close(14)
end subroutine



end program thermalization
