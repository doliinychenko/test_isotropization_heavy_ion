program thermalization
 implicit none

 real dt, dx, dz
 real gs_sigma, many_sigma_sqr, gauss_denom
 integer nt, nx, nz
 integer max_sort

 real, dimension(:,:,:,:,:,:), allocatable :: Tmn, TmnL ! mu, nu, sort, t,x,z
 real, dimension(:,:,:,:,:), allocatable :: jBmu, jSmu, umu ! mu, sort, t,x,z
 real, dimension(:,:), allocatable :: total_p ! mu, t
 integer, dimension(:), allocatable :: total_B, total_S !t

 max_sort = 7
 dx = 1.0
 dz = 1.0
 nt = 10
 nx = 10
 nz = 10
 gs_sigma = 1.0
 many_sigma_sqr = 4 * 4 * gs_sigma * gs_sigma
 gauss_denom = 2 * gs_sigma * gs_sigma

 call init_arrays()
 call Tmn_from_f14('test.f14')
 call print_conserved('conserved_quantities.txt')

contains

subroutine init_arrays()
  allocate(Tmn(0:3, 0:3, 1:max_sort,  1:nt, 0:nx, 0:nz))
  allocate(TmnL(0:3, 0:3, 1:max_sort,  1:nt, 0:nx, 0:nz))
  allocate(jBmu(0:3, 1:max_sort,  1:nt, 0:nx, 0:nz))
  allocate(jSmu(0:3, 1:max_sort,  1:nt, 0:nx, 0:nz))
  allocate(umu(0:3, 1:max_sort,  1:nt, 0:nx, 0:nz))
  allocate(total_p(0:3, 1:nt))
  allocate(total_B(1:nt))
  allocate(total_S(1:nt))
  Tmn     = 0.0
  TmnL    = 0.0
  jBmu    = 0.0
  jSmu    = 0.0
  total_p(0:3,:) = 0.0
  total_B(:) = 0
  total_S(:) = 0
  umu(0,:,:,:,:) = 1.0
  umu(1:3,:,:,:,:) = 0.0
end subroutine

subroutine Tmn_from_f14(fname)
 character(len=*), intent(in) :: fname
 real Elab, r(0:3), p(0:3), m, dr(1:3), sf, upart(0:3)
 integer tsteps, ev, Npart, i, ityp, i3, sort, Bpart, Spart
 integer it, ix, iz

 open(unit = 14, file = 'test.f14')
 do ! event cycle
   if (.not. read_f14_event_header(14, Elab, ev, tsteps, dt)) then
     exit
   endif
   print *, "Event number: ", ev, " Elab[AGeV]: ", Elab, "tsteps: ", tsteps
   do it = 1, tsteps
     read(14,*) Npart
     read(14,*)
     do i = 1, Npart
       read(14,*) r(0:3), p(0:3), m, ityp, i3
       if (it > nt) then; cycle; endif
       ! Analysis
       sort = get_sort(ityp, i3)
       Bpart = BfromItyp(ityp)
       Spart = SfromItyp(ityp)
       total_p(0:3, it) = total_p(0:3, it) + p(0:3)
       total_B(it) = total_B(it) + Bpart
       total_S(it) = total_S(it) + Spart

       if (sort < 0) then; cycle; endif
       do ix = 1, nx; do iz = 1, nz ! loop over space grid
         ! dr - comp. frame vector from grid point to particle
         dr(1) = r(1) - ix * dx
         dr(2) = r(2)
         dr(3) = r(3) - iz * dz
         if (too_far(dr(1:3))) then; cycle; endif
         sf = smearing_factor(dr(1:3), p(0:3))
         upart(0) = 1.0
         upart(1:3) = p(1:3)/p(0)
         do nu = 0,3
           Tmn(0:3, nu, sort, it, ix, iz) = Tmn(0:3, nu, sort, it, ix, iz) +&
                                            upart(0:3) * p(nu) * sf
         end do
         jBmu(0:3) = jBmu(0:3) + upart(0:3) * sf * Bpart
         jSmu(0:3) = jSmu(0:3) + upart(0:3) * sf * Spart
       end do; end do ! end loop over space grid
     end do
   end do
 end do ! end event cycle
 close(14)

end subroutine

real function smearing_factor(dr, p)
  real, intent(in) :: dr(1:3), p(0:3)
  smearing_factor = 1.0 !dummy return, TODO: write the function
end function smearing_factor

subroutine print_conserved(fname)
  character(len=*) fname
  integer it

  open(unit = 8, file = fname)
  write(8,*)"# t[fm/c] E[GeV] px[GeV/c] py[GeV/c] pz[GeV/c]    B     S"
  do it = 1, nt
    write(8,*) it*dt, total_p(0:3, it), total_B(it), total_S(it)
  end do
  close(8)
end subroutine

integer function get_sort(ityp, i3)
  integer, intent(in) :: ityp, i3
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

logical function too_far(dr)
  real, intent(in) :: dr(1:3)
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
  real, intent(out) :: Elab, dt
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


end program thermalization
