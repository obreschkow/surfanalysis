module module_asymmetries

use shared_module_core
use shared_module_arguments
use shared_module_hdf5
use shared_module_lapack
use shared_module_vectors
use shared_module_maths
use shared_module_constants
use shared_module_sort
use module_global
use module_io
use module_gethalo
use module_processhalo
use module_showhalo

implicit none

private
public   :: task_asymmetries

type type_properties
   integer*8   :: id
   integer*8   :: ntot,ngas,ncdm ! number of particles
   real*4      :: mass
   real*4      :: fbaryon
   real*4      :: asymmetry ! dimensionless asymmetry parameter at initial snapshot
   real*4      :: ax,ay,az
   real*4      :: radius_rms_initial,radius_rms_final ! median radius at both snapshots
   real*4      :: lambda_initial,lambda_final ! spin parameter at both snapshots
   real*4      :: jcdm_initial,jcdm_final ! norm of CDM AM at both snapshots
   real*4      :: jgas_initial,jgas_final ! norm of gas AM at both snapshots
end type type_properties

real*4,parameter     :: fraction = 1.0 ! fraction of halos to analyze
real*4,parameter     :: h = 0.6751
real*4,parameter     :: OmegaM = 0.3121
real*4,parameter     :: OmegaL = 1-OmegaM
real*4,parameter     :: fb = 0.1573 ! baryon fraction relative to all matter
integer*4            :: initial_snapshot = 000 ! snapshot at which to evaluate initial properties
real*4               :: mass_min ! [Msun/h]
real*4,parameter     :: mass_max = 1e16 ! [Msun/h]
   
real*4               :: mpart ! [Msun/h] particle mass
logical              :: subhalos


contains

subroutine task_asymmetries

   implicit none
   
   ! read options
   call get_option_value(subhalos,'-subhalos',.true.)
   call get_option_value(initial_snapshot,'-initial',0)
   call get_option_value(mass_min,'-massmin',1e12)
   call require_no_options_left
   
   mpart = para%L**3/real(para%N,4)**3*OmegaM*2.774495e11 ! [Msun/h] particle mass (in the case of hydro runs, this is the mass of 1 CDM+1 gas particle)
   
   call tic()
   call compute_asymmetry_parameters
   call toc()
   
end subroutine task_asymmetries

subroutine compute_asymmetry_parameters
   
   implicit none
   
   real*4                              :: mpartgas,mpartcdm
   integer*4                           :: nmin,nmax
   type(type_halo),allocatable         :: halo(:)
   integer*4                           :: nhalos
   integer*4                           :: npart,ipart
   type(type_properties),allocatable   :: store(:)
   integer*4                           :: i,j,k
   integer*4                           :: ntot
   real*4                              :: f_complete,f_complete_next
   integer*4                           :: ngas,ncdm
   real*4                              :: mgas,mcdm,mtot ! [Msun/h] total masses
   integer*4                           :: ix,iy ! number of CDM and gas particles
   real*4,allocatable                  :: x(:,:),y(:,:) ! positions of CDM and gas particles
   real*4,allocatable                  :: vx(:,:),vy(:,:) ! velocities of CDM and gas particles
   integer*4                           :: mode,dim
   real*4                              :: a(3)
   real*4                              :: cm(3) ! center of mass
   real*4                              :: lambda ! spin parameter
   real*4                              :: scale ! scale factor of loaded snapshot
   real*4                              :: z,E,Hubble,Delta,lambda_factor
   real*4                              :: collapse_factor
   real*4                              :: radius_rms
   real*4                              :: jcdm(3),jgas(3),jmean
   character(len=255)                  :: fn
   character(len=3)                    :: snstr
   
   ! make particle masses
   mpartgas = fb*mpart     ! [Msun/h]
   mpartcdm = (1-fb)*mpart ! [Msun/h]
   
   ! make allowed range of number or particles
   nmin = floor(mass_min/mpartCDM)     ! minimum possible number of particles in halo of the chosen mass range
   nmax = ceiling(mass_max/mpartgas)   ! maximum possible number of particles in halo of the chosen mass range
   
   ! load global properties of all haloes
   call load_all_halo_properties(halo)
   nhalos = size(halo)
   
   ! subselect halos
   nhalos = int(nhalos*fraction,4)
   
   ! compute total number of particles
   npart = sum(halo(1:nhalos)%npart)+sum(halo(1:nhalos)%npartsub)
   
   allocate(store(nhalos))
      
   k = 0
   ipart = 0
   f_complete_next = 0
   
   do i = 1,nhalos
   
      ntot = halo(i)%npart+halo(i)%npartsub
      ipart = ipart+ntot
      f_complete = real(real(ipart,8)/real(npart,8),4)
      
      if (f_complete>=f_complete_next) then
         call progress(f_complete)
         f_complete_next = f_complete_next+0.01
      end if
      
      if (halo(i)%parentid==-1) then ! selecting parent haloes
         
         if ((ntot>=nmin).and.(ntot<=nmax)) then
   
            call load_halo_particles(i,subhalos)
            
            ! number of particles
            ngas = count(p%species==1)
            ncdm = count(p%species==2)
            
            ! masses
            mgas = mpartgas*ngas ! [Msun/h]
            mcdm = mpartCDM*ncdm ! [Msun/h]
            mtot = mgas+mcdm     ! [Msun/h]
            
            if ((mtot>=mass_min).and.(mtot<=mass_max)) then
            
               if ((ngas>=20).and.(ncdm>=20)) then ! require at least 20 particles of each species
               
                  k = k+1
               
                  ! save halo id
                  store(k)%id = i
                  
                  ! save number of particles
                  store(k)%ngas = ngas
                  store(k)%ncdm = ncdm
                  store(k)%ntot = ngas+ncdm
                  
                  ! compute and save masses of group
                  store(k)%mass = mtot ! [Msun/h]
                  store(k)%fbaryon = mgas/mtot
                  
                  if (allocated(x)) deallocate(x); allocate(x(ncdm,3))
                  if (allocated(y)) deallocate(y); allocate(y(ngas,3))
                  if (allocated(vx)) deallocate(vx); allocate(vx(ncdm,3))
                  if (allocated(vy)) deallocate(vy); allocate(vy(ngas,3))
                  
                  do mode = 1,2
                  
                     if (mode==1) then
                     
                        scale = scalefactor(para%snapshot)
                  
                     else
                     
                        scale = scalefactor(initial_snapshot)
                        call load_halo_particles_at_different_snapshot(i,subhalos,initial_snapshot,.true.)
                  
                     end if
                  
                     ! center particles to geometric centre
                     call center_particles ! necessary, even if particles are centered again later, because of wrapping

                     ! separate CDM and gas particles
                     ix = 0
                     iy = 0
                     do j = 1,int(nparticles,4)
                        if (p(j)%species==2) then ! DM particles
                           ix = ix+1
                           x(ix,:) = p(j)%x
                           vx(ix,:) = p(j)%v
                        else ! gas particles
                           iy = iy+1
                           y(iy,:) = p(j)%x
                           vy(iy,:) = p(j)%v
                        end if
                     end do
                  
                     ! center particles and velocities to CM
                     do dim = 1,3
                        cm = real(sum(real(x(:,dim),8))*mpartcdm+sum(real(y(:,dim),8))*mpartgas,4)/mtot
                        x(:,dim) = x(:,dim)-cm
                        y(:,dim) = y(:,dim)-cm
                        cm = real(sum(real(vx(:,dim),8))*mpartcdm+sum(real(vy(:,dim),8))*mpartgas,4)/mtot
                        vx(:,dim) = vx(:,dim)-cm
                        vy(:,dim) = vy(:,dim)-cm
                     end do
                     
                     ! convert comoving to physical units
                     x = x*scale/h*1e3    ! [pkpc]
                     y = y*scale/h*1e3    ! [pkpc]
                     vx = vx*sqrt(scale)  ! [pkm/s] see Gadget-2 manual
                     vy = vy*sqrt(scale)  ! [pkm/s]
                     
                     ! compute angular momenta
                     jcdm(1) = sum(x(:,2)*vx(:,3)-x(:,3)*vx(:,2))/ncdm ! [pkpc km/s]
                     jcdm(2) = sum(x(:,3)*vx(:,1)-x(:,1)*vx(:,3))/ncdm ! [pkpc km/s]
                     jcdm(3) = sum(x(:,1)*vx(:,2)-x(:,2)*vx(:,1))/ncdm ! [pkpc km/s]
                     jgas(1) = sum(y(:,2)*vy(:,3)-y(:,3)*vy(:,2))/ngas ! [pkpc km/s]
                     jgas(2) = sum(y(:,3)*vy(:,1)-y(:,1)*vy(:,3))/ngas ! [pkpc km/s]
                     jgas(3) = sum(y(:,1)*vy(:,2)-y(:,2)*vy(:,1))/ngas ! [pkpc km/s]
                     
                     ! compute mass-weighted rms radius [pkpc]
                     radius_rms = sqrt((mcdm/ncdm*sum(x(:,1)**2+x(:,2)**2+x(:,3)**2)+ &
                                      & mgas/ngas*sum(y(:,1)**2+y(:,2)**2+y(:,3)**2))/mtot)
                     
                     ! compute approximate spin parameter (following eq. (4) of Obreschkow 2015b)
                     jmean = (mcdm*sqrt(sum(jcdm**2))+mgas*sqrt(sum(jgas**2)))/mtot ! [pkpc km/s]
                     z = 1/scale-1 ! redshift
                     E = sqrt(OmegaM*(1+z)**3+OmegaL)
                     Hubble = 100*1e3/unit%Mpc*h*E ! [Hz] Hubble constant at z
                     Delta = 18*const%pi**2+82*OmegaL/E**2-39*OmegaL**2/E**4 ! virial overdensity at z
                     lambda_factor = Hubble**(1.0/3.0)*Delta**(1.0/6.0)/(2*const%G)**(2.0/3.0) ! [SI]
                     lambda = real(lambda_factor*(jmean*unit%kpc*1e3)/(store(k)%mass*real(unit%Msun,8)/h)**(2.0/3.0),4)
   
                     ! store values
                     if (mode==1) then
                     
                        store(k)%radius_rms_final = radius_rms
                        store(k)%lambda_final = lambda
                        store(k)%jcdm_final = sqrt(sum(jcdm**2))
                        store(k)%jgas_final = sqrt(sum(jgas**2))
                        
                     else
                     
                        store(k)%radius_rms_initial = radius_rms
                        store(k)%lambda_initial = lambda
                        store(k)%jcdm_initial = sqrt(sum(jcdm**2))
                        store(k)%jgas_initial = sqrt(sum(jgas**2))
                        
                     end if
                     
                  end do
                     
                  ! compute asymmetry parameter
                  call get_dimensions(x,y,vx,vy,a) ! size of initial halo in its proper coordinates
                  collapse_factor = store(k)%radius_rms_initial/store(k)%radius_rms_final
                  store(k)%ax = a(1)
                  store(k)%ay = a(2)
                  store(k)%az = a(3)
                  store(k)%asymmetry = abs(a(1)**2-a(2)**2)/a(3)**2*max(1.0,collapse_factor)*store(k)%lambda_initial
                                    
               end if
               
            end if
      
         end if
         
      end if
      
   end do
   
   ! save data
   write(snstr,'(I0.3)') initial_snapshot
   if (subhalos) then
      fn = trim(para%path_analysis)//trim(snfile(para%snapshot))//'_'//trim(snstr)//'_asymmetry_subhalos1'
   else
      fn = trim(para%path_analysis)//trim(snfile(para%snapshot))//'_'//trim(snstr)//'_asymmetry_subhalos0'
   end if
   call hdf5_create(trim(fn)//'.hdf') ! create HDF5 file
   call hdf5_open(trim(fn)//'.hdf',.true.) ! open HDF5 file
      
   ! Group "simulation"
   call hdf5_add_group('simulation')
   call hdf5_write_data('simulation/name',trim(para%simulation),'simulation name')
   call hdf5_write_data('simulation/box_l',para%L,'[simulation units] box side length')
   call hdf5_write_data('simulation/box_n',para%N,'cubic root of particle number')
   call hdf5_write_data('simulation/snapshot_final',para%snapshot,'index of final snapshot')
   call hdf5_write_data('simulation/snapshot_initial',initial_snapshot,'index of initial snapshot')
   
   ! Group "surfanalysis"
   call hdf5_add_group('surfanalysis')
   call hdf5_write_data('surfanalysis/timestamp',timestamp(),'surfanalysis timestamp')
   call hdf5_write_data('surfanalysis/version',trim(version),'version of surfanalysis used to extract the halo')
   call hdf5_write_data('surfanalysis/developer','Danail Obreschkow; danail.obreschkow@icrar.org')
   call hdf5_write_data('surfanalysis/subhalos',log2int(subhalos), &
   & 'logical flag (0/1 = subhalo particles are included/excluded)')
   
   ! Group "halos"
   call hdf5_add_group('halos')
   call hdf5_write_data('halos/id',store(1:k)%id,'unique halo id')
   call hdf5_write_data('halos/ntot',store(1:k)%ntot,'total number of particles')
   call hdf5_write_data('halos/ngas',store(1:k)%ngas,'number of gas particles')
   call hdf5_write_data('halos/ncdm',store(1:k)%ncdm,'number of dark matter particles')
   call hdf5_write_data('halos/mass',store(1:k)%mass,'[Msun/h] total mass of this group')
   call hdf5_write_data('halos/fbaryon',store(1:k)%fbaryon,'baryon mass fraction')
   call hdf5_write_data('halos/radius_rms_initial',store(1:k)%radius_rms_initial, &
   & '[pkpc] mass-weighted rms radius of halo at initial snapshot')
   call hdf5_write_data('halos/lambda_initial',store(1:k)%lambda_initial, &
   & 'spin parameter at initial snapshot')
   call hdf5_write_data('halos/jcdm_initial',store(1:k)%jcdm_initial, &
   & '[pkpc km/s] norm of specific AM of CDM at initial snapshot')
   call hdf5_write_data('halos/jgas_initial',store(1:k)%jgas_initial, &
   & '[pkpc km/s] norm of specific AM of gas at initial snapshot')
   call hdf5_write_data('halos/radius_rms_final',store(1:k)%radius_rms_final, &
   & '[pkpc] mass-weighted rms radius of halo at final snapshot')
   call hdf5_write_data('halos/lambda_final',store(1:k)%lambda_final, &
   & 'spin parameter at final snapshot')
   call hdf5_write_data('halos/jcdm_final',store(1:k)%jcdm_final, &
   & '[pkpc km/s] norm of specific AM of CDM at final snapshot')
   call hdf5_write_data('halos/jgas_final',store(1:k)%jgas_final, &
   & '[pkpc km/s] norm of specific AM of gas at final snapshot')
   call hdf5_write_data('halos/asymmetry',store(1:k)%asymmetry, &
   & 'dimensionless asymmetry parameter at initial snapshot')
   call hdf5_write_data('halos/ax',store(1:k)%ax, &
   & '[pkpc] mass-weighted rms length along major axis perpendicular to spin for initial snapshot')
   call hdf5_write_data('halos/ay',store(1:k)%ay, &
   & '[pkpc] mass-weighted rms length along minor axis perpendicular to spin for initial snapshot')
   call hdf5_write_data('halos/az',store(1:k)%az, &
   & '[pkpc] mass-weighted rms length the spin axis for initial snapshot')
   
   call hdf5_close() ! close HDF5 file
   
end subroutine compute_asymmetry_parameters

subroutine get_dimensions(x_,y_,vx_,vy_,a)

   ! returns the characteristic sizes of the halo along its three proper axes, defined such that
   ! the z-axis is aligned with the AM vector, and the x and y axes are quadupole eigenvectors in the
   ! plane orthogonal to the AM vector

   implicit none
   real*4,allocatable,intent(in)    :: x_(:,:),y_(:,:)   ! center of mass positions
   real*4,allocatable,intent(in)    :: vx_(:,:),vy_(:,:) ! center of mass velocities
   real*4,intent(out)               :: a(3)              ! size of the halo along its three proper axes
   real*4,allocatable               :: x(:,:),y(:,:)     ! rotated center of mass positions
   real*4,allocatable               :: vx(:,:),vy(:,:)   ! rotated center of mass velocities
   integer*4                        :: n,nx,ny,i,dim
   real*4                           :: j(3)
   real*4                           :: R(3,3)
   real*4,parameter                 :: ex(3) = (/1.0,0.0,0.0/)
   real*4,parameter                 :: ez(3) = (/0.0,0.0,1.0/)
   real*4                           :: axis(3),angle
   real*4                           :: Q(2,2),values(2),vectors(2,2)
   
   ! number of particles
   nx = size(x_,1)
   ny = size(y_,1)
   n = nx+ny
   
   ! initialise variable positions and velocities
   allocate(x(nx,3),y(ny,3),vx(nx,3),vy(ny,3))
   x = x_
   y = y_
   vx = vx_
   vy = vy_
   
   ! rotate such that AM points along the z-axis
   j(1) = sum(x(:,2)*vx(:,3)-x(:,3)*vx(:,2))*(1-fb)+sum(y(:,2)*vy(:,3)-y(:,3)*vy(:,2))*fb
   j(2) = sum(x(:,3)*vx(:,1)-x(:,1)*vx(:,3))*(1-fb)+sum(y(:,3)*vy(:,1)-y(:,1)*vy(:,3))*fb
   j(3) = sum(x(:,1)*vx(:,2)-x(:,2)*vx(:,1))*(1-fb)+sum(y(:,1)*vy(:,2)-y(:,2)*vy(:,1))*fb
   axis = cross_product(j,ez)
   angle = acos(j(3)/sqrt(sum(j**2)))
   R = rotation_matrix(axis,angle)
   do i=1,nx
      x(i,:) = matmul(R,x(i,:))
      vx(i,:) = matmul(R,vx(i,:))
   end do
   do i=1,ny
      y(i,:) = matmul(R,y(i,:))
      vy(i,:) = matmul(R,vy(i,:))
   end do
   
   ! inertia tensor in the proper x-y-plane
   Q(1,1) = (1-fb)*sum(x(:,2)**2)+fb*sum(y(:,2)**2)
   Q(2,2) = (1-fb)*sum(x(:,1)**2)+fb*sum(y(:,1)**2)
   Q(1,2) = -(1-fb)*sum(x(:,1)*x(:,2))-fb*sum(y(:,1)*y(:,2))
   Q(2,1) = Q(1,2)
   call eigen(Q,values,vectors) ! vectors(:,1) and vectors(:,2) are the eigenvectors
   
   ! rotate about z-axis to align x and y with the principal axes of the quadrupole
   angle = -atan2(vectors(2,1),vectors(1,1))
   R = rotation_matrix(ez,angle)
   do i=1,nx
      x(i,:) = matmul(R,x(i,:))
      vx(i,:) = matmul(R,vx(i,:))
   end do
   do i=1,ny
      y(i,:) = matmul(R,y(i,:))
      vy(i,:) = matmul(R,vy(i,:))
   end do
   
   ! alignment checks
   ! j(1) = sum(x(:,2)*vx(:,3)-x(:,3)*vx(:,2))*(1-fb)+sum(y(:,2)*vy(:,3)-y(:,3)*vy(:,2))*fb
   ! j(2) = sum(x(:,3)*vx(:,1)-x(:,1)*vx(:,3))*(1-fb)+sum(y(:,3)*vy(:,1)-y(:,1)*vy(:,3))*fb
   ! j(3) = sum(x(:,1)*vx(:,2)-x(:,2)*vx(:,1))*(1-fb)+sum(y(:,1)*vy(:,2)-y(:,2)*vy(:,1))*fb
   ! write(*,*) 'check 1:',j/sqrt(sum(j**2)) ! should be 0,0,1
   ! Q(1,1) = (1-fb)*sum((2*x(:,1)**2-x(:,2)**2-x(:,3)**2))+fb*sum((2*y(:,1)**2-y(:,2)**2-y(:,3)**2))
   ! Q(2,2) = (1-fb)*sum((2*x(:,2)**2-x(:,3)**2-x(:,1)**2))+fb*sum((2*y(:,2)**2-y(:,3)**2-y(:,1)**2))
   ! Q(1,2) = (1-fb)*sum((3*x(:,1)*x(:,2)))+fb*sum((3*y(:,1)*y(:,2)))
   ! Q(2,1) = Q(1,2)
   ! call eigen(Q,values,vectors) ! vectors(:,1) and vectors(:,2) are the eigenvectors
   ! write(*,*) 'check 2:',vectors(:,1) ! should be 1,0
   
   ! size vector
   do dim = 1,3
      a(dim) = sqrt(((1-fb)*sum(x(:,dim)**2)+fb*sum(y(:,dim)**2))/((1-fb)*nx+fb*ny))
   end do
   
   !if (a(2)>a(1)) then
   !   write(*,*) a
   !   open(1,file='/Users/00077389/Desktop/test.bin',action='write',form='unformatted',status='replace',access='stream')
   !   write(1) 1.0*nx,1.0*ny,x_,y_,vx_,vy_
   !   close(1)
   !   stop
   !end if
    
end subroutine get_dimensions

real*4 function median(vector)

   implicit none
   real*4,intent(in)       :: vector(:)
   real*4,allocatable      :: x(:)
   integer*4               :: n
   integer*4,allocatable   :: index(:)
   
   n = size(vector)
   allocate(x(n),index(n))
   x = vector
   index = 0
   
   call sort(x,index)
   
   median = x((n+1)/2)
   
end function median
   
end module module_asymmetries