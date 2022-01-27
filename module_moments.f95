module module_moments

use shared_module_core
use shared_module_arguments
use shared_module_hdf5
use shared_module_lapack
use shared_module_maths
use module_global
use module_io
use module_gethalo
use module_processhalo

implicit none

private
public   :: task_moments

integer*4,parameter  :: lmax = 4

type type_properties
   integer*8   :: id
   real*4      :: mass
   real*4      :: darkness
   real*4      :: mu(0:lmax)
   real*4      :: nu(0:lmax)
   real*4      :: jcdm(3),jcdm_norm
   real*4      :: jgas(3),jgas_norm
   real*4      :: inertia(3)
end type type_properties

real*4,parameter     :: fraction = 1 ! fraction of halos to analyze
real*4,parameter     :: mass_min = 1e11 ! [Msun/h]
real*4,parameter     :: mass_max = 1e16 ! [Msun/h]
real*4,parameter     :: h = 0.6751
real*4,parameter     :: OmegaM = 0.3121
real*4,parameter     :: OmegaL = 1-OmegaM
real*4,parameter     :: fb = 0.1573 ! baryon fraction relative to all matter
   
real*4               :: mpart ! [Msun/h] particle mass
logical              :: subhalos

contains

subroutine task_moments

   implicit none
   
   ! read options
   call get_option_value(subhalos,'-subhalos',.true.)
   call require_no_options_left
   
   mpart = para%L**3/real(para%N,4)**3*OmegaM*2.774495e11 ! [Msun/h] particle mass (in the case of hydro runs, this is the mass of 1 CDM+1 gas particle)
   
   call tic()
   call make_moment_indices
   call toc()
   
end subroutine task_moments

subroutine make_moment_indices
   
   implicit none
   
   integer*4                           :: i,j,ix,iy,ntot,ntyp(2),nmin,nmax
   character(len=255)                  :: fn
   type(type_halo),allocatable         :: halo(:)
   real*4                              :: mtyp(2),mpartgas,mpartCDM
   type(type_properties),allocatable   :: store(:)
   real*4,allocatable                  :: x(:,:),y(:,:),vx(:,:),vy(:,:)
   integer*4                           :: k,nhalos
   real*4                              :: darkness
   real*4                              :: mu(0:lmax)
   real*4                              :: nu(0:lmax)
   real*4                              :: t(3,3) ! tensor of intertia
   integer*8                           :: npart,ipart
   real*4                              :: f_complete,f_complete_next
   
   mpartgas = fb*mpart     ! [Msun/h]
   mpartCDM = (1-fb)*mpart ! [Msun/h]
   nmin = floor(mass_min/mpartCDM)
   nmax = ceiling(mass_max/mpartgas)
   
   call load_all_halo_properties(halo)
   nhalos = size(halo)
   nhalos = int(nhalos*fraction,4)
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
            
            ntyp(1) = count(p%species==1)
            ntyp(2) = count(p%species==2)
            mtyp(1) = mpartgas*ntyp(1) ! [Msun/h]
            mtyp(2) = mpartCDM*ntyp(2) ! [Msun/h]
            
            if ((sum(mtyp)>=mass_min).and.(sum(mtyp)<=mass_max)) then
            
               if ((ntyp(1)>=20).and.(ntyp(2)>=20)) then
               
                  k = k+1
               
                  ! save particle id
                  store(k)%id = i
                  
                  ! compute and save mass of group
                  store(k)%mass = sum(mtyp)        ! [Msun/h]
                  
                  ! center particles
                  call center_particles ! necessary, even if particles are centered again later, because of wrapping
                  
                  ! separate CDM and gas particles
                  if (allocated(x)) deallocate(x); allocate(x(ntyp(2),3))
                  if (allocated(y)) deallocate(y); allocate(y(ntyp(1),3))
                  if (allocated(vx)) deallocate(vx); allocate(vx(ntyp(2),3))
                  if (allocated(vy)) deallocate(vy); allocate(vy(ntyp(1),3))
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
      
                  ! compute and save moment indices
                  call shape(x,y,darkness,mu,nu)
                  store(k)%darkness = darkness
                  store(k)%mu = mu
                  store(k)%nu = nu
                  
                  ! compute angular momenta
                  store(k)%jcdm(1) = sum(x(:,2)*vx(:,3)-x(:,3)*vx(:,2))
                  store(k)%jcdm(2) = sum(x(:,3)*vx(:,1)-x(:,1)*vx(:,3))
                  store(k)%jcdm(3) = sum(x(:,1)*vx(:,2)-x(:,2)*vx(:,1))
                  store(k)%jcdm_norm = sqrt(store(k)%jcdm(1)**2+store(k)%jcdm(2)**2+store(k)%jcdm(3)**2)
                  store(k)%jgas(1) = sum(y(:,2)*vy(:,3)-y(:,3)*vy(:,2))
                  store(k)%jgas(2) = sum(y(:,3)*vy(:,1)-y(:,1)*vy(:,3))
                  store(k)%jgas(3) = sum(y(:,1)*vy(:,2)-y(:,2)*vy(:,1))
                  store(k)%jgas_norm = sqrt(store(k)%jgas(1)**2+store(k)%jgas(2)**2+store(k)%jgas(3)**2)
                  
                  ! compute tensor of inertia
                  t(1,1) = mpartCDM*sum(x(:,2)**2+x(:,3)**2)+mpartgas*sum(y(:,2)**2+y(:,3)**2)
                  t(2,2) = mpartCDM*sum(x(:,1)**2+x(:,3)**2)+mpartgas*sum(y(:,1)**2+y(:,3)**2)
                  t(3,3) = mpartCDM*sum(x(:,1)**2+x(:,2)**2)+mpartgas*sum(y(:,1)**2+y(:,2)**2)
                  t(1,2) = -mpartCDM*sum(x(:,1)*x(:,2))-mpartgas*sum(y(:,1)*y(:,2))
                  t(1,3) = -mpartCDM*sum(x(:,1)*x(:,3))-mpartgas*sum(y(:,1)*y(:,3))
                  t(2,3) = -mpartCDM*sum(x(:,2)*x(:,3))-mpartgas*sum(y(:,2)*y(:,3))
                  t(2,1) = t(1,2)
                  t(3,1) = t(1,3)
                  t(3,2) = t(2,3)
                  call eigen(t,store(k)%inertia)
                  
               end if
               
            end if
      
         end if
         
      end if
      
   end do
   
   
   
   ! save data
   fn = trim(para%path_analysis)//trim(snfile(para%snapshot))//'_momentindex'
   call hdf5_create(trim(fn)//'.hdf') ! create HDF5 file
   call hdf5_open(trim(fn)//'.hdf',.true.) ! open HDF5 file
      
   ! Group "simulation"
   call hdf5_add_group('simulation')
   call hdf5_write_data('simulation/name',trim(para%simulation),'simulation name')
   call hdf5_write_data('simulation/box_l',para%L,'[simulation units] box side length')
   call hdf5_write_data('simulation/box_n',para%N,'cubic root of particle number')
   call hdf5_write_data('simulation/snapshot',para%snapshot,'snapshot index')
   
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
   call hdf5_write_data('halos/mass',store(1:k)%mass,'[simulation units] total mass of this group')
   call hdf5_write_data('halos/darkness',store(1:k)%darkness,'species index mu(-1)')
   call hdf5_write_data('halos/mu0',store(1:k)%mu(0), &
   & 'moment index 0 (monopole difference between CDM and baryons)')
   call hdf5_write_data('halos/mu1',store(1:k)%mu(1), &
   & 'moment index 1 (dipole of difference between CDM and baryons)')
   call hdf5_write_data('halos/mu2',store(1:k)%mu(2), &
   & 'moment index 2 (quadrupole difference between CDM and baryons)')
   call hdf5_write_data('halos/mu3',store(1:k)%mu(3), &
   & 'moment index 3 (octupole difference between CDM and baryons)')
   call hdf5_write_data('halos/mu4',store(1:k)%mu(4), &
   & 'moment index 4 (hexadecapole difference between CDM and baryons)')
   call hdf5_write_data('halos/nu0',store(1:k)%nu(0), &
   & 'moment index 0 (monopole of CDM-to-baryons contrast)')
   call hdf5_write_data('halos/nu1',store(1:k)%nu(1), &
   & 'moment index 1 (dipole of CDM-to-baryons contrast)')
   call hdf5_write_data('halos/nu2',store(1:k)%nu(2), &
   & 'moment index 2 (quadrupole of CDM-to-baryons contrast)')
   call hdf5_write_data('halos/nu3',store(1:k)%nu(3), &
   & 'moment index 3 (octupole of CDM-to-baryons contrast)')
   call hdf5_write_data('halos/nu4',store(1:k)%nu(4), &
   & 'moment index 4 (hexadecapole of CDM-to-baryons contrast)')
   call hdf5_write_data('halos/jcdmx',store(1:k)%jcdm(1),'x-component of specific AM of CDM')
   call hdf5_write_data('halos/jcdmy',store(1:k)%jcdm(2),'y-component of specific AM of CDM')
   call hdf5_write_data('halos/jcdmz',store(1:k)%jcdm(3),'z-component of specific AM of CDM')
   call hdf5_write_data('halos/jcdm',store(1:k)%jcdm_norm,'norm of specific AM of CDM')
   call hdf5_write_data('halos/jgasx',store(1:k)%jgas(1),'x-component of specific AM of gas')
   call hdf5_write_data('halos/jgasy',store(1:k)%jgas(2),'y-component of specific AM of gas')
   call hdf5_write_data('halos/jgasz',store(1:k)%jgas(3),'z-component of specific AM of gas')
   call hdf5_write_data('halos/jgas',store(1:k)%jgas_norm,'norm of specific AM of gas')
   call hdf5_write_data('halos/inertia_min',store(1:k)%inertia(1),'smallest eigenvalue of the moment of inertia tensor')
   call hdf5_write_data('halos/inertia_mid',store(1:k)%inertia(2),'intermediate eigenvalue of the moment of inertia tensor')
   call hdf5_write_data('halos/inertia_max',store(1:k)%inertia(3),'largest eigenvalue of the moment of inertia tensor')

   call hdf5_close() ! close HDF5 file
   
end subroutine make_moment_indices

subroutine shape(x_,y_,d,mu,nu)

   implicit none
   real*4,allocatable,intent(in)    :: x_(:,:),y_(:,:)
   real*4,allocatable               :: x(:,:),y(:,:) ! coordinates relative to geometric centre
   real*4,intent(out)               :: d,mu(0:lmax),nu(0:lmax) ! moments of the difference, difference of the moments
   integer*4                        :: n,nx,ny
   real*4                           :: rxmean,rymean,cg
   real*4,allocatable               :: rx(:),ry(:)
   complex*8,allocatable            :: fx(:),fy(:)
   real*4                           :: kx,ky,dk,prefactor
   integer*4                        :: dim,l,m,i ! loop counters
   
   ! initialize moment array
   mu = 0
   nu = 0
   
   ! species index
   nx = size(x_,1)
   ny = size(y_,1)
   n = nx+ny
   d = (nx-ny)/real(n)
   
   ! center particles to mid-point between CM of mass of DM and gas
   allocate(x(nx,3),y(ny,3))
   do dim = 1,3
      cg = real(sum(real(x_(:,dim),8))/real(nx,8)+sum(real(y_(:,dim),8))/real(ny,8),4)*0.5
      x(:,dim) = x_(:,dim)-cg
      y(:,dim) = y_(:,dim)-cg
   end do
  
   ! radius vector
   allocate(rx(nx),ry(ny))
   rx = sqrt(x(:,1)**2+x(:,2)**2+x(:,3)**2)
   ry = sqrt(y(:,1)**2+y(:,2)**2+y(:,3)**2)
   
   ! moment analysis (core part of this function)
   rxmean = sum(rx)/real(nx)
   rymean = sum(ry)/real(ny)
   
   do l = 0,lmax
      if (allocated(fx)) deallocate(fx); allocate(fx(-l:l)); fx=0
      if (allocated(fy)) deallocate(fy); allocate(fy(-l:l)); fy=0
      do m = -l,l
         do i = 1,nx
            fx(m) = fx(m)+spherical_harmonic(l,m,x(i,:))*rx(i)
         end do
         do i = 1,ny
            fy(m) = fy(m)+spherical_harmonic(l,m,y(i,:))*ry(i)
         end do
      end do
      
      fx = fx/real(nx)
      fy = fy/real(ny)
      
      prefactor = sqrt(4*pi/(2.0*l+1.0))
      kx = sqrt(sum(abs(fx)**2))
      ky = sqrt(sum(abs(fy)**2))
      dk = sqrt(sum(abs(fx-fy)**2))
      
      mu(l) = prefactor*(kx-ky)/max(rxmean,rymean)
      nu(l) = prefactor*dk/(rxmean+rymean)
            
    end do
    
end subroutine shape
   
end module module_moments