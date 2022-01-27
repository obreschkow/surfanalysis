module module_showhalo

   use shared_module_core
   use shared_module_arguments
   use shared_module_graphics
   use module_global
   use module_io
   use module_gethalo
   use module_trackhalo
   use module_processhalo

   implicit none
   
   private
   public   :: task_showhalo
   public   :: load_halo_particles_at_different_snapshot

   type type_display_options
      integer*4   :: mode
      integer*4   :: npixels
      real*4      :: sidelength ! [simulation units]
      real*4      :: smoothinglength  ! [simulation units]
      real*4      :: gamma
      real*4      :: lum
      integer*4   :: projection
   end type type_display_options

   type(type_display_options) :: graph
   real*4,allocatable         :: rgb(:,:,:)

contains

subroutine task_showhalo

   implicit none
   
   integer*4                  :: haloid
   character(len=255)         :: outputfile
   logical                    :: subhalos
   integer*4                  :: at
   
   call get_task_value(haloid)
   
   ! handle options
   call get_option_value(outputfile,'-outputfile','')
   call get_option_value(subhalos,'-subhalos',.false.)
   call get_option_value(at,'-at',para%snapshot,min=0)
   call get_option_value(graph%mode,'-mode',0,min=0,max=2)
   call get_option_value(graph%npixels,'-npixels',800,min=1,max=4000) 
   call get_option_value(graph%sidelength,'-sidelength',2.0,min=0.0) 
   call get_option_value(graph%smoothinglength,'-smoothinglength',0.02,min=0.0) 
   call get_option_value(graph%lum,'-lum',1.0,min=0.0) 
   call get_option_value(graph%gamma,'-gamma',0.6,min=0.0) 
   call get_option_value(graph%projection,'-projection',1,min=1,max=3) 
   call require_no_options_left
   
   ! additional checks
   if ((graph%mode==1).and.(graph%smoothinglength>=graph%sidelength)) call error('smoothinglength must be smaller than sidelength')
   
   call load_halo_particles(haloid,subhalos)
   
   if (at==para%snapshot) then
      call load_halo_particles(haloid,subhalos)
   else
      call load_halo_particles_at_different_snapshot(haloid,subhalos,at)
   end if
   
   call center_particles
   call center_velocities
   
   call out('Render halo ',haloid)
   call raster_halo
   
   if (isempty(outputfile)) then
      call raster2bitmap(rgb,'.tmp.bmp')
      call system('open .tmp.bmp')
   else
      call raster2bitmap(rgb,outputfile)
   end if

end subroutine task_showhalo

subroutine load_halo_particles_at_different_snapshot(haloid,include_subhalos,at,supressverbose)

   ! this calls load_halo_evolving_particles for a single snapshot "at"

   implicit none
   integer*4,intent(in)          :: haloid
   logical,intent(in)            :: include_subhalos
   integer*4,intent(in)          :: at
   logical,optional,intent(in)   :: supressverbose
   integer*4                     :: d
   real*4,allocatable            :: x(:,:,:)
   real*4,allocatable            :: v(:,:,:)

   call load_halo_evolving_particles(haloid,include_subhalos,.false.,at,at,x,v,supressverbose)
   do d = 1,3
      p(:)%x(d) = x(:,d,at)
      p(:)%v(d) = v(:,d,at)
   end do

end subroutine load_halo_particles_at_different_snapshot

subroutine raster_halo

   implicit none
   real*4,allocatable            :: f(:,:,:),x(:),y(:),vx(:),vy(:)
   real*4,allocatable            :: kernel(:,:,:)
   integer*4                     :: nsmooth,npx
   integer*4                     :: kernelsteps = 50
   integer*4                     :: j,k,ix,iy,n
   real*4                        :: q0,q,factor,x0,y0,dx,dy,dt,density,add
   real*8                        :: rrms,vrms
   integer*8                     :: i,length
   
   ! rotation
   allocate(x(nparticles),y(nparticles),vx(nparticles),vy(nparticles))
   if (graph%projection==1) then
      x = p%x(1)
      y = p%x(2)
      vx = p%v(1)
      vy = p%v(2)
   else if (graph%projection==2) then
      x = p%x(2)
      y = p%x(3)
      vx = p%v(2)
      vy = p%v(3)
   else if (graph%projection==3) then
      x = p%x(3)
      y = p%x(1)
      vx = p%v(3)
      vy = p%v(1)
   end if
   
   ! determine dt
   rrms = sqrt(sum(real(x**2+y**2,8))/nparticles)
   vrms = sqrt(sum(real(vx**2+vy**2,8))/nparticles)
   dt = real(rrms/vrms*0.3,4)
   density = real(nparticles/(rrms/graph%sidelength*graph%npixels)**2,4)
   npx = max(1,nint(dt*vrms/graph%sidelength*graph%npixels))
   
   ! raster channels
   if (graph%mode == 0) then
   
      allocate(f(1:graph%npixels,1:graph%npixels,2))
      f = 0
   
      do i = 1,nparticles
         ix = nint((x(i)/graph%sidelength+0.5)*graph%npixels)
         if ((ix>=1).and.(ix<=graph%npixels)) then
            iy = nint((y(i)/graph%sidelength+0.5)*graph%npixels)
            if ((iy>=1).and.(iy<=graph%npixels)) then
               f(ix,iy,p(i)%species) = f(ix,iy,p(i)%species)+1
            end if
         end if
      end do
   
   else if (graph%mode == 1) then
   
      ! make smoothing kernel
      nsmooth = int(graph%npixels*graph%smoothinglength/graph%sidelength*2)
      allocate(kernel(-nsmooth:nsmooth,-nsmooth:nsmooth,kernelsteps))
      do i = -nsmooth,nsmooth
         do j = -nsmooth,nsmooth
            q0 = sqrt(real(i**2+j**2))/nsmooth*2
            do k = 1,kernelsteps
               factor = 1+0.015*(k-1)**2.0
               q = q0*factor
               if (q<1) then
                  kernel(i,j,k) = (1-1.5*q**2+0.75*q**3)
               else if (q<2) then
                  kernel(i,j,k) = (0.25*(2-q)**3)
               else
                  kernel(i,j,k) = 0
               end if
            end do
         end do
      end do
      do k = 1,kernelsteps
         kernel(:,:,k) = kernel(:,:,k)/sum(kernel(:,:,k))
      end do
      
      allocate(f(1-2*nsmooth:graph%npixels+2*nsmooth,1-2*nsmooth:graph%npixels+2*nsmooth,2))
      f = 0
   
      ! smooth all particles
      do i = 1,nparticles
         ix = nint((x(i)/graph%sidelength+0.5)*graph%npixels)
         if ((ix>=1-nsmooth).and.(ix<=graph%npixels+nsmooth)) then
            iy = nint((y(i)/graph%sidelength+0.5)*graph%npixels)
            if ((iy>=1-nsmooth).and.(iy<=graph%npixels+nsmooth)) then
               k = min(nint(f(ix,iy,p(i)%species)/density**1*100+1),kernelsteps)
               f(ix-nsmooth:ix+nsmooth,iy-nsmooth:iy+nsmooth,p(i)%species) = &
               &  f(ix-nsmooth:ix+nsmooth,iy-nsmooth:iy+nsmooth,p(i)%species)+kernel(:,:,k)
            end if
         end if
      end do
      
   else if (graph%mode == 2) then
   
      allocate(f(1:graph%npixels,1:graph%npixels,2))
      f = 0
   
      length = 0
      do i = 1,nparticles
         x0 = (x(i)/graph%sidelength+0.5)*graph%npixels
         y0 = (y(i)/graph%sidelength+0.5)*graph%npixels
         dx = vx(i)*dt/graph%sidelength*graph%npixels
         dy = vy(i)*dt/graph%sidelength*graph%npixels
         n = nint(sqrt(dx**2+dy**2))+1
         add = 1/real(n,4)
         length = length+n
         do k = -n,n
            ix = nint(x0+real(k)/n/2*dx)
            if ((ix>=1).and.(ix<=graph%npixels)) then
               iy = nint(y0+real(k)/n/2*dy)
               if ((iy>=1).and.(iy<=graph%npixels)) f(ix,iy,p(i)%species) = f(ix,iy,p(i)%species)+add
            end if
         end do
      end do
      f = f/real(npx)*real(length)/real(nparticles)
   
   end if
   
   ! contrast
   f = f/density*2
   f = f**graph%gamma*graph%lum
   f = 1.0-exp(-f)
   if ((graph%mode == 0).or.(graph%mode == 2)) where(f>0) f=0.1+0.9*f
   
   ! convert to rgb
   allocate(rgb(graph%npixels,graph%npixels,3))
   rgb = 0
   rgb(:,:,1) = f(1:graph%npixels,1:graph%npixels,1)+0.2*f(1:graph%npixels,1:graph%npixels,2)
   rgb(:,:,2) = 0.2*f(1:graph%npixels,1:graph%npixels,2)
   rgb(:,:,3) = f(1:graph%npixels,1:graph%npixels,2)
   where (rgb>1) rgb = 1
   
end subroutine raster_halo

end module module_showhalo