module module_processhalo

use shared_module_core
use shared_module_arguments
use module_global

contains

subroutine center_particles(shift1,shift2)

   implicit none
   real*4,optional,intent(out)   :: shift1(3),shift2(3)
   integer*4                     :: d
   real*4                        :: h,c
   
   if (present(shift1)) shift1=0
   if (present(shift2)) shift2=0
   
   do d = 1,3
   
      ! wrap coordinates
      h = para%L/2
      if (maxval(p%x(d))-minval(p%x(d))>h) then
         p%x(d) = mod(p%x(d)+h,para%L)
         if (maxval(p%x(d))-minval(p%x(d))>h) call error('halo larger than half the simulation box.')
         if (present(shift1)) shift1(d)=h
      end if

      ! reset center of mass
      c = real(sum(real(p%x(d),8))/nparticles,4)
      p%x(d) = p%x(d)-c
      if (present(shift2)) shift2(d)=-c
   
   end do

end subroutine center_particles

subroutine center_velocities

   implicit none
   integer*4                     :: d
   real*4                        :: c
   
   do d = 1,3
   
      ! reset center of mass velocity
      c = real(sum(real(p%v(d),8))/nparticles,4)
      p%v(d) = p%v(d)-c
   
   end do

end subroutine center_velocities

end module module_processhalo