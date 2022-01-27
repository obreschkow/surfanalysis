module module_getparticle

use shared_module_core
use shared_module_arguments
use module_global
use module_io

implicit none

contains

subroutine task_getparticle

   implicit none
   integer*8            :: particleid
   
   call get_task_value(particleid)
   call write_particle(particleid)
   
end subroutine task_getparticle

subroutine write_particle(id)

   implicit none
   integer*8,intent(in) :: id
   character(len=255)   :: txt
   type(type_particle)  :: particle
   integer*4            :: ifile
   
   call get_particle(id,particle,ifile)
   call hline
   call out('Simulation:  ',trim(para%simulation))
   call out('Particle:    ',particle%id)
   call out('Species:     ',particle%species)
   call out('File #:      ',ifile)
   write(txt,'(A,F11.4,A,F11.4,A,F11.4)') 'Position: ',particle%x(1),',',particle%x(2),',',particle%x(3)
   call out(txt)
   write(txt,'(A,F11.4,A,F11.4,A,F11.4)') 'Velocity: ',particle%v(1),',',particle%v(2),',',particle%v(3)
   call out(txt)
   call hline
   
end subroutine write_particle

subroutine get_particle(id,particle,ifile)

   implicit none
   integer*8,intent(in)             :: id
   type(type_particle),intent(out)  :: particle
   integer*4,intent(out),optional   :: ifile
   integer*8                        :: position
   integer*8                        :: file_size
   character(255)                   :: fn
   
   if (id<1) call error('Particle ID below allowed range.')
   
   ifile = int((id-1)/nparticles_per_sorted_file,4)
   position = (id-ifile*nparticles_per_sorted_file-1)*bytes_per_particle+1_8
   
   fn = trim(filename(ifile,para%path_surfsuite,snfile(para%snapshot),para%ext_sorted))
   
   if (.not.exists(fn)) then
      fn = trim(filename(0,para%path_surfsuite,snfile(para%snapshot),para%ext_sorted))
      if (.not.exists(fn)) then
         write(*,*) fn
         call error('Sorted particle files do not exist for this simulation. Consider running the task "sortparticles".')
      else
         call error('Particle ID outside allowed range or file lost.')
      end if
   end if
   
   inquire(file=trim(fn), size=file_size)
   
   if (position>file_size) call error('Particle ID outside allowed range.')
   
   open(1,file=trim(fn), &
   & action='read',form='unformatted',status='old',access='stream')
   read(1,pos=position) particle
   close(1)

end subroutine get_particle
   
end module module_getparticle