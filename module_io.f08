module module_io

use shared_module_core
use shared_module_parameters
use module_global

contains

function snfile(sn) result(fn)

   implicit none
   integer,intent(in)   :: sn ! snapshot number
   character(len=255)   :: fn
   
   write(fn,'('//trim(para%snapshot_fmt)//')') trim(para%snapshot_prefix),sn

end function snfile

function filename(index,part1,part2,part3,multi) result(fn)

   ! returns the concatenated filename
   ! [part1][part2][part3].[index]
   ! if this file exists, otherwise (and only if index==0)
   ! [part1][part2][part3]
   ! if it exists
   ! the first or second format can be forced by setting multi

   implicit none
   integer,intent(in)               :: index
   character(*),intent(in)          :: part1
   character(*),intent(in),optional :: part2
   character(*),intent(in),optional :: part3
   logical,intent(in),optional      :: multi
   character(len=255)               :: filebase,filebasedot,fn
   character(len=5)                 :: str
   
   filebase = trim(part1)
   if (present(part2)) filebase = trim(filebase)//trim(part2)
   if (present(part3)) filebase = trim(filebase)//trim(part3)
   
   write(str,'(I5)') index
   filebasedot = trim(filebase)//'.'//adjustl(trim(str))
   
   if (present(multi)) then
      if (multi) then
         fn = trim(filebasedot)
      else
         if (index>0) call error('Error in function filename(): single file with index>0')
         fn = trim(filebase)
      end if
   else
      if (exists(trim(filebase)//'.0').or.(index>0)) then
         fn = trim(filebasedot)
      else
         fn = trim(filebase)
      end if
   end if
   
end function filename

function get_number_of_subfiles(filebase) result(nfiles)

   ! counts the number of files with filebase, i.e. [filebase].0, [filebase].1, ...

   implicit none
   character(*),intent(in) :: filebase
   character(len=255)      :: fn
   logical                 :: single_file_exists
   logical                 :: multiple_file_exists
   integer*4               :: i,nfiles
   
   ! determine if the snapshot is stored in one or multiple files
   single_file_exists = exists(trim(filename(0,filebase,multi=.false.)))
   multiple_file_exists = exists(trim(filename(0,filebase,multi=.true.)))
   if ((.not.single_file_exists).and.(.not.multiple_file_exists)) then
      call error('could not find any files for the base '//trim(filebase))
   else if (single_file_exists.and.multiple_file_exists) then
      call error('ambiguous, since single and multiple files exist for the base '//trim(filebase))
   else if (single_file_exists) then
      nfiles = 1
   else
      i = 1
      do
         fn = filename(i,filebase,multi=.true.)
         inquire(file=trim(fn),exist=multiple_file_exists)
         if (.not.multiple_file_exists) exit
         i = i+1
      end do
      if (i==1) call error('there is only one sub-file in the base '//trim(filebase))
      nfiles = i
   end if
   
end function get_number_of_subfiles

subroutine load_parameters

   implicit none
   
   call read_parameters
   
   call get_parameter_value(para%simulation,'simulation')
   call get_parameter_value(para%snapshot,'snapshot')
   call get_parameter_value(para%L,'L',min=1e-10)
   call get_parameter_value(para%N,'N',min=1)
   call get_parameter_value(para%path_gadget,'path_gadget')
   call get_parameter_value(para%path_velociraptor,'path_velociraptor')
   call get_parameter_value(para%path_surfsuite,'path_surfsuite')
   call get_parameter_value(para%path_analysis,'path_analysis')
   call get_parameter_value(para%ext_groups,'ext_groups')
   call get_parameter_value(para%ext_particles,'ext_particles')
   call get_parameter_value(para%ext_sorted,'ext_sorted')
   call get_parameter_value(para%ext_halos,'ext_halos')
   call get_parameter_value(para%ext_halolist,'ext_halolist')
   call get_parameter_value(para%snapshot_fmt,'snapshot_fmt')
   call get_parameter_value(para%snapshot_prefix,'snapshot_prefix')
   call get_parameter_value(para%file_scalefactors,'file_scalefactors')
   call require_no_parameters_left
   
   ! fix paths
   para%path_gadget = dir(para%path_gadget,ispath=.true.)
   para%path_velociraptor = dir(para%path_velociraptor,ispath=.true.)
   para%path_surfsuite = dir(para%path_surfsuite,ispath=.true.)
   para%path_analysis = dir(para%path_analysis,ispath=.true.)
   
   ! check paths and files
   if (.not.exists(para%path_gadget)) then
      call out('ERROR: Could not find directory specified by path_gadget:')
      call out(trim(para%path_gadget))
      call out('Consider specifying a different parameterset using "-parameterset" or')
      call out('changing the paths in the parameter file.')
      stop
   end if
   
   ! load scale factors
   call load_scalefactors
   
   ! make paths
   call make_path(para%path_surfsuite)
   call make_path(para%path_analysis)
   
   ! check permissions
   call check_file(para%path_surfsuite,'r')
   
end subroutine load_parameters

subroutine load_scalefactors

   implicit none
   
   integer*4   :: snapshot,snapshot_min,snapshot_max,snapshot_expected
   integer*4   :: status
   real*4      :: a,aold
   
   call check_file(para%file_scalefactors)
   
   ! determine minimum and maximum snapshot
   snapshot_min = huge(snapshot_min)
   snapshot_max = -huge(snapshot_max)
   open(1,file=trim(para%file_scalefactors),action='read',form='formatted')
   do
      read(1,*,IOSTAT=status) snapshot,a
      if (status>0) then
         call error('unknown format of scale factor file')
      else if (status<0) then
         exit ! end of file
      else
         snapshot_min = min(snapshot_min,snapshot)
         snapshot_max = max(snapshot_max,snapshot)
      end if
   end do
   close(1)

   ! read scale factors
   if (allocated(scalefactor)) deallocate(scalefactor)
   allocate(scalefactor(snapshot_min:snapshot_max))
   snapshot_expected = snapshot_min
   aold = 0.0
   open(1,file=trim(para%file_scalefactors),action='read',form='formatted')
   do
      read(1,*,IOSTAT=status) snapshot,a
      if (status>0) then
         call error('unknown format of scale factor file')
      else if (status<0) then
         exit ! end of file
      else
         if (snapshot.ne.snapshot_expected) call error('snapshots in scale factor file must be in increasing order')
         if (a<aold) call error('scale factors must be increasing with increasing snapshot index')
         scalefactor(snapshot) = a
         snapshot_expected = snapshot+1
      end if
   end do
   close(1)

end subroutine load_scalefactors

subroutine save_particles_sorted_format(index)

   implicit none
   integer,intent(in)   :: index
   character(len=255)   :: fn
   integer*8            :: i

   fn = filename(index,para%path_surfsuite,snfile(para%snapshot),para%ext_sorted)
   open(1,file=trim(fn),action='write',form='unformatted',status='replace',access='stream')
   write(1) (p(i),i=1,nparticles)
   close(1)
   
end subroutine save_particles_sorted_format

subroutine load_particles_sorted_format(index)

   implicit none
   integer,intent(in)   :: index
   character(len=255)   :: fn
   integer*8            :: i,file_size
   
   ! determine number of particles from file size
   fn = filename(index,para%path_surfsuite,snfile(para%snapshot),para%ext_sorted)
   inquire(file=trim(fn), size=file_size)
   nparticles = file_size/bytes_per_particle
   
   if (nparticles>0) then
   
      call out('Number of particles: ',nparticles)
      if (allocated(p)) deallocate(p)
      allocate(p(nparticles))
      open(1,file=trim(fn),action='read',form='unformatted',status='old',access='stream')
      read(1) (p(i),i=1,nparticles)
      close(1)
      
   else
      
      call error('Empty file '//trim(fn))
      deallocate(p)
      
   end if

end subroutine load_particles_sorted_format

end module module_io