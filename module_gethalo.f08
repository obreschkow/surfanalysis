module module_gethalo

   use shared_module_core
   use shared_module_arguments
   use shared_module_hdf5
   use module_global
   use module_io
   use module_processhalo

   private 
   public   :: task_gethalo
   public   :: nhalos
   public   :: load_halo_properties ! only load structure "halo"
   public   :: load_all_halo_properties ! load array of structure "halo" for all halos
   public   :: load_halo_particles  ! load particle data of a halo

contains

subroutine task_gethalo

   implicit none
   integer*4                  :: haloid
   character(len=255)         :: outputfile
   integer*4                  :: outputformat ! 1=binary, 2=ascii, 3=hdf5
   logical                    :: subhalos,center
   type(type_halo)            :: halo
   
   call get_task_value(haloid)
      
   ! handle options
   call get_option_value(outputfile,'-outputfile','')
   call get_option_value(outputformat,'-outputformat',1,min=1,max=3)
   call get_option_value(subhalos,'-subhalos',.false.)
   call get_option_value(center,'-center',.false.)
   call require_no_options_left
   
   ! load halo
   call load_halo_properties(haloid,halo)
   call load_halo_particles(haloid,subhalos)
   if (center) call center_particles
   
   ! write/save halo
   if (isempty(outputfile)) then
      call hline
      call write_halo_properties(haloid,halo)
      call write_halo_particle_stats
      call hline
   else
      call save_halo(outputfile,outputformat,haloid,halo)
   end if
   
end subroutine task_gethalo

integer*4 function nhalos()

   implicit none
   character(len=255)   :: fn
   integer*8            :: filesize
    
   fn = trim(para%path_surfsuite)//trim(snfile(para%snapshot))//trim(para%ext_halolist)
   call check_file(fn)
   inquire(file=trim(fn),size=filesize)
   nhalos = int(filesize/int(bytes_per_halo,8),4)

end function nhalos

subroutine load_halo_particles(haloid,include_subhalos)

   implicit none
   integer*4,intent(in)          :: haloid
   logical,intent(in)            :: include_subhalos
   
   call load_particles(haloid,include_subhalos)
   
   contains

   recursive subroutine load_particles(haloid,include_subhalos,recursive_call)

      implicit none
      integer*4,intent(in)          :: haloid
      logical,intent(in)            :: include_subhalos
      logical,intent(in),optional   :: recursive_call ! this should only be true when called recursively from within this subroutine
      character(len=255)            :: fn
      integer*8                     :: position
      type(type_halo)               :: halo
      integer*4                     :: i
      
      call load_halo_properties(haloid,halo)
      
      if (.not.present(recursive_call)) then
         if (include_subhalos) then
            nparticles = halo%npart+halo%npartsub
         else
            nparticles = halo%npart
         end if
         if (allocated(p)) deallocate(p)
         allocate(p(nparticles))
         nparticles = 0
      end if
      
      if (halo%npart>0) then
         ! note: there are some halos that have 0 particles, but children with particles
         ! it's a weird feature of the velociraptor outputs
         position = bytes_per_particle*halo%offset+1
         fn = filename(halo%file,para%path_surfsuite,snfile(para%snapshot),para%ext_halos)
         open(1,file=trim(fn),action='read',form='unformatted',status='old',access='stream')
         read(1,pos=position) p(nparticles+1)
         do i = 2,halo%npart
            read(1) p(nparticles+i)
         end do
         nparticles = nparticles+halo%npart
      end if
   
      if (include_subhalos) then
         if (halo%firstchildid>-1) call load_particles(halo%firstchildid,include_subhalos,.true.)
         if (present(recursive_call).and.(halo%siblingid>-1)) call load_particles(halo%siblingid,include_subhalos,.true.)
      end if
   
      if (.not.present(recursive_call)) then
         if (nparticles.ne.size(p)) call error('in load_halo_particles: n.ne.nparticles')
      end if
   
   end subroutine load_particles

end subroutine load_halo_particles

subroutine load_halo_properties(haloid,halo)

   implicit none
   integer*4,intent(in)                   :: haloid
   type(type_halo),intent(out)            :: halo
   character(len=255)                     :: fn
   character(len=20)                      :: idstr
   integer*8                              :: position
   
   write(idstr,*) haloid
   if (haloid<1) call error('HaloID below 1: '//trim(idstr))
   if (haloid>nhalos()) call error('HaloID too large: '//trim(idstr))
   
   fn = trim(para%path_surfsuite)//trim(snfile(para%snapshot))//trim(para%ext_halolist)
   position = int(bytes_per_halo,8)*int(haloid-1,8)+1_8
   open(1,file=trim(fn),action='read',form='unformatted',status='old',access='stream')
   read(1,pos=position) halo
   close(1)

end subroutine load_halo_properties

subroutine load_all_halo_properties(halo)

   implicit none
   type(type_halo),allocatable,intent(out)   :: halo(:)
   character(len=255)                        :: fn
   
   allocate(halo(nhalos()))
   
   fn = trim(para%path_surfsuite)//trim(snfile(para%snapshot))//trim(para%ext_halolist)
   open(1,file=trim(fn),action='read',form='unformatted',status='old',access='stream')
   read(1) halo
   close(1)

end subroutine load_all_halo_properties

subroutine save_halo(outputfile,outputformat,haloid,halo)

   implicit none
   character(len=255),intent(in) :: outputfile
   integer*4,intent(in)          :: outputformat ! 3=binary, 2=ascii, 1=hdf5
   integer*4,intent(in)          :: haloid
   type(type_halo),intent(in)    :: halo
   integer*8                     :: i
   
   if (outputformat==3) then
   
      open(1,file=trim(outputfile),action='write',form='unformatted',status='replace',access='stream')
      write(1) (p(i),i=1,nparticles)
      close(1)
      
   else if (outputformat==2) then
   
      open(1,file=trim(outputfile),action='write',form='formatted',status='replace')
      write(1,'(I11,I2,3E17.10,3E14.6)') (p(i),i=1,nparticles)
      close(1)
      
   else if (outputformat==1) then
      
      call hdf5_create(trim(outputfile)) ! create HDF5 file
      call hdf5_open(trim(outputfile),.true.) ! open HDF5 file
   
      ! Group "simulation"
      call hdf5_add_group('simulation')
      call hdf5_write_data('simulation/parameterset',trim(para%parameterset))
      call hdf5_write_data('simulation/name',trim(para%simulation),'simulation name')
      call hdf5_write_data('simulation/box_l',para%L,'[simulation units] box side length')
      call hdf5_write_data('simulation/box_n',para%N,'cubic root of particle number')
      call hdf5_write_data('simulation/snapshot',para%snapshot,'snapshot index')
      
      ! Group "surfsuite"
      call hdf5_add_group('surfsuite')
      call hdf5_write_data('surfsuite/timestamp',timestamp(),'surfsuite timestamp')
      call hdf5_write_data('surfsuite/version',trim(version),'version of surfsuite used to extract the halo')
      call hdf5_write_data('surfsuite/developer','Danail Obreschkow; danail.obreschkow@icrar.org')
      
      ! Group "halo"
      call hdf5_add_group('halo')
      call hdf5_write_data('halo/id',haloid,'unique halo id')
      call hdf5_write_data('halo/parentid',halo%parentid,'id of parent halo')
      call hdf5_write_data('halo/n_subhalos',halo%nchildren,'number of subhalos')
      call hdf5_write_data('halo/file',halo%file,'surfsuite file number')
      call hdf5_write_data('halo/n_particles_main',halo%npart,'number of particles in main halo')
      call hdf5_write_data('halo/n_particles',nparticles,'number of particles in file')
      call hdf5_write_data('halo/includes_subhalos',log2int(nparticles.ne.halo%npart), &
      & 'logical flag (0/1 = subhalo particles are excluded/included)')
      
      ! Group "particles"
      call hdf5_add_group('particles')
      call hdf5_write_data('particles/id',p%id,'unique particle id')
      call hdf5_write_data('particles/species',p%species,'particle species (1 = gas, 2 = dark matter)')
      call hdf5_write_data('particles/rx',p%x(1),'[simulation units] x-coordinate of position')
      call hdf5_write_data('particles/ry',p%x(2),'[simulation units] y-coordinate of position')
      call hdf5_write_data('particles/rz',p%x(3),'[simulation units] z-coordinate of position')
      call hdf5_write_data('particles/vx',p%v(1),'[simulation units] x-coordinate of velocity')
      call hdf5_write_data('particles/vy',p%v(2),'[simulation units] y-coordinate of velocity')
      call hdf5_write_data('particles/vz',p%v(3),'[simulation units] z-coordinate of velocity')
      
      call hdf5_close() ! close HDF5 file
      
   else
      call error('unknown outputformat.')
   end if
end subroutine save_halo

subroutine write_halo_properties(haloid,halo)
   implicit none
   integer*4,intent(in)                   :: haloid
   type(type_halo),intent(in)             :: halo
   call out('Parameterset:  ',trim(para%parameterset))
   call out('Simulation:    ',trim(para%simulation))
   call out('Snapshot:      ',para%snapshot)
   call out('Halo:          ',haloid)
   call out('File:          ',halo%file)
   call out('Offset:        ',halo%offset)
   call out('Npart:         ',halo%npart)
   call out('Npart sub:     ',halo%npartsub)
   call out('Parent ID:     ',halo%parentid)
   call out('Nchidren:      ',halo%nchildren)
   call out('1st child ID:  ',halo%firstchildid)
   call out('Sibling ID:    ',halo%siblingid)
end subroutine write_halo_properties

subroutine write_halo_particle_stats
   implicit none
   character(len=255)   :: txt
   character(len=12)    :: wrap
   real*4               :: h
   if (nparticles>0) then
      
      write(txt,'(A,I0,A,I0,A,I0,A)')     'ID range:      ', &
      &minval(p%id),' to ',maxval(p%id),' (first: ',p(1)%id,')'
      call out(txt)
      write(txt,'(A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A)') &
      & '# particles:   ',nparticles,' (',count(p%species==1),',',count(p%species==2),',',count(p%species==3),',',&
      & count(p%species==4),',',count(p%species==5),',',count(p%species==6),')'
      call out(txt)
      
      h = para%L/2
      if (maxval(p%x(1))-minval(p%x(1))>h) then
         wrap = '   (wrapped)'
      else
         wrap = ''
      end if
      write(txt,'(A,F0.4,A,F0.4,A)') 'x-range:       ',minval(p%x(1)),' to ',maxval(p%x(1)),wrap
      call out(txt)
      if (maxval(p%x(2))-minval(p%x(2))>h) then
         wrap = '   (wrapped)'
      else
         wrap = ''
      end if
      write(txt,'(A,F0.4,A,F0.4,A)') 'y-range:       ',minval(p%x(2)),' to ',maxval(p%x(2)),wrap
      call out(txt)
      if (maxval(p%x(3))-minval(p%x(3))>h) then
         wrap = '   (wrapped)'
      else
         wrap = ''
      end if
      write(txt,'(A,F0.4,A,F0.4,A)') 'z-range:       ',minval(p%x(3)),' to ',maxval(p%x(3)),wrap
      call out(txt)
      write(txt,'(A,F0.4,A,F0.4,A)') 'vx-range:      ',minval(p%v(1)),' to ',maxval(p%v(1)),wrap
      call out(txt)
      write(txt,'(A,F0.4,A,F0.4,A)') 'vy-range:      ',minval(p%v(2)),' to ',maxval(p%v(2)),wrap
      call out(txt)
      write(txt,'(A,F0.4,A,F0.4,A)') 'vz-range:      ',minval(p%v(3)),' to ',maxval(p%v(3)),wrap
      call out(txt)
      
   end if
end subroutine write_halo_particle_stats
   
end module module_gethalo