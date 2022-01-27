module module_global

   type type_para
      character(len=255)   :: parameterfile
      character(len=255)   :: parameterset
      character(len=255)   :: simulation
      integer*4            :: snapshot
      real*4               :: L ! box side length in simulation units
      integer*4            :: N ! number of particles on a side
      character(len=255)   :: path_gadget
      character(len=255)   :: path_velociraptor
      character(len=255)   :: path_surfsuite
      character(len=255)   :: path_analysis
      character(len=255)   :: ext_groups
      character(len=255)   :: ext_particles
      character(len=255)   :: ext_sorted
      character(len=255)   :: ext_halos
      character(len=255)   :: ext_halolist
      character(len=255)   :: snapshot_fmt
      character(len=255)   :: snapshot_prefix
      character(len=255)   :: file_scalefactors
   end type type_para

   type(type_para)         :: para

   type type_halo
      integer*4   :: file
      integer*4   :: offset
      integer*4   :: npart    ! number of particles in halo without subhalos
      integer*4   :: npartsub ! number or particles in subhalos (without parent)
      integer*4   :: parentid
      integer*4   :: nchildren
      integer*4   :: firstchildid
      integer*4   :: siblingid
   end type type_halo
   
   type type_particle
      integer*8            :: id          ! unique identifiers
      integer*4            :: species     ! particle species
      real*4               :: x(3)        ! positions
      real*4               :: v(3)        ! velocities
   end type type_particle

   integer*4,parameter     :: bytes_per_particle = 36 ! bytes per instance of type_particle
   integer*4,parameter     :: bytes_per_halo = 32 ! bytes per instance of type_halo
   integer*8,parameter     :: nparticles_per_sorted_file = 100000000_8

   ! currently loaded particles (can be a subvolume, a single halo, substructure, etc.)
   type(type_particle),allocatable  :: p(:)        ! array of particle properties
   integer*8                        :: nparticles  ! number of particles stored in p(:)
   
   ! scale factors of different snapshots
   real*4,allocatable      :: scalefactor(:)
   
   ! physical/mathematical constants
   real*4,parameter        :: pi = 3.14159265358979
   
end module module_global