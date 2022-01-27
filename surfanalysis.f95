program surfanalysis

   use shared_module_core
   use shared_module_arguments
   use module_io
   use module_moments
   use module_asymmetries

   implicit none
   
   ! start user interface
   call set_version('0.24')
   call handle_arguments(require_task=.true.)
   call start_output
   
   ! handle general options
   call get_option_value(para%parameterfile,'-parameterfile','parameters.txt')
   call set_parameterfile(para%parameterfile)
   call get_option_value(para%parameterset,'-parameterset','')
   call set_parameterset(para%parameterset)
   
   ! load parameters
   call load_parameters
   
   ! overwrite default snapshot, if option -snapshot given
   call get_option_value(para%snapshot,'-snapshot',int(para%snapshot,4))
   call check_file(para%path_analysis,'rw')

   ! run task
   if (istask('moments',.false.)) then
      call task_moments
   else if (istask('asymmetries',.false.)) then
      call task_asymmetries ! for Jie's work
   else
      call unknown_task
   end if
   
   ! finalize output on screen/logfile
   call stop_output
    
end program surfanalysis