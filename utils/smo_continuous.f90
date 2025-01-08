program smo_continuous
  implicit none

  real, allocatable:: n_k(:), n_k_old(:)

  integer:: n_particles
  real:: s_dt, s_dt_0
  real:: cluster_medio
  
  integer:: count_stable
  real:: time_stable, k_mean
  
  character:: arg_char*50, param*50, param_value*50
  integer:: arg, num_args, eq_position, letter

  real:: time, time_init, time_final, delta_image, next_image
  integer::count

  integer:: mass_i, mass_j, mass_k
  real:: coag_creation, coag_destruction, coag_const
  real:: frag_creation, frag_destruction, frag_const

  real:: area_expo, diff_expo
  
  real,allocatable::coagulation_kernel(:,:), fragmentation_kernel(:,:)
  
  num_args = command_argument_count()

  n_particles  = 2000
  s_dt_0       = 0.0001
  time_init    = 1
  time_final   = 1000
  coag_const   = 0.001
  frag_const   = 0.0
  diff_expo    = -1.0
  area_expo    = 0.0 
      
  do arg = 1, num_args
     call get_command_argument(arg, arg_char)
     eq_position = 0
     param = " "
     param_value = " "
     do letter = 1, 50
         if( arg_char(letter:letter) .eq. '=')then
            eq_position = letter
         end if
      end do
      param = arg_char(1:eq_position-1)
      param_value = arg_char(eq_position+1:50)
      !write(*, *) param, param_value
      if(param .eq. "n_particles") read(param_value,*) n_particles 
      if(param .eq. "s_dt_0")      read(param_value,*) s_dt_0
      if(param .eq. "time_init")   read(param_value,*) time_init
      if(param .eq. "time_final")  read(param_value,*) time_final
      if(param .eq. "coag_const")  read(param_value,*) coag_const
      if(param .eq. "frag_const")  read(param_value,*) frag_const
      if(param .eq. "diff_expo")  read(param_value,*) diff_expo
      if(param .eq. "area_expo")  read(param_value,*) area_expo
      
   end do

   write(*,*) "n_particles", n_particles
   write(*,*) "coag_const", coag_const
   write(*,*) "frag_const", frag_const
   
   
   allocate(n_k(n_particles), n_k_old(n_particles))
   allocate(coagulation_kernel(n_particles,n_particles), fragmentation_kernel(n_particles,n_particles))

   ! initial condition
   n_k = 0
   n_k_old = 0
  
   n_k(1) = n_particles
   n_k_old(1) = n_particles

   s_dt    = s_dt_0

   time_stable   = 0
   count_stable  = 0    

   k_mean        = cluster_medio(n_k, n_particles)
    
   time = time_init
    
   delta_image   = 0.1
   next_image    = 1

   count       = 0

   ! kernels

   do mass_i = 1, n_particles
      do mass_j = 1, n_particles
         coagulation_kernel(mass_i, mass_j) = (((mass_i**diff_expo + mass_j**diff_expo)/2.0) &
&              * ((mass_i**area_expo + mass_j**area_expo)/2.0)) 
         fragmentation_kernel(mass_i, mass_j) = 1.0
      end do
   end do

   open(10, file="smo_result.dat")
   
   do while (time < time_final)
      do mass_k = 1, n_particles
         coag_creation    = 0.0
         frag_destruction = 0.0
         do mass_i = 1, mass_k-1
            coag_creation    = coag_creation + n_k_old(mass_i) * n_k_old(mass_k - mass_i) * &
&                 coagulation_kernel(mass_i, mass_k - mass_i)
            frag_destruction = frag_destruction - n_k_old(mass_k) * &
&            fragmentation_kernel(mass_i, mass_k - mass_i) 
         end do
         coag_creation = coag_creation * coag_const
         frag_destruction = frag_destruction * frag_const
         
         coag_destruction = 0.0
         frag_creation    = 0.0
         do mass_j = 1, n_particles - mass_k
            coag_destruction = coag_destruction + n_k_old(mass_j) * n_k_old(mass_k) * coagulation_kernel(mass_j, mass_k)
            frag_creation    = frag_creation - n_k_old(mass_k + mass_j) * fragmentation_kernel(mass_j, mass_k)
         end do
         coag_destruction = coag_destruction * coag_const
         frag_creation    = frag_creation    * frag_const
         n_k(mass_k) = (0.5*(coag_creation + frag_destruction) - (coag_destruction + frag_creation)) * s_dt + n_k_old(mass_k) 
         
      end do
      n_k_old = n_k
      if(time > next_image)then
         k_mean    =  cluster_medio(n_k, n_particles)
         write(10, *) time, k_mean
         !write(*, '(A, F12.4, 4x, F10.3)', ADVANCE='NO') "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b",time, k_mean ! '(f15.4, 4x, f10.3)'
      end if
      
      if (k_mean > (n_particles - 1) .and. time_stable .eq. 0) then
         time_stable = time
      else if(k_mean > (n_particles - 1) .and. (time_stable .ne.  0))then
         count_stable = count_stable + 1
      end if
      if (count_stable > 5) exit

      if(count == 10)    s_dt = s_dt * 2
      if(count == 20)    s_dt = s_dt * 5
      if(count == 30)    s_dt = s_dt * 10
      if(count == 40)    s_dt = s_dt * 10
      if(count == 50)    s_dt = s_dt * 10
      if(count == 100)   s_dt = s_dt * 10
      if(count == 200)   s_dt = s_dt * 2
      if(count == 500)   s_dt = s_dt * 5
      if(count == 700)   s_dt = s_dt * 2
      if(count == 1000)  s_dt = s_dt * 5
      if(count == 5000)  s_dt = s_dt * 2
      if(count == 10000) s_dt = s_dt * 5

      if(mod(count,100) .eq. 0)then
         write(*, *) time, k_mean
         close(10)
         open(10, file="smo_result.dat", status="old", position="append")
      end if
      count = count + 1
      
      time = time + s_dt
   end do
close(10)  

end program smo_continuous



real function cluster_medio(n_k, n_particles)
  real:: n_k(n_particles)
  cluster_medio = n_particles/sum(n_k)
  
end function cluster_medio
