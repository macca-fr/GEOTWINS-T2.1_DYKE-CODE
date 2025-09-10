      Program main_DYKE_PROPAGATION
      Use PROPAGATION,    Only : SET_PROP_PAR, START_PROP
      Use COMP_FIELD,     Only : SET_OUTPUT, COMPUTE_EX_STRESS_DIRECTIONS, COMPUTE_LITH_P_PROFILE
      Use BELEMENT,       Only : SET_LITH_PAR, SET_IDR_PAR, SET_TK_PAR, SET_IN_CRACK_GEOM
      Use EXTERNAL_FIELD, Only : SET_EX_STRESS_MATRIX, DEALLOC_EX_STRESS_MATRIX

      IMPLICIT NONE

      REAL(8) :: pi=4.D0*datan(1.D0)

      REAL(8) :: delta_out,scl
      LOGICAL :: comp_dyke_stress,dikeRF_output,movie_
       
      REAL(8) :: mu1,mu2,nu1,nu2
      REAL(8) :: rhoR1,rhoR2,z_dns,z_ref
      LOGICAL :: Lt_stress_on
      
      REAL(8) :: rhoM,Kappa,A0,DeA0,melt_frac,melt_in_dyke
      
      REAL(8) :: SigXX,SigZZ,mx_SXX,mz_SXX,mx_SZZ,mz_SZZ
      LOGICAL :: Tk_stress_on
      
      REAL(8), ALLOCATABLE :: x0z0_in(:,:),delta_in(:),Hd_in(:)
      REAL(8) :: delta0,Hd0,x_top,z_top,xf_top,zf_top
      INTEGER :: n_dis_in, n_dis
      LOGICAL :: delta_optimal

      INTEGER :: n_iter
      INTEGER :: n_dir
      REAL(8) :: alfa
      REAL(8) :: x_stop_min,x_stop_max,z_stop
      REAL(8) :: E_prop

      INTEGER :: i,run,dyke_run,n_dykes

      REAL(8), ALLOCATABLE :: X_start(:), Z_start(:), A_0(:), DeA_0(:), delta_0(:)

      CHARACTER*50 :: inpf00,inpf01,inpf02,outf00,outf01,input_stress
      
      REAL(8) :: x_min,x_max,z_min,z_max,Dz_snap,out_res,grd_1,grd_2
      
      inpf00 = './input/input_field.dat'
      inpf01 = './input/input_BE.dat'
      Open(40,file = inpf00,access='sequential',form='formatted')
      Open(41,file = inpf01,access='sequential',form='formatted')

!~       outf00 = './output/1_arrival.dat'
!~       outf01 = './output/2_stopped.dat'
!~       Open(101,file = outf00,access='sequential',form='formatted',status='new')
!~       Open(102,file = outf01,access='sequential',form='formatted',status='new')

!-------------------------------------------------------------------!
! SET OUTPUT PARAMETERS
      Read(40,*)                              ! first row of input_stress file
      Read(40,*) input_stress                 ! name of file with the interpolated background stress field
      Read(40,*)                              !
      Read(40,*) x_min, x_max, z_min, z_max   ! spatial domain for the input stress file [km] or [m]
      Read(40,*) grd_1				          ! grid step for the input stress file [km] or [m]
      Read(40,*)                              !
      Read(40,*) Dz_snap, comp_dyke_stress    ! z-interval (depth-interval) for calculation of dyke shape and stress-displacements fields / logical variable: if .false. the stress and displacement due to the dyke is not printed in output
      Read(40,*)                              !
      Read(40,*) dikeRF_output, out_res, grd_2! if dikeRF is .true. the dyke stress is computed in a ref frame with z parallel to the dike; out_res is the spatial resolution for output stress and displ calculation; grd_2 is the grid step for plotting the principal stress [km] or [m]
      Read(40,*)                              !
      Read(40,*) scl                          ! scaling factor for plotting the dyke shape

      movie_ = .false.!.true.!
      delta_out = 0.D0

      Call SET_OUTPUT(x_min, z_min, x_max, z_max,out_res,movie_,comp_dyke_stress,dikeRF_output,delta_out,scl)
      Call SET_EX_STRESS_MATRIX(x_min,x_max,z_min,z_max,grd_1,input_stress)
      Call COMPUTE_EX_STRESS_DIRECTIONS(x_min,x_max,z_min,z_max,grd_2)
!-------------------------------------------------------------------!

!-------------------------------------------------------------------!
! SET ROCK PARAMETERS

      Read(41,*)                     ! first  row of input_BE file
      Read(41,*)                     ! row 2
      Read(41,*) mu1, nu1, rhoR1     ! row 3 [Mpa] rigidity of the upper lower layer (z>0) / Poisson's ratio (z>0) / [kg/dm^3] rock density in the lower layer (z>0)

      Read(41,*)                     ! row 4
      Read(41,*) mu2, nu2, rhoR2     ! row 5 [Mpa] rigidity of the upper upper layer (z<0) / Poisson's ratio (z<0) / [kg/dm^3] rock density in the upper layer (z<0) 
      Read(41,*)                     ! row 6
      Read(41,*) z_dns, z_ref 		 ! row 7 - z_dns [km] depth of density transition / z_ref [km] depth of the reference surface (coincide with elastic free surface only if mu2=0 and z_ref=0,
                                     !         otherwise represent the z at which the lithostatic pressure vanishes) 

      Lt_stress_on = .true. ! if .true. the dike cross section lay on a vertical plane, if .false. the cross section is horizontal and
                            ! the lithostatic pressure is constant and given by z_ref, that became the depth of the dike cross section.

      Call SET_LITH_PAR(mu1,mu2,nu1,nu2,rhoR1,rhoR2,z_dns,z_ref,Lt_stress_on)
      Call COMPUTE_LITH_P_PROFILE(z_min,z_max,grd_1)
!-------------------------------------------------------------------!

!-------------------------------------------------------------------!
! SET THE MAGMA PARAMETERS

      Read(41,*)             ! row 8
      Read(41,*) rhoM, Kappa ! row 9 [kg/dm^3] magma density / [MPa] magma Bulk modulus


!-------------------------------------------------------------------!

!-------------------------------------------------------------------!
! SET THE TECTONIC STRESS PARAMETERS 
 
      Tk_stress_on = .true.!.false.     ! if .false. the tectonic stress is not considered

      Read(41,*)                                         ! row 10
      Read(41,*) SigXX,SigZZ,mx_SXX,mz_SXX,mx_SZZ,mz_SZZ ! row 11 stress gradient parameters (MPa)

       Call SET_TK_PAR(SigXX,mx_SXX,mz_SXX,SigZZ,mx_SZZ,mz_SZZ,Tk_stress_on)
!~       Call SET_TK_PAR(0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,Tk_stress_on)
!-------------------------------------------------------------------!

!-------------------------------------------------------------------!
      Read(41,*)                                              ! row 12
      Read(41,*) n_dykes                                      ! row 13 number of dikes for the current run
      Read(41,*)                                              ! row 14

      ALLOCATE (X_start(n_dykes),Z_start(n_dykes),A_0(n_dykes),DeA_0(n_dykes),delta_0(n_dykes))
      Do i=1,n_dykes
        Read(41,*) X_start(i), Z_start(i), A_0(i), DeA_0(i), delta_0(i) ! row 15 (x,z) initial positions of the dyke tip h_melt, cross sectional area, dip angle delta_0 (degrees from vertical positive counter-clockwise, negative clockwise)
      EndDo
      Read(41,*)                                              ! row n_dykes+15
      Read(41,*) n_dis_in, Hd0                                ! row n_dykes+16 initial number of elementary dislocations / elements length
      Read(41,*)                                              ! row n_dykes+17
      Read(41,*) delta_optimal                                ! row n_dykes+18 if true the dyke start optimally oriented if false start with angle delta0 (degrees from vertical positive counter-clockwise, negative clockwise)
      Read(41,*)                                              ! row n_dykes+19
      Read(41,*) n_iter                                       ! row n_dykes+20 max number of iterations per dyke
      Read(41,*)                                              ! row n_dykes+21
      Read(41,*) n_dir, alfa                                  ! row n_dykes+22 number of tested direction / angle between tested directions
      Read(41,*)                                              ! row n_dykes+23
      Read(41,*) x_stop_min, x_stop_max, z_stop               ! row n_dykes+24 boundaries for the dyke propagation paths
      Read(41,*)                                              ! row n_dykes+25
      Read(41,*) E_prop                                       ! row n_dykes+26 energy threashold for propagation

      alfa = alfa*pi/180.D0

      dyke_run = 1 
      Do run = 1,n_dykes ! cycle for performing more runs (in this case same parameters but different starting points and A0)
!-------------------------------------------------------------------!
!-------------------------------------------------------------------!
! SET THE INITIAL GEOMETRY OF THE CRACK
      
      ALLOCATE (x0z0_in(n_dis_in,2),delta_in(n_dis_in),Hd_in(n_dis_in))
 
      x_top = X_start(run)                         ! initial x of the upper tip of the crack
      z_top = Z_start(run)                         ! initial z of the upper tip of the crack
      A0    = A_0(run)                             ! reference cross sectional area of the crack 
      DeA0  = DeA_0(run)                           ! increment of cross sectional area of the crack 
      delta_in = delta_0(run)*pi/180.D0

      Hd_in    = Hd0

      Call SET_IDR_PAR(rhoM,Kappa,A0)
      Call SET_IN_CRACK_GEOM(x0z0_in,delta_in,Hd_in,n_dis_in,x_top,z_top,delta_optimal)
!-------------------------------------------------------------------!

!-------------------------------------------------------------------!
! SET THE PARAMETERS NEEDED FOR THE DIKE PROPAGATION

      Call SET_PROP_PAR(n_iter,n_dir,alfa,x_stop_min,x_stop_max,z_stop, &
     &                  Dz_snap,E_prop,E_prop,E_prop)
!-------------------------------------------------------------------!

! START DIKE PROPAGATION
      n_dis = n_dis_in
      
      If ( (x_top.lt.x_min)                                             &
     &           .or.                                                   &
     &     (x_top.gt.x_max)                                             &
     &           .or.                                                   &
     &     (z_top.lt.z_min)                                             &
     &           .or.                                                   &
     &     (z_top.gt.z_max) ) Then
     
        Write(*,*) '--------------------------------------------------------------------'
        Write(*,*) '--------------------------- WARNING --------------------------------'
        Write(*,*) '--------------------------------------------------------------------'
        Write(*,*) 'The dyke number ', run, ' of ', n_dykes, ' will not be simulated'
        Write(*,*) 'because its starting point is not located within the elastic domain.'
        Write(*,*) 'The actual total number of dyke will be ', n_dykes-(run-dyke_run+1)
        Write(*,*) '--------------------------------------------------------------------'
        Write(*,*) '------------------ PRESS ANY KEY TO CONTINUE -----------------------'
        Write(*,*) '--------------------------------------------------------------------'
        Read(*,*)

      Else

        CALL START_PROP(x0z0_in,delta_in,Hd_in,A0,DeA0,n_dis,dyke_run,xf_top,zf_top)
        dyke_run = dyke_run + 1

      EndIf

!-------------------------------------------------------------------!

      DEALLOCATE (x0z0_in,delta_in,Hd_in)
      DEALLOCATE (X_start,Z_start,A_0,DeA_0,delta_0)

!-------------------------------------------------------------------!
      End Do

      Call DEALLOC_EX_STRESS_MATRIX()

      End Program main_DYKE_PROPAGATION
