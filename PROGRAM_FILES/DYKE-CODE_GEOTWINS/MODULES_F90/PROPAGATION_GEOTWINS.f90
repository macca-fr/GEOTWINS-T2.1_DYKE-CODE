      Module PROPAGATION
      Use COMP_FIELD
      Use BELEMENT

      INTEGER, PRIVATE :: n_iter
      INTEGER, PRIVATE :: n_dir_in
      REAL(8), PRIVATE :: alfa
      REAL(8), PRIVATE :: x_stop_min,x_stop_max,z_stop,Dz_snap
      REAL(8), PRIVATE :: E_threshold_1,E_threshold_2,E_threshold_int
      LOGICAL, PRIVATE :: int_dir,check_int
      
      REAL(8), ALLOCATABLE, PRIVATE :: x0z0_old(:,:),delta_old(:),Hd_old(:)
      REAL(8), PRIVATE :: E_old
      INTEGER, PRIVATE :: n_dis_old

      REAL(8), PRIVATE :: pi=4.D0*datan(1.D0)
      
      contains

!~ ---------------------------------------------------------------------
!~ ---------------------------------------------------------------------

      Subroutine SET_PROP_PAR(n_iter_,n_dir_,alfa_,x_stop_min_,x_stop_max_,z_stop_,   &
     &                        Dz_snap_,E_threshold_1_,E_threshold_2_,E_threshold_int_)

      IMPLICIT NONE
      
      INTEGER :: n_iter_
      INTEGER :: n_dir_
      REAL(8) ::  alfa_
      REAL(8) :: x_stop_min_,x_stop_max_,z_stop_,Dz_snap_
      REAL(8) :: E_threshold_1_,E_threshold_2_,E_threshold_int_

      n_iter  = n_iter_

      n_dir_in = n_dir_
      alfa  = alfa_

      x_stop_min  = x_stop_min_
      x_stop_max  = x_stop_max_
      z_stop      = z_stop_
      Dz_snap     = Dz_snap_

      E_threshold_1   = E_threshold_1_
      E_threshold_2   = E_threshold_2_
      E_threshold_int = E_threshold_int_

      int_dir = .false.
      If ((E_threshold_1.eq.E_threshold_2).and.(E_threshold_1.eq.E_threshold_int)) Then
        check_int = .false.
      Else
        check_int = .true.
      EndIf

      End Subroutine SET_PROP_PAR

!~ ---------------------------------------------------------------------
!~ ---------------------------------------------------------------------

      Subroutine START_PROP(x0z0_in,delta_in,Hd_in,A0,DeA0,n_dis,run,xf_top,zf_top)
!~ x0z0_in  INPUT
!~ delta_in INPUT
!~ Hd_in    INPUT
!~ A0       INPUT
!~ n_dis    INPUT/OUTPUT
!~ run      INPUT
!~ xf_top  OUTPUT
!~ zf_top  OUTPUT
      IMPLICIT NONE

      REAL(8) :: x0z0_in(:,:),delta_in(:),Hd_in(:)
      REAL(8) :: A0,DeA0,xf_top,zf_top
      INTEGER :: n_dis,run

      REAL(8), ALLOCATABLE :: Ut(:),    Ub(:),    x0z0(:,:),    delta(:),    Hd_(:)
      REAL(8), ALLOCATABLE :: Ut_old(:),Ub_old(:)
      REAL(8) :: A,A_split,A_tot,A_old,E,DeE,OverP_av,rho_iter,Ub_max,Ub_min
      REAL(8) :: x_top,z_top,x_start,z_start,z_snap
!~ DELTA -------------------------------------------------------------
      REAL(8) :: L_norm,dip_av,dip_ndis,dip_max,Sig_tip(2,2),           &
     &           Sig_xz_av,Sig_xz_n,Sig_xz_max,DELTA1,DELTA2,DELTA3
!~ -------------------------------------------------------------------
!~ Plit sig_T sig_N ------------------------------------------------
      REAL(8) :: Plit,sig_N,sig_T
!~ -------------------------------------------------------------------
      INTEGER :: n_dir,iter,i,l,PROP_STEP
      LOGICAL :: esc,esc_old,dyke_split
      CHARACTER(50) :: outf01,outf02,outf03,outf04,outf05
      CHARACTER(2)  :: n_outf

      REAL(8), ALLOCATABLE :: P_lit(:),OverP(:)

      Write(n_outf,20) run
      outf01 = './output/XtZtDeE_'//n_outf//'.dat'
      outf02 = './output/RUNINFO_'//n_outf//'.dat'
!~ DELTA -------------------------------------------------------------
!~       outf03 = './output/XtZtDELTAINF_'//n_outf//'.txt'
!~       outf04 = './output/XtZtDELTA_'//n_outf//'.dat'
!~ -------------------------------------------------------------------
!~ M-FLUX -------------------------------------------------------------
!~       outf03 = './output/XtZtDELTAINF_'//n_outf//'.txt'
      outf04 = './output/M-FLUX_'//n_outf//'.dat'
!~ -------------------------------------------------------------------
      outf05 = './output/XtZtPlitSigNSigT_'//n_outf//'.dat'

      Open(23,file = outf01,access='sequential',form='formatted',status='new')
      Open(22,file = outf02,access='sequential',form='formatted',status='new')
!~ DELTA -------------------------------------------------------------
!~       Open(23,file = outf03,access='sequential',form='formatted',status='new')
!~       Open(24,file = outf04,access='sequential',form='formatted',status='new')
!~ -------------------------------------------------------------------
!~ M-FLUX -------------------------------------------------------------
!~       Open(23,file = outf03,access='sequential',form='formatted',status='new')
      Open(24,file = outf04,access='sequential',form='formatted',status='new')
      Write (24,*) 'ITER,  PROP_STEP,  A0,  A_iter, rho_iter,  DeE,  n_dis,  DP,  Xt,  Zt,  PROP(Y/N)'
!~ END M-FLUX ---------------------------------------------------------

!~ -------------------------------------------------------------------
      Open(25,file = outf05,access='sequential',form='formatted',status='new')

      E = 0.D0
      A = A0

      ALLOCATE (Ut(n_dis),Ub(n_dis),x0z0(n_dis,2),delta(n_dis),Hd_(n_dis))

      x0z0  = x0z0_in
      delta = delta_in
      Hd_   = Hd_in

      Write (22,*) 'n_dis input  =', n_dis
      Write (*,*)  'n_dis input  =', n_dis
      Write (22,*) 'A0  =', A
      Write (*,*)  'A0  =', A

      Call SOLVE_CRACK(Ut,Ub,x0z0,delta,Hd_,A,n_dis,dyke_split,esc)
      Call COMP_ENERGY(Ut,Ub,x0z0,delta,Hd_,A,n_dis,E)
      Write (22,*) 'Etot =',E
      
      x_top = x0z0(n_dis,1) - Hd_(n_dis)*sin(delta(n_dis))*0.5d0
      z_top = x0z0(n_dis,2) - Hd_(n_dis)*cos(delta(n_dis))*0.5d0
      
      x_start = x_top
      z_start = z_top

      z_snap = z_top - Dz_snap

      Do i=1,n_dis
!        Write (21,21) x0z0(i,1),x0z0(i,2),0.D0
        Write (23,21) x0z0(i,1),x0z0(i,2),0.D0
      EndDo
      Write (22,*) 'n_dis output =', n_dis
      Write (*,*)  'n_dis output =', n_dis
      Write (22,*) 'A  =', A
      Write (*,*)  'A  =', A
      Write (22,*) 'delta_in  =', delta(n_dis)/pi*180.D0
      Write (*,*)  'delta_in  =', delta(n_dis)/pi*180.D0
      Write (22,*) '(x_start;z_start)=(',x_start,';',z_start,')'
      Write (*,*)  '(x_start;z_start)=(',x_start,';',z_start,')'
      
      Call DISPL_FIELD(Ut,Ub,x0z0,delta,Hd_,n_dis,run,0)
      Call STRSS_FIELD(Ut,Ub,x0z0,delta,Hd_,n_dis,run,0)
      Call CRACK_SHAPE(Ut,Ub,x0z0,delta,Hd_,n_dis,run,0)

!~ CALCOLO OVERPRESSURE AND MAGMA DENSITY AT DEPTH -------------------
      ALLOCATE (OverP(n_dis))
      Call COMP_OverP(Ut,Ub,x0z0,delta,Hd_,A,n_dis,OverP,OverP_av,rho_iter)
      Write (22,*) 'OverP_av  =', OverP_av
      Write (*,*)  'OverP_av  =', OverP_av
      Write (22,*) 'OverP_max =', OverP(n_dis)
      Write (*,*)  'OverP_max =', OverP(n_dis)
      Write (24,*) 0, 0, A0, A, rho_iter, '  XXXXXXXXXXXXXXXXXXXXXXX', n_dis, OverP(n_dis), x_top, z_top, 'N'
      DEALLOCATE (OverP)
!~ FINE CALCOLO OVERPRESSURE -----------------------------------------

!~ DELTA -------------------------------------------------------------
!~         Call COMP_DELTA_iter(x_top,z_top,delta,Hd_,n_dis,z_start,       &
!~      &                           OverP_av,L_norm,                       &
!~      &                           dip_av,dip_ndis,dip_max,               &
!~      &                           Sig_tip,Sig_xz_av,Sig_xz_n,Sig_xz_max, &
!~      &                           DELTA1,DELTA2,DELTA3)
!~         Write (24,22) x_top,z_top,OverP_av,L_norm,                      &
!~        &              Sig_xz_av,Sig_xz_n,Sig_xz_max,                    &
!~        &              DELTA1,DELTA2,DELTA3
!~         Write (23,*) '--------------------------------'
!~         Write (23,*) '---------- ITER = 0 -------'
!~         Write (23,*) '--------------------------------'
!~         Write (23,*) 'x_top, z_top    = ',x_top,z_top
!~         Write (23,*) 'S_xx, S_zz, Sxz = ',Sig_tip(1,1),Sig_tip(2,2),Sig_tip(1,2)
!~         Write (23,*) 'dip_av          = ',dip_av
!~         Write (23,*) 'dip_ndis        = ',dip_ndis
!~         Write (23,*) 'dip_max         = ',dip_max
!~         Write (23,*) 'S_xz_av         = ',Sig_xz_av
!~         Write (23,*) 'S_xz_n          = ',Sig_xz_n
!~         Write (23,*) 'S_xz_max        = ',Sig_xz_max
!~         Write (23,*) 'OverP_av        = ',OverP_av
!~         Write (23,*) 'L_norm          = ',L_norm
!~         Write (23,*) 'DELTA1, DELTA2, DELTA3  = ',DELTA1,DELTA2,DELTA3
!~         Write (23,*) ''
!~ END DELTA ---------------------------------------------------------


!~ COMPUTE P_lit sig_N sig_T -----------------------------------------
        Call COMP_LITH_P (Plit,x0z0(n_dis,2))
        Call COMP_sigN_sigT (sig_N,sig_T,x0z0(n_dis,1),x0z0(n_dis,2),delta(n_dis)) 
        Write (25,21) x0z0(n_dis,1),x0z0(n_dis,2),Plit,sig_N,sig_T
!~ END COMPUTE P_lit sig_N sig_T -------------------------------------

      PROP_STEP = 0

      Do iter = 1,n_iter

        PROP_STEP = PROP_STEP + 1
        
        Write (22,*) '**************************************************'
        Write (*,*)  '**************************************************'
        Write (22,*) 'ITER =',iter,'/',n_iter
        Write (*,*)  'ITER =',iter,'/',n_iter

!~         If (iter.le.1) Then
!~           n_dir = 61
!~         Else
!~           If (iter.le.10) Then
!~             n_dir = 21
!~           Else
            n_dir = n_dir_in
!~           EndIf
!~         EndIf

        ALLOCATE (Ut_old(n_dis),Ub_old(n_dis),x0z0_old(n_dis,2),delta_old(n_dis),Hd_old(n_dis))
      
        Ut_old(1:n_dis)     = Ut(1:n_dis)
        Ub_old(1:n_dis)     = Ub(1:n_dis)
        x0z0_old(1:n_dis,:) = x0z0(1:n_dis,:)
        delta_old(1:n_dis)  = delta(1:n_dis)
        Hd_old(1:n_dis)     = Hd_(1:n_dis)
        
        A_old     = A
        E_old     = E
        n_dis_old = n_dis

        DEALLOCATE (Ut,Ub,x0z0,delta,Hd_)
        ALLOCATE (Ut(n_dis+1),Ub(n_dis+1),x0z0(n_dis+1,2),delta(n_dis+1),Hd_(n_dis+1))
 
        x0z0(1:n_dis,:) = x0z0_old(1:n_dis,:)
        delta(1:n_dis)  = delta_old(1:n_dis)
        Hd_(1:n_dis)    = Hd_old(1:n_dis)

!~         If (check_int.eqv..false.) Then
        Call TEST_DIR(Ut,Ub,x0z0,delta,Hd_,A,n_dis,E,n_dir,alfa,dyke_split,esc)
!~         Else
!~           If (sign(1.D0,(z_top-Hd_(n_dis)*sin(delta(n_dis))))/sign(1.D0,z_top).lt.0.D0) Then
!~             Call TEST_DIR_INT(Ut,Ub,x0z0,delta,Hd_,A,n_dis,E,n_dir,alfa,esc)
!~           Else
!~             Call TEST_DIR(Ut,Ub,x0z0,delta,Hd_,A,n_dis,E,n_dir,alfa,esc)
!~           EndIf
!~         EndIf

        x_top = x0z0(n_dis,1) - Hd_(n_dis)*sin(delta(n_dis))*0.5d0
        z_top = x0z0(n_dis,2) - Hd_(n_dis)*cos(delta(n_dis))*0.5d0

!~         If ((delta(n_dis).eq.pi*0.D0).and.(int_dir.eqv..true.)) Then
!~           DeE = (E_old - E)/Hd_(n_dis) - E_threshold_int ! MPa*m ( E = U+W [MPa*km^2*10^3] )
!~         Else
!~           If (z_top.gt.0.D0) Then
            DeE = (E_old - E)/Hd_(n_dis) - E_threshold_1 ! MPa*m ( E = U+W [MPa*km^2*10^3] )
!~           Else
!~             DeE = (E_old - E)/Hd_(n_dis) - E_threshold_2 ! MPa*m ( E = U+W [MPa*km^2*10^3] )
!~           EndIf
!~         EndIf

!~ ****** INCREASING A0 if DeE<0

      If (DeE.le.0.D0) Then

        PROP_STEP = PROP_STEP - 1

        DEALLOCATE (Ut,Ub,x0z0,delta,Hd_)
      
        n_dis = n_dis_old
      
        ALLOCATE (Ut(n_dis),Ub(n_dis),x0z0(n_dis,2),delta(n_dis),Hd_(n_dis))

        x0z0  = x0z0_old
        delta = delta_old
        Hd_   = Hd_old
             
        A0 = A0 + DeA0*A0
        A  = A  + DeA0*A0*(A/A0)
        Call RESET_A0(A0)
      
        Call DEALLOC_Fold_Fbuf()
        Call SOLVE_CRACK(Ut,Ub,x0z0,delta,Hd_,A,n_dis,dyke_split,esc)
        Call COMP_ENERGY(Ut,Ub,x0z0,delta,Hd_,A,n_dis,E)

        x_top = x0z0(n_dis,1) - Hd_(n_dis)*sin(delta(n_dis))*0.5d0
        z_top = x0z0(n_dis,2) - Hd_(n_dis)*cos(delta(n_dis))*0.5d0
        
        ALLOCATE (OverP(n_dis))
        Call COMP_OverP(Ut,Ub,x0z0,delta,Hd_,A,n_dis,OverP,OverP_av,rho_iter)
        Write (24,*) iter, PROP_STEP, A0, A, rho_iter, DeE, n_dis, OverP(n_dis), x_top, z_top, 'N'
        DEALLOCATE (OverP)
!~ ****** END OF INCREASING A0 if DeE<0
      Else ! ****** if DeE>0 write some output

        Write (22,*) '--------------------------------------------------'
        Write (*,*)  '--------------------------------------------------'
        Write (22,*) 'n_dis =',n_dis
        Write (22,*) 'DeE  (iter=',iter, ') =',DeE
        Write (22,*) 'delta(     ',n_dis,') =',delta(n_dis)*180.D0/pi

        ALLOCATE (P_lit(n_dis))
        Call LITH_P(P_lit,x0z0,n_dis)
        Write (22,*) 'Plit (     ',n_dis,') =',P_lit(n_dis)
        Write (22,*) 'Plit (     ',n_dis-1,') =',P_lit(n_dis-1)
        Write (22,*) '[...]'
        Write (22,*) 'Plit (     ',2,') =',P_lit(2)
        Write (22,*) 'Plit (     ',1,') =',P_lit(1)
        DEALLOCATE (P_lit)

        Write (22,*) 'x0z0(1,1) = ',x0z0(1,1)
        Write (22,*) 'x0z0(1,2) = ',x0z0(1,2)

        Write (22,*) 'A_iter =',A
        Write (22,*) 'x_top  =',x_top
        Write (22,*) 'z_top  =',z_top
        Write (*,*)  'n_dis=',n_dis
        Write (*,*)  'delta(',n_dis,') =',delta(n_dis)*180.D0/pi
        Write (*,*)  'A(iter=',iter,') =',A
        Write (*,*)  'x_top =',x_top,'  z_top =',z_top

        ALLOCATE (OverP(n_dis))
        Call COMP_OverP(Ut,Ub,x0z0,delta,Hd_,A,n_dis,OverP,OverP_av,rho_iter)
        Write (22,*) 'OverP_av =', OverP_av
        Write (*,*)  'OverP_av =', OverP_av
        Write (24,*) iter, PROP_STEP, A0, A, rho_iter, DeE, n_dis, OverP(n_dis), x_top, z_top, 'Y'
        DEALLOCATE (OverP)

!~ DELTA -------------------------------------------------------------
!~         Call COMP_DELTA_iter(x_top,z_top,delta,Hd_,n_dis,z_start,       &
!~      &                           OverP_av,L_norm,                       &
!~      &                           dip_av,dip_ndis,dip_max,               &
!~      &                           Sig_tip,Sig_xz_av,Sig_xz_n,Sig_xz_max, &
!~      &                           DELTA1,DELTA2,DELTA3)
!~ -------------------------------------------------------------------

        Write (22,*) 'Ut_max =', maxval(Ut), 'Ut_last_el =', Ut(1)
        Ub_max = maxval(Ub)
        Ub_min = minval(Ub)
        If (abs(Ub_max).lt.abs(Ub_min)) Ub_max = Ub_min
        Write (22,*) 'Ub_max =', Ub_max, 'Ub_last_el =', Ub(1)

        If ( (esc.eqv..true.).or.(DeE.le.0.D0)) Then
          If (esc.eqv..true.) write(*,*) '***** TEST DISLOCATIONS CLOSED *****'
          If (DeE.le.0.D0)    write(*,*) '** ENERGY RELEASE UNDER THRESHOLD **'
            
          Call DISPL_FIELD(Ut_old,Ub_old,x0z0_old,delta_old,Hd_old,n_dis_old,run,iter)
          Call STRSS_FIELD(Ut_old,Ub_old,x0z0_old,delta_old,Hd_old,n_dis_old,run,iter)
          Call CRACK_SHAPE(Ut_old,Ub_old,x0z0_old,delta_old,Hd_old,n_dis_old,run,iter)

          DEALLOCATE (Ut_old,Ub_old,x0z0_old,delta_old,Hd_old)
          xf_top = x_top
          zf_top = z_top
          Call DEALLOC_Fold_Fbuf()
          EXIT
        EndIf

!        Write (21,21) x_top,z_top,DeE
        Write (23,21) x_top,z_top,DeE
!~ DELTA -------------------------------------------------------------
!~         Write (24,22) x_top,z_top,OverP_av,L_norm,                      &
!~        &              Sig_xz_av,Sig_xz_n,Sig_xz_max,                    &
!~        &              DELTA1,DELTA2,DELTA3
        
!~         Write (23,*) '--------------------------------'
!~         Write (23,*) '------- ITER = ',iter,' --------'
!~         Write (23,*) '--------------------------------'
!~         Write (23,*) 'x_top, z_top    = ',x_top,z_top
!~         Write (23,*) 'S_xx, S_zz, Sxz = ',Sig_tip(1,1),Sig_tip(2,2),Sig_tip(1,2)
!~         Write (23,*) 'dip_av          = ',dip_av
!~         Write (23,*) 'dip_ndis        = ',dip_ndis
!~         Write (23,*) 'dip_max         = ',dip_max
!~         Write (23,*) 'S_xz_av         = ',Sig_xz_av
!~         Write (23,*) 'S_xz_n          = ',Sig_xz_n
!~         Write (23,*) 'S_xz_max        = ',Sig_xz_max
!~         Write (23,*) 'OverP_av        = ',OverP_av
!~         Write (23,*) 'L_norm          = ',L_norm
!~         Write (23,*) 'DELTA1, DELTA2, DELTA3  = ',DELTA1,DELTA2,DELTA3
!~         Write (23,*) ''
!~ END DELTA ---------------------------------------------------------
!~ COMPUTE P_lit sig_N sig_T -----------------------------------------
        Call COMP_LITH_P (Plit,x0z0(n_dis,2))
        Call COMP_sigN_sigT (sig_N,sig_T,x0z0(n_dis,1),x0z0(n_dis,2),delta(n_dis)) 
        Write (25,21) x0z0(n_dis,1),x0z0(n_dis,2),Plit,sig_N,sig_T
!~ END COMPUTE P_lit sig_N sig_T -------------------------------------

        If (movie.eqv..true.) Call CRACK_SHAPE(Ut,Ub,x0z0,delta,Hd_,n_dis,run,iter)

        EndIf ! ***** ENDIF FOR THE "INCREASE A0" IF CONDITION
      
! ********* DYKE SPLIT **********
! ********* DYKE SPLIT **********
! ********* DYKE SPLIT **********
! ********* DYKE SPLIT **********
! ********* DYKE SPLIT **********
! ********* DYKE SPLIT **********
            If (dyke_split.eqv..true.) Then
              Do i=n_dis,1,-1
                If (Ut(i).le.0.D0) Then

                  A_split    = 0.D0
                  A_tot      = 0.D0
                  Do l=1,n_dis
                     If (Ut(l).gt.0.D0) A_tot = A_tot + Hd_(l)*Ut(l)*1.D-3
                  EndDo
                  
                  n_dis = n_dis - i

                  Do l=1,n_dis
                    If (Ut(l+i).gt.0.D0) A_split = A_split + Hd_(l+1)*Ut(l+i)*1.D-3
                    delta(l)  = delta(l+i)
                    Hd_(l)    = Hd_(l+i)
                    x0z0(l,:) = x0z0(l+i,:)
                  EndDo

                  Do l=1,i
                    delta(n_dis+l)  = 0.D0
                    Hd_(n_dis+l)    = 0.D0
                    x0z0(n_dis+l,:) = 0.D0
                    Ut(n_dis+l)     = 0.D0
                    Ub(n_dis+l)     = 0.D0
                  EndDo

                  A0 = A0*(A_split)/A_tot
                  A  = A_split

                  Call RESET_A0(A0)
                  Call DEALLOC_Fold_Fbuf()
                  Call SOLVE_CRACK(Ut,Ub,x0z0,delta,Hd_,A,n_dis,dyke_split,esc)
                  Call COMP_ENERGY(Ut,Ub,x0z0,delta,Hd_,A,n_dis,E)
                  EXIT
      
                EndIf
              EndDo
            EndIf
! ********* END DYKE SPLIT **********
! ********* END DYKE SPLIT **********
! ********* END DYKE SPLIT **********
! ********* END DYKE SPLIT **********
! ********* END DYKE SPLIT **********

        If ( ( z_top.le.z_stop     + Hd_(n_dis)*cos(delta(n_dis)) ) &
     &              .or.                                            &
     &       ( x_top.le.x_stop_min + Hd_(n_dis)*sin(delta(n_dis)) ) &
     &              .or.                                            &
     &       ( x_top.ge.x_stop_max - Hd_(n_dis)*sin(delta(n_dis)) ) &
     &              .or.                                            &
     &       ( z_top.le.z_snap     + Hd_(n_dis)*cos(delta(n_dis)) ) &
     &              .or.                                            &
     &       ( iter.eq.n_iter )                                     ) Then
!~     &       ( dyke_split.eqv..true. )                              &
!~      &              .or.                                            &
 
          Call DISPL_FIELD(Ut,Ub,x0z0,delta,Hd_,n_dis,run,iter)
          Call STRSS_FIELD(Ut,Ub,x0z0,delta,Hd_,n_dis,run,iter)
          If (movie.eqv..false.) Call CRACK_SHAPE(Ut,Ub,x0z0,delta,Hd_,n_dis,run,iter)

!~ ! ********* DYKE SPLIT **********
!~             If (dyke_split.eqv..true.) Then
!~               Do i=2,n_dis-1
!~                 If (Ut(i).le.0.D0) Then
!~                   write(*,*) '-----------------------------------------------------------------'
!~                   write(*,*) '-----------------------------------------------------------------'
!~                   write(*,*) Ut
!~                   write(*,*) '-----------------------------------------------------------------'
!~                   write(*,*) 'i =' , i
!~                   write(*,*) '-----------------------------------------------------------------'
!~                   write(*,*) 'n_dis =' , n_dis

!~                   A_split    = 0.D0

!~                   n_dis = n_dis - i
!~                   write(*,*) '-----------------------------------------------------------------'
!~                   write(*,*) 'n_dis =' , n_dis

!~                   Do l=1,n_dis
!~                     delta(l)  = delta(l+i)
!~                     Hd_(l)    = Hd_(l+i)
!~                     x0z0(l,:) = x0z0(l+i,:)
!~                     A_split   = A_split + Hd_(l)*Ut(l)*1.D-3
!~                   EndDo

!~                   Do l=1,i
!~                     delta_in(n_dis+l)  = 0.D0
!~                     Hd_in(n_dis+l)     = 0.D0
!~                     x0z0_in(n_dis+l,:) = 0.D0
!~                     Ut(n_dis+l)        = 0.D0
!~                     Ub(n_dis+l)        = 0.D0
!~                   EndDo

!~                   A0 = A0*(A_split)/A
!~                   A  = A_split

!~                   Call RESET_A0(A0)
!~                   Call DEALLOC_Fold_Fbuf()
!~                   Call SOLVE_CRACK(Ut,Ub,x0z0,delta,Hd_,A,n_dis,dyke_split,esc)
!~                   write(*,*) '-----------------------------------------------------------------'
!~                   write(*,*) 'n_dis =' , n_dis
!~                   Call COMP_ENERGY(Ut,Ub,x0z0,delta,Hd_,A,n_dis,E)
!~                   write(*,*) '-----------------------------------------------------------------'
!~                   write(*,*) '-----------------------------------------------------------------'
!~                   write(*,*) Ut
!~                   write(*,*) '-----------------------------------------------------------------'
!~                   read (*,*)
!~                   EXIT
      
!~                 EndIf
!~               EndDo
!~               DEALLOCATE (Ut_old,Ub_old,x0z0_old,delta_old,Hd_old)
!~               CYCLE
!~             EndIf
!~ ! ********* END DYKE SPLIT **********


          If (z_top.le.z_snap+Hd_(n_dis)*cos(delta(n_dis))) Then
            z_snap=z_snap-Dz_snap
            DEALLOCATE (Ut_old,Ub_old,x0z0_old,delta_old,Hd_old)
            CYCLE
          EndIf

          DEALLOCATE (Ut_old,Ub_old,x0z0_old,delta_old,Hd_old)
          xf_top = x_top
          zf_top = z_top
          Call DEALLOC_Fold_Fbuf()
          EXIT
        EndIf

        DEALLOCATE (Ut_old,Ub_old,x0z0_old,delta_old,Hd_old)

      EndDo

      DEALLOCATE (Ut,Ub,x0z0,delta,Hd_)
      xf_top = x_top
      zf_top = z_top
      Call DEALLOC_Fold_Fbuf()

20    Format (i2.2)
21    Format (2x,5d16.6)
22    Format (2x,10d16.6)

!~       Close(21)
      Close(22)
      Close(23)
      Close(24)
      Close(25)

      End Subroutine START_PROP

!~ ---------------------------------------------------------------------
!~ ---------------------------------------------------------------------

      Subroutine TEST_DIR(Ut,Ub,x0z0,delta,Hd_,A,n_dis,E,n_dir,alfa,dyke_split,esc)
!~ !test different directions of propagation and give as output the crack
!~ !opening in the direction that minimise the total energie of the system
!~ !Ut     INPUT/OUTPUT
!~ !Ub     INPUT/OUTPUT
!~ !x0z0   INPUT/OUTPUT
!~ !delta  INPUT/OUTPUT
!~ !Hd_    INPUT/OUTPUT
!~ !A      INPUT/OUTPUT
!~ !n_dis  INPUT/OUTPUT
!~ !E      OUTPUT
!~ !n_dir  INPUT
!~ !alfa   INPUT
!~ !dyke_split OUTPUT
!~ !esc        OUTPUT
      IMPLICIT NONE

      REAL(8) :: Ut(:),Ub(:),x0z0(:,:),delta(:),Hd_(:)
      REAL(8) :: A,E,alfa
      INTEGER :: n_dis,n_dir
      LOGICAL :: dyke_split,esc

      REAL(8), ALLOCATABLE :: Ut_dir(:),Ub_dir(:),x0z0_dir(:,:),delta_dir(:),Hd_dir(:)
      REAL(8) :: A_dir,E_min,angle
      REAL(8), ALLOCATABLE :: E_(:),E_tot(:)
      INTEGER :: n_dis_in,dir,i
      INTEGER :: n_dis_min,delta_n_dis
      INTEGER, ALLOCATABLE :: n_dis_(:)
      LOGICAL, ALLOCATABLE :: dyke_split_dir(:),esc_dir(:)
      LOGICAL :: redo_test_dir

      redo_test_dir = .false.

10    If (n_dir.eq.1) Then
      
      dir=1
      n_dis = n_dis+1

      delta(n_dis)  = delta(n_dis-1) - ((n_dir+1)*0.5D0-dir)*alfa
      Hd_(n_dis)    = Hd_(n_dis-1)
      x0z0(n_dis,1) = x0z0(n_dis-1,1) -                                  &
     &                          0.5D0*Hd_(n_dis-1)*sin(delta(n_dis-1)) - &
     &                          0.5D0*Hd_(n_dis)  *sin(delta(n_dis))
      x0z0(n_dis,2) = x0z0(n_dis-1,2) -                                  &
     &                          0.5D0*Hd_(n_dis-1)*cos(delta(n_dis-1)) - &
     &                          0.5D0*Hd_(n_dis)  *cos(delta(n_dis))
        
      esc = .false.

      Call SOLVE_CRACK(Ut,Ub,x0z0,delta,Hd_,A,n_dis,dyke_split,esc)
      Call COMP_ENERGY(Ut,Ub,x0z0,delta,Hd_,A,n_dis,E)
      
      Else

      Write (22,*) '--------------------------------------------------'
      Write (22,*) '       ** TEST DIRECTIONS OF PROPAGATION **       '
      Write (22,*) '--------------------------------------------------'
      Write (*,*)  '--------------------------------------------------'

      ALLOCATE (E_(n_dir),E_tot(n_dir),n_dis_(n_dir),dyke_split_dir(n_dir),esc_dir(n_dir))

      ALLOCATE (x0z0_dir(n_dis+1,2),delta_dir(n_dis+1),Hd_dir(n_dis+1), &
     &          Ut_dir(n_dis+1),Ub_dir(n_dis+1))
     
      Call SET_Fold_eq_Fbuf

!~       Call omp_set_num_threads(5) 
!~ !$omp parallel do private(A_dir,n_dis_,x0z0_dir,delta_dir,Hd_dir,esc_dir,Ut_dir,Ub_dir,E_,E_tot,angle)
      Do dir=1,n_dir

        A_dir     = A
        n_dis_(dir) = n_dis+1

        x0z0_dir(1:n_dis,:) = x0z0(1:n_dis,:)
        delta_dir(1:n_dis)  = delta(1:n_dis)
        Hd_dir(1:n_dis)     = Hd_(1:n_dis)

        delta_dir(n_dis+1)  = delta_dir(n_dis) - ((n_dir+1)*0.5D0-dir)*alfa
        Hd_dir(n_dis+1)     = Hd_dir(n_dis)
        x0z0_dir(n_dis+1,1) = x0z0_dir(n_dis,1) -                             &
     &                          0.5D0*Hd_dir(n_dis)  *sin(delta_dir(n_dis)) - &
     &                          0.5D0*Hd_dir(n_dis+1)*sin(delta_dir(n_dis+1))
        x0z0_dir(n_dis+1,2) = x0z0_dir(n_dis,2) -                             &
     &                          0.5D0*Hd_dir(n_dis)  *cos(delta_dir(n_dis)) - &
     &                          0.5D0*Hd_dir(n_dis+1)*cos(delta_dir(n_dis+1))
        
        esc_dir(dir) = .false.

        Call SOLVE_CRACK(Ut_dir,Ub_dir,x0z0_dir,delta_dir,Hd_dir,A_dir,n_dis_(dir),dyke_split_dir(dir),esc_dir(dir))
        If (esc_dir(dir).eqv..false.) Then
          Call COMP_ENERGY(Ut_dir,Ub_dir,x0z0_dir,delta_dir,Hd_dir,A_dir,n_dis_(dir),E_(dir))
        Else
          E_(dir) = E_old
        EndIf
        Call SET_Fbuf_eq_Fold

!~         If (int_dir.eqv..false.) Then
!~           If (x0z0_dir(n_dis_(dir),2)-0.5D0*Hd_dir(n_dis_(dir))*cos(delta_dir(n_dis_(dir))).gt.0.D0) Then
        E_tot(dir) = E_(dir) + E_threshold_1*Hd_dir(n_dis_(dir))
!~           Else
!~             E_tot(dir) = E_(dir) + E_threshold_2*Hd_dir(n_dis_(dir))
!~           EndIf
!~         Else
!~           E_tot(dir) = E_(dir) + E_threshold_int*Hd_dir(n_dis_(dir))
!~         EndIf

        angle=delta_dir(n_dis_(dir))*180.D0/pi
        Write(22,*)'Etot(',dir,'/',n_dir,')=',E_(dir),' delta=',angle,           &
     &             ' n_dis=',n_dis_(dir),' dyke_split_dir=',dyke_split_dir(dir), &
     &             ' esc_dir=',esc_dir(dir)
        Write(*,*) 'test direction (',dir,'/',n_dir,') =',delta_dir(n_dis_(dir))*180.D0/pi

      EndDo
!~ !$omp end parallel do

      DEALLOCATE (x0z0_dir,delta_dir,Hd_dir,Ut_dir,Ub_dir)

      E_min = E_tot(1)
      Do i=1,n_dir
        If (E_tot(i).le.E_min) Then
          E     = E_(i)
          E_min = E_tot(i)
          dir   = i
        EndIf
      EndDo

!~    IF THE NUMBER OF DISL. EL. CHANGES FOR DIFFERENT TEST DIRECTION THE ENERGY RELEASE CAN BE AFFECTED BY THAT,
!~    THEREFORE I RE-TEST THE DIRECTIONS USING THE MIN NUMBER OF DISL. EL. FOUND TESTING THE DIRECTIONS OF PROPAGATION
      If (minval(n_dis_).ne.maxval(n_dis_)) Then
      If ( ( (dir.ne.1).and.(dir.ne.n_dir).and.(n_dis_(dir-1).ne.n_dis_(dir+1)) ) &
     &                                    .or.                                    &
     &                        ( (dir.eq.1).and.(n_dis_(dir).ne.n_dis_(dir+1)) )   &
     &                                    .or.                                    &
     &                    ( (dir.eq.n_dir).and.(n_dis_(dir).ne.n_dis_(dir-1)) ) ) Then
     
        If (redo_test_dir.eqv..true.) Then
          Write (*,*) 'ERROR IN SUBROUTINE TEST_DIR - MODULE PROPAGATION'
          Stop
        EndIf
      
        redo_test_dir = .true.
      
        n_dis_min   = minval(n_dis_)
        delta_n_dis = (n_dis_old+1) - n_dis_min
      
        n_dis = n_dis_min-1
        x0z0(1:n_dis,:) = x0z0_old(delta_n_dis+1:n_dis_old,:)
        delta(1:n_dis)  = delta_old(delta_n_dis+1:n_dis_old)
        Hd_(1:n_dis)    = Hd_old(delta_n_dis+1:n_dis_old)
        
        Call DEALLOC_Fold_Fbuf()
        DEALLOCATE (E_,E_tot,n_dis_,dyke_split_dir,esc_dir)
        Call SOLVE_CRACK(Ut,Ub,x0z0,delta,Hd_,A,n_dis,dyke_split,esc)
        
        GoTo 10
        
      EndIf
      EndIf
      
      If (redo_test_dir.eqv..true.) Then
        n_dis           = n_dis_old
        x0z0(1:n_dis,:) = x0z0_old(1:n_dis,:)
        delta(1:n_dis)  = delta_old(1:n_dis)
        Hd_(1:n_dis)    = Hd_old(1:n_dis)
        Call DEALLOC_Fold_Fbuf()
      EndIf      
!~ ------------------------------------------------------------------------------------------      
      
      n_dis = n_dis+1

      delta(n_dis)  = delta(n_dis-1) - ((n_dir+1)*0.5D0-dir)*alfa
      Hd_(n_dis)    = Hd_(n_dis-1)
      x0z0(n_dis,1) = x0z0(n_dis-1,1) -                               &
     &                          0.5D0*Hd_(n_dis-1)*sin(delta(n_dis-1)) - &
     &                          0.5D0*Hd_(n_dis)  *sin(delta(n_dis))
      x0z0(n_dis,2) = x0z0(n_dis-1,2) -                               &
     &                          0.5D0*Hd_(n_dis-1)*cos(delta(n_dis-1)) - &
     &                          0.5D0*Hd_(n_dis)  *cos(delta(n_dis))
        
      esc = .false.

      Call SOLVE_CRACK(Ut,Ub,x0z0,delta,Hd_,A,n_dis,dyke_split,esc)

!~ ------------------------------------------------------------------------------------------      
      If ((redo_test_dir.eqv..true.).and.(esc.eqv..false.)) Then
        Call COMP_ENERGY(Ut,Ub,x0z0,delta,Hd_,A,n_dis,E)
      EndIf      
!~ ------------------------------------------------------------------------------------------      

      DEALLOCATE (E_,E_tot,n_dis_,dyke_split_dir,esc_dir)
      
      EndIf

      End Subroutine TEST_DIR

!~ ---------------------------------------------------------------------
!~ ---------------------------------------------------------------------

!~       Subroutine TEST_DIR_INT(Ut,Ub,x0z0,delta,Hd_,A,n_dis,E,n_dir,alfa,esc)
!~ !!~ test different directions of propagation and give as output the crack
!~ !!~ opening in the direction that minimise the total energie of the system
!~ !!~ Ut     INPUT/OUTPUT
!~ !!~ Ub     INPUT/OUTPUT
!~ !!~ x0z0   INPUT/OUTPUT
!~ !!~ delta  INPUT/OUTPUT
!~ !!~ Hd_    INPUT/OUTPUT
!~ !!~ A      INPUT/OUTPUT
!~ !!~ n_dis  INPUT/OUTPUT
!~ !!~ E      OUTPUT
!~ !!~ n_dir  INPUT
!~ !!~ alfa   INPUT
!~ !!~ esc    OUTPUT
!~       IMPLICIT NONE
!~ 
!~       REAL(8) :: Ut(:),Ub(:),x0z0(:,:),delta(:),Hd_(:)
!~       REAL(8) :: A_old,E_old,A,E,alfa
!~       INTEGER :: n_dis,n_dir
!~       LOGICAL :: esc
!~ 
!~       REAL(8), ALLOCATABLE :: Ut_dir(:),Ub_dir(:),x0z0_dir(:,:),delta_dir(:),Hd_dir(:)
!~       REAL(8), ALLOCATABLE :: x0z0_in(:,:),delta_in(:),Hd_in(:)
!~       REAL(8) :: A_dir,A_in,E_dir
!~       REAL(8) :: c_1,c_2,Ddelta,delta_1,E_,E_dir_,h_n_dir
!~       INTEGER :: n_dis_dir,n_dis_in,dir
!~       LOGICAL :: esc_dir
!~ 
!~       esc=.true.
!~ 
!~       ALLOCATE (x0z0_dir(n_dis+1,2),delta_dir(n_dis+1),Hd_dir(n_dis+1),Ut_dir(n_dis+1),Ub_dir(n_dis+1))
!~       ALLOCATE (x0z0_in(n_dis,2),   delta_in(n_dis),   Hd_in(n_dis))
!~ 
!~       A_in = A
!~       n_dis_in = n_dis
!~       x0z0_in(:,:) = x0z0(1:n_dis,:)
!~       delta_in(:)  = delta(1:n_dis)
!~       Hd_in(:)     = Hd_(1:n_dis)
!~ 
!~       h_n_dir = (n_dir-1)*0.5D0
!~ 
!~       If (tan(delta(n_dis)).gt.0.D0) Then
!~         Ddelta =-alfa*2.D0
!~         c_1 = delta(n_dis)
!~         c_2 = delta(n_dis)-h_n_dir*alfa*2.D0
!~       Else
!~         Ddelta = alfa*2.D0
!~         c_1 = delta(n_dis)+h_n_dir*alfa*2.D0
!~         c_2 = delta(n_dis)
!~       EndIf
!~       
!~       If (sin(delta(n_dis)).gt.0.D0) Then
!~         delta_1 = 90.D0*pi/180.D0
!~       Else
!~         delta_1 =-90.D0*pi/180.D0
!~       EndIf
!~ 
!~       Write (22,*) '--------------------------------------------------'
!~       Write (22,*) '** TEST INTERFACE DIRECTION **'
!~       dir = 0
!~       Do While (c_1.gt.c_2)
!~ 
!~         dir = dir+1
!~ 
!~         A_dir     = A_in
!~         n_dis_dir = n_dis_in+1
!~ 
!~         x0z0_dir(1:n_dis_in,:) = x0z0_in(:,:)
!~         delta_dir(1:n_dis_in)  = delta_in(:)
!~         Hd_dir(1:n_dis_in)     = Hd_in(:)
!~ 
!~         delta_dir(n_dis_dir)  = delta_1 + Ddelta*(dir-1)
!~         Hd_dir(n_dis_dir)     = Hd_dir(n_dis_dir-1)
!~         x0z0_dir(n_dis_dir,1) = x0z0_dir(n_dis_dir-1,1) -                               &
!~      &                          0.5D0*Hd_dir(n_dis_dir-1)*sin(delta_dir(n_dis_dir-1)) - &
!~      &                          0.5D0*Hd_dir(n_dis_dir)  *sin(delta_dir(n_dis_dir))
!~         x0z0_dir(n_dis_dir,2) = x0z0_dir(n_dis_dir-1,2) -                               &
!~      &                          0.5D0*Hd_dir(n_dis_dir-1)*cos(delta_dir(n_dis_dir-1)) - &
!~      &                          0.5D0*Hd_dir(n_dis_dir)  *cos(delta_dir(n_dis_dir))
!~         
!~         esc_dir = .false.
!~ 
!~         Write (22,*) '--------------------------------------------------'
!~         Write (*,*)  '--------------------------------------------------'
!~ 
!~         Call SOLVE_CRACK(Ut_dir,Ub_dir,x0z0_dir,delta_dir,Hd_dir,A_dir,n_dis_dir,esc_dir)
!~         Call COMP_ENERGY(Ut_dir,Ub_dir,x0z0_dir,delta_dir,Hd_dir,A_dir,n_dis_dir,E_dir)
!~ 
!~         If (dir.eq.1) Then
!~           E_dir_ = E_dir + E_threshold_int*Hd_dir(n_dis_dir)
!~         Else
!~           If (x0z0_dir(n_dis_dir,2)-0.5D0*Hd_dir(n_dis_dir)*cos(delta_dir(n_dis_dir)).gt.0.D0) Then
!~             E_dir_ = E_dir + E_threshold_1*Hd_dir(n_dis_dir)
!~           Else
!~             E_dir_ = E_dir + E_threshold_2*Hd_dir(n_dis_dir)
!~           EndIf
!~         EndIf
!~ 
!~         Write (22,*) 'Etot (',dir,'/',n_dir,')=',E_dir_,' @ delta=',delta_dir(n_dis_dir)*180.D0/pi
!~         Write (*,*)  'Etot (',dir,'/',n_dir,')=',E_dir_,' @ delta=',delta_dir(n_dis_dir)*180.D0/pi
!~ 
!~         If(((dir.eq.1).or.(E_dir_.lt.E_)).and.(esc_dir.eqv..false.)) Then
!~           Ut    = Ut_dir
!~           Ub    = Ub_dir
!~           x0z0  = x0z0_dir
!~           delta = delta_dir
!~           Hd_   = Hd_dir
!~           A     = A_dir
!~           n_dis = n_dis_dir
!~           E     = E_dir
!~           E_    = E_dir_
!~           esc   = esc_dir
!~         EndIf
!~ 
!~         If (tan(delta(n_dis)).gt.0.D0) Then
!~           c_1 = delta_dir(n_dis_dir)
!~         Else
!~           c_2 = delta_dir(n_dis_dir)
!~         EndIf
!~  
!~       EndDo
!~ 
!~       DEALLOCATE (Ut_dir,Ub_dir,x0z0_dir,delta_dir,Hd_dir)
!~       
!~       If (delta(n_dis).eq.delta_1) int_dir=.true.
!~ 
!~       End Subroutine TEST_DIR_INT

!~ ---------------------------------------------------------------------
!~ ---------------------------------------------------------------------

      Subroutine COMP_ENERGY(Ut,Ub,x0z0,delta,Hd_,A,n_dis,E)
!~ compute total energy of the system
!~ Ut    INPUT
!~ Ub    INPUT
!~ x0z0  INPUT
!~ delta INPUT
!~ Hd_   INPUT
!~ A     INPUT
!~ n_dis INPUT
!~ E     OUTPUT
      IMPLICIT NONE
      
      REAL(8) :: Ut(:),Ub(:),x0z0(:,:),delta(:),Hd_(:)
      REAL(8) :: A,E
      INTEGER :: n_dis

      REAL(8) :: U,W

      Call COMP_U(Ut,Ub,x0z0,delta,Hd_,A,n_dis,U)
      Call COMP_W(Ut,Ub,x0z0,delta,Hd_,A,n_dis,W)
      
      E = U + W ! [MPa*km^2*10^3] it will be devided by the disl lenght in km so that the energy release will be in MPa*m

      End Subroutine COMP_ENERGY

!~ ---------------------------------------------------------------------
!~ ---------------------------------------------------------------------

      End Module PROPAGATION
