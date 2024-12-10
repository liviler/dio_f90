!==============================================================================!
! MODULE Expectation Rotation                                                   !
!                                                                              !
! This module calculates the                                                   !
!                                                                              !
! List of routines and functions:                                              !
! - subroutine                                                                 !
!==============================================================================!
Module ExpectationRotation
use Constants, only: r64
use Globals, only: constraint
implicit none
real(r64) :: Jx2_block(2),Jy2_block(2),Jz2_block(2)
real(r64) :: IM_OA_term1,IM_OA_term2,IM_OA_term3,IM_OA_term4,IM_OA_term5,&
            IM_OA_term1_T,IM_OA_term2_T,IM_OA_term3_T,IM_OA_term4_T,IM_OA_term5_T,&
            IM_OA_term1_3,IM_OA_term2_3,IM_OA_term3_3,IM_OA_term4_3, &
            IM_OA_term1_3_T,IM_OA_term2_3_T,IM_OA_term3_3_T,IM_OA_term4_3_T, &
            IM_OA_term1_4,IM_OA_term2_4,IM_OA_term3_4,IM_OA_term4_4
contains
subroutine calculate_rotational_correction_energy_DIR
    !-----------------------------------------------------
    !   Erot = -<J^2>/(2I) 
    !-----------------------------------------------------
    use Globals,only: expectations,option
    real(r64),dimension(3) :: JSquare,Jxyz
    real(r64),dimension(3) :: IM
    real(r64) :: Erot!,ErotC
    integer :: case
    case = option%crankCase
    if(case==1) then 
        call JxyzSquare(JSquare)
        call inertia_moment_Belyaev(IM)
        ! Erot = -JSquare(1)/(2*IMB(1))-JSquare(2)/(2*IMB(2))-JSquare(3)/(2*IMB(3))
        Erot = -(JSquare(1)+JSquare(2)+JSquare(3))/(2*IM(1)+0.0001)
    else if (case==2) then
        call JxyzSquare(JSquare)
        call inertia_moment_Nilsson_Prior(IM)
        Erot = -(JSquare(1)+JSquare(2)+JSquare(3))/(2*IM(1)+0.0001)
    else if (case==3) then
        call Jxy_Odd_Mass(Jxyz)
        call JxyzSquare_Odd_Mass(JSquare)
        call inertia_moment_Odd_Mass(IM)
        Erot = -((JSquare(1)+JSquare(2)+JSquare(3)) - (Jxyz(1)**2+Jxyz(2)**2+Jxyz(3)**2))/(2*IM(1)+0.0001)
    endif
    expectations%Erot = Erot
    call write_Erot
    call store_rotational_energy_related_terms
    contains
    subroutine write_Erot
        use Globals,only: outputfile,constraint
        character(len=*), parameter ::  format1 = "(2(a5,2x),10(a9,5x))", &
                                        format2 = "(2(f5.2,2x),10(f12.8,2x))"
        if(constraint%index==1) then
            write(outputfile%u_rotationalE,format1) "beta2 ","beta3 ","Erot  ",&
                                                    "I_x  ","I_y  ","I_z  ","  <J^2_x>","  <J^2_y>","  <J^2_z>",&
                                                    "  <J_x> ","  <J_y> ","  <J_z> "
        endif
        write(outputfile%u_rotationalE,format2) constraint%betac(constraint%index),constraint%bet3c(constraint%index),&
              Erot,IM(1),IM(2),IM(3),JSquare(1),JSquare(2),JSquare(3),Jxyz(1),Jxyz(2),Jxyz(3)
    end subroutine
end subroutine

subroutine  store_rotational_energy_related_terms
    use Globals, only: constraint,expectations,outputfile
    real(r64),dimension(3) :: JSqure_EE ,JSqure_OA,Jxyz_OA
    real(r64),dimension(3) :: IM_Belyaev, IM_Nilsson, IM_OA,IM_OA_T,IM_OA_3,IM_OA_3_T,IM_OA_4,IM_OA_4_add
    character(len=*), parameter ::  format1 = "(2(a16,2x),4x,17(a10,4x),6(a12,4x),16(a13,4x))",&
                                    format2 = "(2(f5.2,14x),17(f12.8,2x),6(f12.8,4x),16(f12.8,5x))"
    integer :: start_clock, end_clock, clocks_per_second
    real :: total_time
 
    call JxyzSquare(JSqure_EE)
    call JxyzSquare_Odd_Mass(JSqure_OA)
    call Jxy_Odd_Mass(Jxyz_OA)
    ! Inertia
    call inertia_moment_Belyaev(IM_Belyaev)
    call inertia_moment_Nilsson_Prior(IM_Nilsson)
    call inertia_moment_Odd_Mass(IM_OA)
    call inertia_moment_Odd_Mass_T(IM_OA_T)
    call inertia_moment_Odd_Mass_3(IM_OA_3)
    call inertia_moment_Odd_Mass_3_T(IM_OA_3_T)
    call inertia_moment_Odd_Mass_4(IM_OA_4)
    call inertia_moment_Odd_Mass_4_add(IM_OA_4_add)
    if(constraint%index==1) then
        write(outputfile%u_rotationalTerms,format1) "beta2_constraint","beta2_calculate",&
             "<J^2_x>_EE","<J^2_y>_EE","<J^2_z>_EE","<J^2_x>_OA","<J^2_y>_OA","<J^2_z>_OA","<J_x>_OA","<J_y>_OA","<J_z>_OA",&
             "I_Belyaev","I_Nilsson","I_OA","I_OA_term1","I_OA_term2","I_OA_term3","I_OA_term4","I_OA_term5",&
             "I_OA_T","I_OA_term1_T","I_OA_term2_T","I_OA_term3_T","I_OA_term4_T","I_OA_term5_T",&
             "I_OA_3","I_OA_term1_3","I_OA_term2_3","I_OA_term3_3","I_OA_term4_3", &
             "I_OA_3_T","I_OA_term1_3_T","I_OA_term2_3_T","I_OA_term3_3_T","I_OA_term4_3_T", &
             "I_OA_4","I_OA_term1_4","I_OA_term2_4","I_OA_term3_4","I_OA_term4_4",&
             "I_OA_4_add"
    endif
    write(outputfile%u_rotationalTerms,format2) constraint%betac(constraint%index),expectations%beta2,&
            JSqure_EE,JSqure_OA,Jxyz_OA,&
            IM_Belyaev(1),IM_Nilsson(1),IM_OA(1),IM_OA_term1,IM_OA_term2,IM_OA_term3,IM_OA_term4,IM_OA_term5,&
            IM_OA_T(1),IM_OA_term1_T,IM_OA_term2_T,IM_OA_term3_T,IM_OA_term4_T,IM_OA_term5_T,&
            IM_OA_3(1),IM_OA_term1_3,IM_OA_term2_3,IM_OA_term3_3,IM_OA_term4_3,&
            IM_OA_3_T(1),IM_OA_term1_3_T,IM_OA_term2_3_T,IM_OA_term3_3_T,IM_OA_term4_3_T,&
            IM_OA_4(1),IM_OA_term1_4,IM_OA_term2_4,IM_OA_term3_4,IM_OA_term4_4,&
            IM_OA_4_add(1)
    
end subroutine

subroutine JxyzSquare(JSquare)
    !----------------------------------------------------------------------------
    !   This subroutine is to calculate the expectation value of J^2_x, J^2_y, J^2_z
    ! of the ground state.
    !
    !   The ground state is BCS wave function, and the single-particle wave function 
    ! expanded in the spherical harmonic oscillator basis.
    !
    !   Return:
    !       JSquare(1): <J^2_x>
    !       JSquare(2): <J^2_y>
    !       JSquare(3): <J^2_z>
    !----------------------------------------------------------------------------
    use Constants, only: zero
    use Globals, only: dirac,pairing
    real(r64),dimension(3) :: JSquare
    real(r64) :: Jx2, Jy2, Jz2,vk,vl,uk,ul
    integer :: it,ib,k1,k2,k,l
    real(r64), dimension(2,2) :: JxME,JzME
    complex(16), dimension(2,2) :: JyME
    Jx2 = zero ! <J^2_x>
    Jy2 = zero ! <J^2_y>
    Jz2 = zero ! <J^2_z>
    do it = 1,2
        ib = 1 ! no symmetry is considered
        k1 = 1
        k2 = dirac%nk(it) ! number of energies
        ! Jx2_block(it) = zero
        ! Jy2_block(it) = zero
        ! Jz2_block(it) = zero
        do k = k1,k2 ! k (different single particle state)
            if(pairing%vv(k,it) .lt. 1.d-10) cycle
            do l = k1,k2 ! k'
                call Jxyz_Matrix_Element(it,ib,k,l,JxME,JyME,JzME)
                vk = dsqrt(pairing%vv(k,it)/2)
                vl = dsqrt(pairing%vv(l,it)/2)
                uk = dsqrt(1.d0 - pairing%vv(k,it)/2)
                ul = dsqrt(1.d0 - pairing%vv(l,it)/2)
                Jx2 = Jx2 + (JxME(1,1)**2 + JxME(2,1)**2 + JxME(1,2)**2 + JxME(2,2)**2 )*vk**2*(1-vl**2) &
                    + 2*vk*uk*vl*ul*(JxME(1,1)*JxME(2,2) - JxME(1,2)*JxME(2,1))
                Jy2 = Jy2 + (-JyME(1,1)**2 -JyME(2,1)**2 -JyME(1,2)**2 -JyME(2,2)**2 )*vk**2*(1-vl**2) &
                    + 2*vk*uk*vl*ul*(JyME(1,1)*JyME(2,2) - JyME(1,2)*JyME(2,1))
                Jz2 = Jz2 + (JzME(1,1)**2 + JzME(2,1)**2 + JzME(1,2)**2 + JzME(2,2)**2 )*vk**2*(1-vl**2) &
                    + 2*vk*uk*vl*ul*(JzME(1,1)*JzME(2,2) - JzME(1,2)*JzME(2,1))
                ! if(k==pairing%block_level(it) .or. l==pairing%block_level(it)) then
                !     Jx2_block(it) = Jx2_block(it) + (JxME(1,1)**2 + JxME(2,1)**2 + JxME(1,2)**2 + JxME(2,2)**2 ) &
                !                     *vk**2*(1-vl**2) + 2*vk*uk*vl*ul*(JxME(1,1)*JxME(2,2) - JxME(1,2)*JxME(2,1))
                !     Jy2_block(it) = Jy2_block(it) + (-JyME(1,1)**2 -JyME(2,1)**2 - JyME(1,2)**2 - JyME(2,2)**2 ) &
                !                     *vk**2*(1-vl**2) + 2*vk*uk*vl*ul*(JyME(1,1)*JyME(2,2) - JyME(1,2)*JyME(2,1))
                !     Jz2_block(it) = Jz2_block(it) + (JzME(1,1)**2 + JzME(2,1)**2 + JzME(1,2)**2 + JzME(2,2)**2 ) &
                !                     *vk**2*(1-vl**2) + 2*vk*uk*vl*ul*(JzME(1,1)*JzME(2,2) - JzME(1,2)*JzME(2,1))
                ! endif
            enddo 
        enddo
    enddo
    JSquare(1) = Jx2
    JSquare(2) = dble(Jy2)
    JSquare(3) = Jz2

end subroutine

subroutine JxyzSquare_Odd_Mass(JSquare)
    !------------------------------------------------------------------------------------------
    !   This subroutine is to calculate the expectation value of J^2_x, J^2_y, J^2_z
    ! of the odd A nuclear wave function |OA>.
    !
    !   The odd A nuclear wave function |OA> is 
    !      |OA> = \alpha_{k_b}^{\dagger} |BCS>,
    ! and the single-particle wave function expanded in the spherical harmonic oscillator basis.
    !
    !   Return:
    !       JSquare(1): <OA|J^2_x|OA>
    !       JSquare(2): <OA|J^2_y|OA>
    !       JSquare(3): <OA|J^2_z|OA>
    !-----------------------------------------------------------------------------------------
    use Constants, only: zero
    use Globals, only: dirac,pairing,OddA
    real(r64),dimension(3) :: JSquare
    real(r64) :: Jx2, Jy2, Jz2,vk,vl,uk,ul,vk_block,uk_block,tmpx,tmpy,tmpz
    integer :: it,ib,k1,k2,k,l,k_block
    real(r64), dimension(2,2) :: JxME,JzME,JxME1,JzME1,JxME2,JzME2,JxME3,JzME3,JxME4,JzME4
    complex(16), dimension(2,2) :: JyME,JyME1,JyME2,JyME3,JyME4
    Jx2 = zero ! <OA|J^2_x|OA>
    Jy2 = zero ! <OA|J^2_y|OA>
    Jz2 = zero ! <OA|J^2_z|OA>
    do it = 1,2
        ib = 1 ! no symmetry is considered
        tmpx = zero
        tmpy = zero
        tmpz = zero
        if(it.eq.1) then
            k_block = OddA%qusiparticle_state_neutron !pairing%block_level(it)
            vk_block = dsqrt(pairing%vv(k_block,it)/2)
            uk_block = dsqrt(1.d0 - pairing%vv(k_block,it)/2)
            call Jxyz_Matrix_Element(it,ib,k_block,k_block,JxME1,JyME1,JzME1)
            ! Jx2 = 2*JxME1(1,2)*JxME1(2,1)*uk_block**2*vk_block**2 - (JxME1(1,1)*uk_block**2 )**2 - (JxME1(2,2)*vk_block**2)**2
            ! Jy2 = 2*JyME1(1,2)*JyME1(2,1)*uk_block**2*vk_block**2 - (JyME1(1,1)*uk_block**2 )**2 - (JyME1(2,2)*vk_block**2)**2
            ! Jz2 = 2*JzME1(1,2)*JzME1(2,1)*uk_block**2*vk_block**2 - (JzME1(1,1)*uk_block**2 )**2 - (JzME1(2,2)*vk_block**2)**2
            ! tmpx = tmpx + uk_block**2 * JxME1(1,1) - vk_block**2 * JxME1(2,2)
            ! tmpy = tmpy + uk_block**2 * JyME1(1,1) - vk_block**2 * JyME1(2,2)
            ! tmpz = tmpz + uk_block**2 * JzME1(1,1) - vk_block**2 * JzME1(2,2)
        else
            k_block = 0
        endif
        k1 = 1
        k2 = dirac%nk(it) ! number of energies
        do k = k1,k2 ! k (different single particle state)
            if(pairing%vv(k,it) .lt. 1.d-10) cycle
            vk = dsqrt(pairing%vv(k,it)/2)
            uk = dsqrt(1.d0 - pairing%vv(k,it)/2)
            do l = k1,k2 ! k'
                call Jxyz_Matrix_Element(it,ib,k,l,JxME,JyME,JzME)
                vl = dsqrt(pairing%vv(l,it)/2)
                ul = dsqrt(1.d0 - pairing%vv(l,it)/2)
                Jx2 = Jx2 + (JxME(1,1)**2 + JxME(2,1)**2 + JxME(1,2)**2 + JxME(2,2)**2 )*vk**2*(1-vl**2) &
                    + 2*vk*uk*vl*ul*(JxME(1,1)*JxME(2,2) - JxME(1,2)*JxME(2,1))
                Jy2 = Jy2 + (-JyME(1,1)**2 -JyME(2,1)**2 -JyME(1,2)**2 -JyME(2,2)**2 )*vk**2*(1-vl**2) &
                    + 2*vk*uk*vl*ul*(JyME(1,1)*JyME(2,2) - JyME(1,2)*JyME(2,1))
                Jz2 = Jz2 + (JzME(1,1)**2 + JzME(2,1)**2 + JzME(1,2)**2 + JzME(2,2)**2 )*vk**2*(1-vl**2) &
                    + 2*vk*uk*vl*ul*(JzME(1,1)*JzME(2,2) - JzME(1,2)*JzME(2,1))
            enddo 
            if(it==1) then
                call Jxyz_Matrix_Element(it,ib,k,k,JxME2,JyME2,JzME2)
                call Jxyz_Matrix_Element(it,ib,k,k_block,JxME3,JyME3,JzME3)
                call Jxyz_Matrix_Element(it,ib,k_block,k,JxME4,JyME4,JzME4)
                ! Jx2 = Jx2 + (JxME3(1,1)**2 + JxME3(2,1)**2)*uk_block**2*(1-2*vk**2) &
                !     - (JxME3(1,2)**2 + JxME3(2,2)**2 )*vk_block**2*(1-2*vk**2) &
                !     + (JxME2(1,1) + JxME2(2,2))*JxME1(1,1)*2*vk**2*uk_block**2 &
                !     - (JxME2(1,1) + JxME2(2,2))*JxME1(2,2)*2*vk**2*vk_block**2 &
                !     + (JxME3(1,2)*JxME3(2,1) + JxME4(1,2)*JxME4(2,1) - JxME3(1,1)*JxME3(2,2) - JxME4(1,1)*JxME4(2,2)) &
                !       *2*uk*vk*uk_block*vk_block
                Jx2 = Jx2 + (JxME3(1,1)*JxME4(1,1) + JxME3(2,1)*JxME4(1,2))*uk_block**2*(1-2*vk**2) &
                    - (JxME3(1,2)*JxME4(2,1) + JxME3(2,2)*JxME4(2,2))*vk_block**2*(1-2*vk**2) &
                    + (JxME2(1,1) + JxME2(2,2))*JxME1(1,1)*2*vk**2*uk_block**2 &
                    - (JxME2(1,1) + JxME2(2,2))*JxME1(2,2)*2*vk**2*vk_block**2 &
                    + (JxME3(1,2)*JxME3(2,1) + JxME4(1,2)*JxME4(2,1) - JxME3(1,1)*JxME3(2,2) - JxME4(1,1)*JxME4(2,2)) &
                      *2*uk*vk*uk_block*vk_block
                Jy2 = Jy2 + (JyME3(1,1)*JyME4(1,1) + JyME3(2,1)*JyME4(1,2))*uk_block**2*(1-2*vk**2) &
                    - (JyME3(1,2)*JyME4(2,1) + JyME3(2,2)*JyME4(2,2))*vk_block**2*(1-2*vk**2) &
                    + (JyME2(1,1) + JyME2(2,2))*JyME1(1,1)*2*vk**2*uk_block**2 &
                    - (JyME2(1,1) + JyME2(2,2))*JyME1(2,2)*2*vk**2*vk_block**2 &
                    + (JyME3(1,2)*JyME3(2,1) + JyME4(1,2)*JyME4(2,1) - JyME3(1,1)*JyME3(2,2) - JyME4(1,1)*JyME4(2,2)) &
                      *2*uk*vk*uk_block*vk_block
                Jz2 = Jz2 + (JzME3(1,1)*JzME4(1,1) + JzME3(2,1)*JzME4(1,2))*uk_block**2*(1-2*vk**2) &
                    - (JzME3(1,2)*JzME4(2,1) + JzME3(2,2)*JzME4(2,2))*vk_block**2*(1-2*vk**2) &
                    + (JzME2(1,1) + JzME2(2,2))*JzME1(1,1)*2*vk**2*uk_block**2 &
                    - (JzME2(1,1) + JzME2(2,2))*JzME1(2,2)*2*vk**2*vk_block**2 &
                    + (JzME3(1,2)*JzME3(2,1) + JzME4(1,2)*JzME4(2,1) - JzME3(1,1)*JzME3(2,2) - JzME4(1,1)*JzME4(2,2)) &
                      *2*uk*vk*uk_block*vk_block
                ! Jy2 = Jy2 + (-JyME2(1,1)**2 - JyME2(2,1)**2)*uk_block**2*(1-2*vk**2) &
                !     + (-JyME2(1,2)**2 - JyME2(2,2)**2 )*vk_block**2*(2*vk**2-1)
                ! Jz2 = Jz2 + (JzME2(1,1)**2 + JzME2(2,1)**2)*uk_block**2*(1-2*vk**2) &
                !     + (JzME2(1,2)**2 + JzME2(2,2)**2 )*vk_block**2*(2*vk**2-1)
                
                tmpx = tmpx + vk**2 * (JxME2(1,1) + JxME2(2,2))
                tmpy = tmpy + vk**2 * (JyME2(1,1) + JyME2(2,2))
                tmpz = tmpz + vk**2 * (JzME2(1,1) + JzME2(2,2))
            endif
        enddo
        Jx2 = Jx2 + tmpx**2
        Jy2 = Jy2 + tmpy**2
        Jz2 = Jz2 + tmpz**2
    enddo
    JSquare(1) = Jx2
    JSquare(2) = dble(Jy2)
    JSquare(3) = Jz2
end subroutine

subroutine Jxy_Odd_Mass(Jxyz)
    !------------------------------------------------------------------------------------------
    !   This subroutine is to calculate the expectation value of J_x, J_y, J_z
    ! of the odd A nuclear wave function |OA>.
    !
    !   The odd A nuclear wave function |OA> is 
    !      |OA> = \alpha_{k_b}^{\dagger} |BCS>,
    ! and the single-particle wave function expanded in the spherical harmonic oscillator basis.
    !
    !   Return:
    !       Jxyz(1): <OA|J_x|OA>
    !       Jxyz(2): <OA|J_y|OA>
    !       Jxyz(3): <OA|J_z|OA>
    !-----------------------------------------------------------------------------------------
    use Constants, only: zero
    use Globals, only: dirac,pairing,OddA
    real(r64),dimension(3) :: Jxyz
    integer :: it,ib,k,k1,k2,k_block
    real(r64) :: uk_block,vk_block,uk,vk,Jx,Jy,Jz
    real(r64), dimension(2,2) :: JxME1,JzME1,JxME2,JzME2
    complex(16), dimension(2,2) :: JyME1,JyME2
    Jx = zero
    Jy = zero
    Jz = zero
    do it = 1,2
        ib = 1
        if(it.eq.1) then
            k_block = OddA%qusiparticle_state_neutron !pairing%block_level(it)
            vk_block = dsqrt(pairing%vv(k_block,it)/2)
            uk_block = dsqrt(1.d0 - pairing%vv(k_block,it)/2)
            call Jxyz_Matrix_Element(it,ib,k_block,k_block,JxME1,JyME1,JzME1)
            Jx = Jx + uk_block**2 * JxME1(1,1) - vk_block**2 * JxME1(2,2)
            Jy = Jy + uk_block**2 * JyME1(1,1) - vk_block**2 * JyME1(2,2)
            Jz = Jz + uk_block**2 * JzME1(1,1) - vk_block**2 * JzME1(2,2)
        else
            k_block = 0
        endif
        k1 = 1
        k2 = dirac%nk(it) ! number of energies
        do k= k1,k2
            if(pairing%vv(k,it) .lt. 1.d-10) cycle
            vk = dsqrt(pairing%vv(k,it)/2)
            uk = dsqrt(1.d0 - pairing%vv(k,it)/2)
            call Jxyz_Matrix_Element(it,ib,k,k,JxME2,JyME2,JzME2)
            Jx = Jx + vk**2 * (JxME2(1,1) + JxME2(2,2))
            Jy = Jy + vk**2 * (JyME2(1,1) + JyME2(2,2))
            Jz = Jz + vk**2 * (JzME2(1,1) + JzME2(2,2))
        enddo
    enddo
    Jxyz(1) = Jx
    Jxyz(2) = dble(Jy)
    Jxyz(3) = Jz
end subroutine

subroutine Jxyz_Matrix_Element(it,ib,k,l,JxME,JyME,JzME)
    !------------------------------------------------------------------------------
    !  This subroutine is to calculate the matrix element of J_x, J_y, J_z. 
    !
    !  Input:
    !     k: state |k> ,and k > 0
    !     l: state |k'> and k > 0
    !     ib: the block of |k> and |k'> state
    !     it: 1 nuetron, 2 proton
    ! 
    !  Return:
    !     JxME : <k|J_x|k'>, <\bar{k}|J_x|k'>, <k|J_x|\bar{k'}>, <\bar{k}|J_x|\bar{k'}>
    !     JyME : <k|J_y|k'>, <\bar{k}|J_y|k'>, <k|J_y|\bar{k'}>, <\bar{k}|J_y|\bar{k'}>
    !     JzME : <k|J_z|k'>, <\bar{k}|J_z|k'>, <k|J_z|\bar{k'}>, <\bar{k}|J_z|\bar{k'}>
    !     where |\bar{k'}> is the time reversal state of |k>
    !--------------------------------------------------------------------------------
    use Constants, only: zero,half
    use Globals, only: BS,gfv
    integer, intent(in) :: k,l,ib,it
    real(r64), dimension(2,2), intent(out):: JxME,JzME
    complex(16), dimension(2,2), intent(out) :: JyME

    real(r64) :: Jx11,Jx12,Jz11,Jz12,hfnm1,hfnm2,hfnm2r,hfnj1,hfnj2,hfnj2r,temp1x,temp2x,temp1z,temp2z
    integer :: ifg,i0sp,ndsp,nfgsp,m1,m2,nr1,nl1,nj1,nm1,nr2,nl2,nj2,nm2,nr2r,nl2r,nj2r,nm2r
    complex(16) :: Jy11,Jy12,temp1y,temp2y
    Jx11 = zero ! <k|J_x|k'>
    Jx12 = zero ! <k|J_x|\bar{k'}>
    Jy11 = zero ! <k|J_y|k'>
    Jy12 = zero ! <k|J_y|\bar{k'}>
    Jz11 = zero ! <k|J_z|k'>
    Jz12 = zero ! <k|J_z|\bar{k'}>
    do ifg = 1,2 ! large and small components
        i0sp = BS%HO_sph%iasp(ib,ifg) ! beginning position of HO basis in ib block
        ndsp = BS%HO_sph%idsp(ib,ifg) ! dimesion of shperical HO basis in ib block
        nfgsp = (ifg-1)*BS%HO_sph%idsp(ib,1)
        do m1 = 1,ndsp ! m (different shperical HO state)
        if(BS%HO_sph%fg(nfgsp+m1,k,it) .eq. 0) cycle
        do m2 = 1,ndsp ! m'
            ! m
            nr1 = BS%HO_sph%nljm(i0sp+m1,1)
            nl1 = BS%HO_sph%nljm(i0sp+m1,2)
            nj1 = BS%HO_sph%nljm(i0sp+m1,3)
            nm1 = BS%HO_sph%nljm(i0sp+m1,4)
            ! m'
            nr2 = BS%HO_sph%nljm(i0sp+m2,1)
            nl2 = BS%HO_sph%nljm(i0sp+m2,2)
            nj2 = BS%HO_sph%nljm(i0sp+m2,3)
            nm2 = BS%HO_sph%nljm(i0sp+m2,4)
            ! time reversal states of m2
            nr2r = nr2
            nl2r = nl2
            nj2r = nj2
            nm2r = 1 - nm2 ! because we're adding 1/2 to all the m in the program, so (nm2r-1/2) = -(nm2-1/2) => nm2r = 1-nm2
            if(nr1.ne.nr2 .or. nl1.ne.nl2 .or. nj1.ne.nj2) cycle
            ! m - 1/2 (because we're adding 1/2 to all the m)
            hfnm1  = nm1 - 0.5
            hfnm2  = nm2 - 0.5
            hfnm2r = nm2r- 0.5
            hfnj1  = nj1 - 0.5
            hfnj2  = nj2 - 0.5
            hfnj2r = nj2r- 0.5

            temp1x = zero
            temp2x = zero
            temp1y = zero
            temp2y = zero
            temp1z = zero
            temp2z = zero
            ! component of <k|J_xyz|k'>
            if(nm1.eq.nm2+1) then
                temp1x = BS%HO_sph%fg(nfgsp+m1,k,it)*BS%HO_sph%fg(nfgsp+m2,l,it)*half*dsqrt((hfnj2-hfnm2)*(hfnj2+hfnm2+1))
                temp1y = BS%HO_sph%fg(nfgsp+m1,k,it)*BS%HO_sph%fg(nfgsp+m2,l,it)*DCMPLX(0.d0,-0.5d0) &
                         *dsqrt((hfnj2-hfnm2)*(hfnj2+hfnm2+1))
            endif
            if(nm1.eq.nm2-1) then
                temp1x = BS%HO_sph%fg(nfgsp+m1,k,it)*BS%HO_sph%fg(nfgsp+m2,l,it)*half*sqrt((hfnj2+hfnm2)*(hfnj2-hfnm2+1))
                temp1y = BS%HO_sph%fg(nfgsp+m1,k,it)*BS%HO_sph%fg(nfgsp+m2,l,it)*DCMPLX(0.d0,0.5d0) &
                         *sqrt((hfnj2+hfnm2)*(hfnj2-hfnm2+1))
            endif
            if(nm1.eq.nm2) then
                temp1z = BS%HO_sph%fg(nfgsp+m1,k,it)*BS%HO_sph%fg(nfgsp+m2,l,it)*hfnm1
            endif
            ! component of <k|J_xyz|\bar{k'}>
            if(nm1.eq.nm2r+1) then
                temp2x= gfv%iv(int(nl2+hfnj2-hfnm2))*gfv%iv(ifg-1)*BS%HO_sph%fg(nfgsp+m1,k,it)*BS%HO_sph%fg(nfgsp+m2,l,it) &
                        *half*dsqrt((hfnj2r-hfnm2r)*(hfnj2r+hfnm2r+1)) 
                temp2y= gfv%iv(int(nl2+hfnj2-hfnm2))*gfv%iv(ifg-1)*BS%HO_sph%fg(nfgsp+m1,k,it)*BS%HO_sph%fg(nfgsp+m2,l,it) &
                        *DCMPLX(0.d0,-0.5d0)*dsqrt((hfnj2r-hfnm2r)*(hfnj2r+hfnm2r+1)) 
            endif
            if(nm1.eq.nm2r-1) then
                temp2x= gfv%iv(int(nl2+hfnj2-hfnm2))*gfv%iv(ifg-1)*BS%HO_sph%fg(nfgsp+m1,k,it)*BS%HO_sph%fg(nfgsp+m2,l,it) &
                        *half*dsqrt((hfnj2r+hfnm2r)*(hfnj2r-hfnm2r+1))
                temp2y= gfv%iv(int(nl2+hfnj2-hfnm2))*gfv%iv(ifg-1)*BS%HO_sph%fg(nfgsp+m1,k,it)*BS%HO_sph%fg(nfgsp+m2,l,it) &
                        *DCMPLX(0.d0,0.5d0)*dsqrt((hfnj2r+hfnm2r)*(hfnj2r-hfnm2r+1))
            endif
            if(nm1.eq.nm2r) then
                temp2z= gfv%iv(int(nl2+hfnj2-hfnm2))*gfv%iv(ifg-1)*BS%HO_sph%fg(nfgsp+m1,k,it)*BS%HO_sph%fg(nfgsp+m2,l,it)*hfnm1
            endif
            Jx11 = Jx11 + temp1x
            Jx12 = Jx12 + temp2x
            Jy11 = Jy11 + temp1y
            Jy12 = Jy12 + temp2y
            Jz11 = Jz11 + temp1z
            Jz12 = Jz12 + temp2z
        enddo 
        enddo
    enddo
    ! JxME
    JxME(1,1) = Jx11  ! <k|J_x|k'>
    JxME(2,1) = Jx12  ! <\bar{k}|J_x|k'>
    JxME(1,2) = Jx12  ! <k|J_x|\bar{k'}>
    JxME(2,2) = -Jx11 ! <\bar{k}|J_x|\bar{k'}>
    !JyME
    JyME(1,1) = Jy11          ! <k|J_y|k'>
    JyME(2,1) = CONJG(Jy12)  ! <\bar{k}|J_y|k'>
    JyME(1,2) = Jy12          ! <k|J_y|\bar{k'}>
    JyME(2,2) = -CONJG(Jy11) ! <\bar{k}|J_y|\bar{k'}>
    !JzME
    JzME(1,1) = Jz11  ! <k|J_z|k'>
    JzME(2,1) = Jz12  ! <\bar{k}|J_z|k'>
    JzME(1,2) = Jz12  ! <k|J_z|\bar{k'}> 
    JzME(2,2) = -Jz11 ! <\bar{k}|J_z|\bar{k'}>
end subroutine 

subroutine inertia_moment_Belyaev(IMB)
    !-------------------------------------------------------------------------
    !   This subroutine is to calculate the moment of inertia by Belyaev formula
    !
    !   Return:
    !       IMB(1): moment of inertia in the x direction
    !       IMB(2): moment of inertia in the y direction
    !       IMB(3): moment of inertia in the z direction
    !-------------------------------------------------------------------------
    use Constants, only: zero
    use Globals, only: dirac,pairing
    real(r64),dimension(3) :: IMB
    real(r64), dimension(2,2) :: JxME,JzME
    complex(16), dimension(2,2) :: JyME
    real(r64) :: IMBx,IMBy,IMBz,vk,vl,uk,ul,Ek,El,fac,Jx11MS,Jy11MS,Jz11MS
    integer :: it,ib,k,l
    IMBx = zero
    IMBy = zero
    IMBz = zero
    do it=1,2
        ib = 1
        do k = 1,dirac%nk(it)
            do l = 1,dirac%nk(it)
                vk = dsqrt(pairing%vv(k,it)/2)
                vl = dsqrt(pairing%vv(l,it)/2)
                uk = dsqrt(1.d0 - pairing%vv(k,it)/2)
                ul = dsqrt(1.d0 - pairing%vv(l,it)/2)
                Ek = sqrt((dirac%ee(k,it)-pairing%ala(it))**2 + (pairing%skk(k,it)*pairing%de(k,it))**2)
                El = sqrt((dirac%ee(l,it)-pairing%ala(it))**2 + (pairing%skk(l,it)*pairing%de(l,it))**2)
                fac = 2*(uk*vl-ul*vk)**2/(Ek+El+1.d-8)
                call Jxyz_Matrix_Element(it,ib,k,l,JxME,JyME,JzME)
                Jx11MS = JxME(1,1)**2       ! |<k|Jx|l>|^2
                ! Jy11MS = abs(JyME(1,1))**2  ! |<k|Jy|l>|^2
                ! Jz11MS = JzME(1,1)**2       ! |<k|Jz|l>|^2
                IMBx = IMBx + fac*Jx11MS
                ! IMBy = IMBy + fac*Jy11MS
                ! IMBz = IMBz + fac*Jz11MS
            enddo
        enddo
    enddo
    IMB(1) = IMBx
    ! IMB(2) = IMBy
    ! IMB(3) = IMBz
end subroutine

subroutine inertia_moment_Nilsson_Prior(IM)
    !-------------------------------------------------------------------------
    !   This subroutine is to calculate the moment of inertia derived by Nilsson
    ! and Prior
    !
    !   Return:
    !       IM(1): moment of inertia in the x direction
    !       IM(2): moment of inertia in the y direction
    !       IM(3): moment of inertia in the z direction
    !-------------------------------------------------------------------------
    use Constants, only: zero
    use Globals, only: dirac,pairing
    real(r64),dimension(3) :: IM
    real(r64), dimension(2,2) :: JxME,JzME
    complex(16), dimension(2,2) :: JyME
    real(r64) :: IMx,IMy,IMz,vk,vl,uk,ul,Ek,El,fac,Jx11MS,Jy11MS,Jz11MS,Jx21MS,Jy21MS,Jz21MS
    integer :: it,ib,k,l
    IMx = zero
    IMy = zero
    IMz = zero
    do it=1,2
        ib = 1
        do k = 1,dirac%nk(it)
            do l = 1,dirac%nk(it)
                vk = dsqrt(pairing%vv(k,it)/2)
                vl = dsqrt(pairing%vv(l,it)/2)
                uk = dsqrt(1.d0 - pairing%vv(k,it)/2)
                ul = dsqrt(1.d0 - pairing%vv(l,it)/2)
                Ek = sqrt((dirac%ee(k,it)-pairing%ala(it))**2 + (pairing%skk(k,it)*pairing%de(k,it))**2)
                El = sqrt((dirac%ee(l,it)-pairing%ala(it))**2 + (pairing%skk(l,it)*pairing%de(l,it))**2)
                fac = 2*(uk*vl-ul*vk)**2/(Ek+El+1.d-8)
                call Jxyz_Matrix_Element(it,ib,k,l,JxME,JyME,JzME)
                Jx11MS = JxME(1,1)**2       ! |<k|Jx|l>|^2
                Jx21MS = JxME(2,1)**2       ! |<\bar{k}|Jx|l>|^2
                ! Jy11MS = abs(JyME(1,1))**2  ! |<k|Jy|l>|^2
                ! Jy21MS = abs(JyME(2,1))**2  ! |<\bar{k}|Jy|l>|^2
                ! Jz11MS = JzME(1,1)**2       ! |<k|Jz|l>|^2
                ! Jz21MS = JzME(2,1)**2       ! |<\bar{k}|Jz|l>|^2
                IMx = IMx + fac*(Jx11MS+Jx21MS)
                ! IMy = IMy + fac*(Jy11MS+Jy21MS)
                ! IMz = IMz + fac*(Jz11MS+Jz21MS)
            enddo
        enddo
    enddo
    IM(1) = IMx
    ! IM(2) = IMy
    ! IM(3) = IMz
end subroutine

subroutine inertia_moment_Odd_Mass(IM)
    !-----------------------------------------------------------------------------
    !   This subroutine is to calculate the moment of inertia for odd Mass nucleus
    !
    !   Return:
    !       IM(1): moment of inertia in the x direction
    !       IM(2): moment of inertia in the y direction
    !       IM(3): moment of inertia in the z direction
    !-----------------------------------------------------------------------------
    use Constants, only: zero
    use Globals, only: dirac,pairing,OddA
    real(r64),dimension(3) :: IM
    real(r64), dimension(2,2) :: JxME,JzME
    complex(16), dimension(2,2) :: JyME
    real(r64) :: IMx,IMy,IMz,vk,vl,uk,ul,Ek,El,fac1,Jx11MS,Jy11MS,Jz11MS,Jx21MS,Jy21MS,Jz21MS,&
                 vk_block,uk_block,Ek_block,fac2,fac3
    integer :: it,ib,k,l,k_block
    IMx = zero
    IMy = zero
    IMz = zero
    IM_OA_term1 = zero
    IM_OA_term2 = zero
    IM_OA_term3 = zero
    IM_OA_term4 = zero
    IM_OA_term5 = zero
    do it=1,2
        if(it==1) then
            k_block = OddA%qusiparticle_state_neutron !pairing%block_level(it)
            vk_block = dsqrt(pairing%vv(k_block,it)/2)
            uk_block = dsqrt(1.d0 - pairing%vv(k_block,it)/2)
            Ek_block = sqrt((dirac%ee(k_block,it)-pairing%ala(it))**2 + (pairing%skk(k_block,it)*pairing%de(k_block,it))**2)
        else
            k_block = 0
        endif
        ib = 1
        do k = 1,dirac%nk(it)
            if(k==k_block) cycle
            vk = dsqrt(pairing%vv(k,it)/2)
            uk = dsqrt(1.d0 - pairing%vv(k,it)/2)
            Ek = sqrt((dirac%ee(k,it)-pairing%ala(it))**2 + (pairing%skk(k,it)*pairing%de(k,it))**2)
            do l = 1,dirac%nk(it)
                if(l==k_block) cycle
                vl = dsqrt(pairing%vv(l,it)/2)
                ul = dsqrt(1.d0 - pairing%vv(l,it)/2)
                El = sqrt((dirac%ee(l,it)-pairing%ala(it))**2 + (pairing%skk(l,it)*pairing%de(l,it))**2)
                fac1 = 2*(uk*vl-ul*vk)**2/(Ek+El+1.d-8)
                call Jxyz_Matrix_Element(it,ib,k,l,JxME,JyME,JzME)
                Jx11MS = JxME(1,1)**2       ! |<k|Jx|l>|^2
                Jx21MS = JxME(2,1)**2       ! |<\bar{k}|Jx|l>|^2
                ! Jy11MS = abs(JyME(1,1))**2  ! |<k|Jy|l>|^2
                ! Jy21MS = abs(JyME(2,1))**2  ! |<\bar{k}|Jy|l>|^2
                ! Jz11MS = JzME(1,1)**2       ! |<k|Jz|l>|^2
                ! Jz21MS = JzME(2,1)**2       ! |<\bar{k}|Jz|l>|^2
                IMx = IMx + fac1*(Jx11MS+Jx21MS)
                ! IMy = IMy + fac*(Jy11MS+Jy21MS)
                ! IMz = IMz + fac*(Jz11MS+Jz21MS)
                IM_OA_term1 = IM_OA_term1 + fac1*Jx11MS
                IM_OA_term2 = IM_OA_term2 + fac1*Jx21MS
            enddo
            if(it .eq. 1) then
                ! if(Ek < Ek_block) 
                fac2 = 2*(uk*vk_block - uk_block*vk)**2/(Ek + Ek_block + 1.d-8) &
                +  2*(uk*uk_block + vk*vk_block)**2/(abs(Ek - Ek_block) + 1.d-8)
                fac3 = 2*(uk*vk_block - uk_block*vk)**2/(Ek + Ek_block + 1.d-8) &
                +  2*(uk*uk_block - vk*vk_block)**2/(abs(Ek - Ek_block) + 1.d-8)
                call Jxyz_Matrix_Element(it,ib,k,k_block,JxME,JyME,JzME)
                Jx11MS = JxME(1,1)**2       ! |<k|Jx|k_b>|^2
                Jx21MS = JxME(2,1)**2       ! |<\bar{k}|Jx|k_b>|^2
                IMx = IMx + fac2*Jx11MS + fac3*Jx21MS
                IM_OA_term3 = IM_OA_term3 + fac2*Jx11MS
                IM_OA_term4 = IM_OA_term4 + fac3*Jx21MS
                IM_OA_term5 = IM_OA_term5 +  2*(uk*uk_block + vk*vk_block)**2/(abs(Ek - Ek_block) + 1.d-8)*Jx11MS
                ! if(k==47) print *, 'term3_level','(',k,')',fac2*Jx11MS,fac3*Jx21MS,fac2,Jx11MS,fac3,Jx21MS
            endif
        enddo
    enddo

    IM(1) = IMx
    IM(2) = IMy
    IM(3) = IMz
end subroutine

subroutine inertia_moment_Odd_Mass_T(IM)
    !-----------------------------------------------------------------------------
    !   This subroutine is to calculate the moment of inertia for odd Mass nucleus
    !
    !   Return:
    !       IM(1): moment of inertia in the x direction
    !       IM(2): moment of inertia in the y direction
    !       IM(3): moment of inertia in the z direction
    !-----------------------------------------------------------------------------
    use Constants, only: zero
    use Globals, only: dirac,pairing,OddA
    real(r64),dimension(3) :: IM
    real(r64), dimension(2,2) :: JxME,JzME
    complex(16), dimension(2,2) :: JyME
    real(r64) :: IMx,IMy,IMz,vk,vl,uk,ul,Ek,El,fac1,Jx11MS,Jy11MS,Jz11MS,Jx21MS,Jy21MS,Jz21MS,&
                 vk_block,uk_block,Ek_block,fac2,fac3
    integer :: it,ib,k,l,k_block
    IMx = zero
    IMy = zero
    IMz = zero
    IM_OA_term1_T = zero
    IM_OA_term2_T = zero
    IM_OA_term3_T = zero
    IM_OA_term4_T = zero
    IM_OA_term5_T = zero
    do it=1,2
        if(it==1) then
            k_block = OddA%qusiparticle_state_neutron !pairing%block_level(it)
            vk_block = dsqrt(pairing%vv(k_block,it)/2)
            uk_block = dsqrt(1.d0 - pairing%vv(k_block,it)/2)
            Ek_block = sqrt((dirac%ee(k_block,it)-pairing%ala(it))**2 + (pairing%skk(k_block,it)*pairing%de(k_block,it))**2)
        else
            k_block = 0
        endif
        ib = 1
        do k = 1,dirac%nk(it)
            if(k==k_block) cycle
            vk = dsqrt(pairing%vv(k,it)/2)
            uk = dsqrt(1.d0 - pairing%vv(k,it)/2)
            Ek = sqrt((dirac%ee(k,it)-pairing%ala(it))**2 + (pairing%skk(k,it)*pairing%de(k,it))**2)
            do l = 1,dirac%nk(it)
                if(l==k_block) cycle
                vl = dsqrt(pairing%vv(l,it)/2)
                ul = dsqrt(1.d0 - pairing%vv(l,it)/2)
                El = sqrt((dirac%ee(l,it)-pairing%ala(it))**2 + (pairing%skk(l,it)*pairing%de(l,it))**2)
                fac1 = 2*(uk*vl-ul*vk)**2/(Ek+El+1.d-8)
                call Jxyz_Matrix_Element(it,ib,k,l,JxME,JyME,JzME)
                Jx11MS = JxME(1,1)**2       ! |<k|Jx|l>|^2
                Jx21MS = JxME(2,1)**2       ! |<\bar{k}|Jx|l>|^2
                ! Jy11MS = abs(JyME(1,1))**2  ! |<k|Jy|l>|^2
                ! Jy21MS = abs(JyME(2,1))**2  ! |<\bar{k}|Jy|l>|^2
                ! Jz11MS = JzME(1,1)**2       ! |<k|Jz|l>|^2
                ! Jz21MS = JzME(2,1)**2       ! |<\bar{k}|Jz|l>|^2
                IMx = IMx + fac1*(Jx11MS+Jx21MS)
                ! IMy = IMy + fac*(Jy11MS+Jy21MS)
                ! IMz = IMz + fac*(Jz11MS+Jz21MS)
                IM_OA_term1_T = IM_OA_term1_T + fac1*Jx11MS
                IM_OA_term2_T = IM_OA_term2_T + fac1*Jx21MS
            enddo
            if(it .eq. 1) then
                ! if(Ek < Ek_block) 
                fac2 = 2*(uk*vk_block - uk_block*vk)**2/(Ek + Ek_block + 1.d-8) &
                +  2*(uk*uk_block + vk*vk_block)**2/(Ek - Ek_block + 1.d-8)
                fac3 = 2*(uk*vk_block - uk_block*vk)**2/(Ek + Ek_block + 1.d-8) &
                +  2*(uk*uk_block - vk*vk_block)**2/(Ek - Ek_block + 1.d-8)
                call Jxyz_Matrix_Element(it,ib,k,k_block,JxME,JyME,JzME)
                Jx11MS = JxME(1,1)**2       ! |<k|Jx|k_b>|^2
                Jx21MS = JxME(2,1)**2       ! |<\bar{k}|Jx|k_b>|^2
                IMx = IMx + fac2*Jx11MS + fac3*Jx21MS
                IM_OA_term3_T = IM_OA_term3_T + fac2*Jx11MS
                IM_OA_term4_T = IM_OA_term4_T + fac3*Jx21MS
                IM_OA_term5_T = IM_OA_term5_T +  2*(uk*uk_block + vk*vk_block)**2/(Ek - Ek_block + 1.d-8)*Jx11MS
            endif
        enddo
    enddo

    IM(1) = IMx
    IM(2) = IMy
    IM(3) = IMz
end subroutine

subroutine inertia_moment_Odd_Mass_3(IM)
    !-----------------------------------------------------------------------------
    !   This subroutine is to calculate the moment of inertia for odd Mass nucleus
    !
    !   Return:
    !       IM(1): moment of inertia in the x direction
    !       IM(2): moment of inertia in the y direction
    !       IM(3): moment of inertia in the z direction
    !-----------------------------------------------------------------------------
    use Constants, only: zero
    use Globals, only: dirac,pairing,OddA
    real(r64),dimension(3) :: IM
    real(r64), dimension(2,2) :: JxME,JzME
    complex(16), dimension(2,2) :: JyME
    real(r64) :: IMx,IMy,IMz,vk,vl,uk,ul,Ek,El,fac1,Jx11MS,Jy11MS,Jz11MS,Jx21MS,Jy21MS,Jz21MS,&
                 vk_block,uk_block,Ek_block,fac2,fac3
    integer :: it,ib,k,l,k_block
    IMx = zero
    IMy = zero
    IMz = zero
    IM_OA_term1_3 = zero
    IM_OA_term2_3 = zero
    IM_OA_term3_3 = zero
    IM_OA_term4_3 = zero
    do it=1,2
        if(it==1) then
            k_block = OddA%qusiparticle_state_neutron !pairing%block_level(it)
            vk_block = dsqrt(pairing%vv(k_block,it)/2)
            uk_block = dsqrt(1.d0 - pairing%vv(k_block,it)/2)
            Ek_block = sqrt((dirac%ee(k_block,it)-pairing%ala(it))**2 + (pairing%skk(k_block,it)*pairing%de(k_block,it))**2)
        else
            k_block = 0
        endif
        ib = 1
        do k = 1,dirac%nk(it)
            if(k==k_block) cycle
            vk = dsqrt(pairing%vv(k,it)/2)
            uk = dsqrt(1.d0 - pairing%vv(k,it)/2)
            Ek = sqrt((dirac%ee(k,it)-pairing%ala(it))**2 + (pairing%skk(k,it)*pairing%de(k,it))**2)
            do l = 1,dirac%nk(it)
                ! if(l==k_block) cycle
                vl = dsqrt(pairing%vv(l,it)/2)
                ul = dsqrt(1.d0 - pairing%vv(l,it)/2)
                El = sqrt((dirac%ee(l,it)-pairing%ala(it))**2 + (pairing%skk(l,it)*pairing%de(l,it))**2)
                fac1 = 2*(uk*vl-ul*vk)**2/(Ek+El+1.d-8)
                call Jxyz_Matrix_Element(it,ib,k,l,JxME,JyME,JzME)
                Jx11MS = JxME(1,1)**2       ! |<k|Jx|l>|^2
                Jx21MS = JxME(2,1)**2       ! |<\bar{k}|Jx|l>|^2
                ! Jy11MS = abs(JyME(1,1))**2  ! |<k|Jy|l>|^2
                ! Jy21MS = abs(JyME(2,1))**2  ! |<\bar{k}|Jy|l>|^2
                ! Jz11MS = JzME(1,1)**2       ! |<k|Jz|l>|^2
                ! Jz21MS = JzME(2,1)**2       ! |<\bar{k}|Jz|l>|^2
                IMx = IMx + fac1*(Jx11MS+Jx21MS)
                ! IMy = IMy + fac*(Jy11MS+Jy21MS)
                ! IMz = IMz + fac*(Jz11MS+Jz21MS)
                IM_OA_term1_3 = IM_OA_term1_3 + fac1*Jx11MS
                IM_OA_term3_3 = IM_OA_term3_3 + fac1*Jx21MS
            enddo
            if(it .eq. 1) then
                ! if(Ek < Ek_block)

                ! fac2 = 2*(uk*uk_block + vk*vk_block)**2/(abs(Ek - Ek_block) + 0.1)
                ! fac3 = 2*(uk*uk_block - vk*vk_block)**2/(abs(Ek - Ek_block) + 0.1)

                fac2 = 2*(uk*uk_block + vk*vk_block)**2/(abs(Ek - Ek_block)+0.0001)
                ! fac3 = 2*(uk*uk_block - vk*vk_block)**2/(abs(Ek - Ek_block)+0.0001)
                call Jxyz_Matrix_Element(it,ib,k,k_block,JxME,JyME,JzME)
                Jx11MS = JxME(1,1)**2       ! |<k|Jx|k_b>|^2
                Jx21MS = JxME(2,1)**2       ! |<\bar{k}|Jx|k_b>|^2
                IMx = IMx + fac2*Jx11MS + fac2*Jx21MS
                IM_OA_term2_3 = IM_OA_term2_3 + fac2*Jx11MS
                IM_OA_term4_3 = IM_OA_term4_3 + fac2*Jx21MS
            endif
        enddo
    enddo

    IM(1) = IMx
    IM(2) = IMy
    IM(3) = IMz
end subroutine

subroutine inertia_moment_Odd_Mass_3_T(IM)
    !-----------------------------------------------------------------------------
    !   This subroutine is to calculate the moment of inertia for odd Mass nucleus
    !
    !   Return:
    !       IM(1): moment of inertia in the x direction
    !       IM(2): moment of inertia in the y direction
    !       IM(3): moment of inertia in the z direction
    !-----------------------------------------------------------------------------
    use Constants, only: zero
    use Globals, only: dirac,pairing,OddA
    real(r64),dimension(3) :: IM
    real(r64), dimension(2,2) :: JxME,JzME
    complex(16), dimension(2,2) :: JyME
    real(r64) :: IMx,IMy,IMz,vk,vl,uk,ul,Ek,El,fac1,Jx11MS,Jy11MS,Jz11MS,Jx21MS,Jy21MS,Jz21MS,&
                 vk_block,uk_block,Ek_block,fac2,fac3
    integer :: it,ib,k,l,k_block
    IMx = zero
    IMy = zero
    IMz = zero
    IM_OA_term1_3_T = zero
    IM_OA_term2_3_T = zero
    IM_OA_term3_3_T = zero
    IM_OA_term4_3_T = zero
    do it=1,2
        if(it==1) then
            k_block = OddA%qusiparticle_state_neutron !pairing%block_level(it)
            vk_block = dsqrt(pairing%vv(k_block,it)/2)
            uk_block = dsqrt(1.d0 - pairing%vv(k_block,it)/2)
            Ek_block = sqrt((dirac%ee(k_block,it)-pairing%ala(it))**2 + (pairing%skk(k_block,it)*pairing%de(k_block,it))**2)
        else
            k_block = 0
        endif
        ib = 1
        do k = 1,dirac%nk(it)
            if(k==k_block) cycle
            vk = dsqrt(pairing%vv(k,it)/2)
            uk = dsqrt(1.d0 - pairing%vv(k,it)/2)
            Ek = sqrt((dirac%ee(k,it)-pairing%ala(it))**2 + (pairing%skk(k,it)*pairing%de(k,it))**2)
            do l = 1,dirac%nk(it)
                ! if(l==k_block) cycle
                vl = dsqrt(pairing%vv(l,it)/2)
                ul = dsqrt(1.d0 - pairing%vv(l,it)/2)
                El = sqrt((dirac%ee(l,it)-pairing%ala(it))**2 + (pairing%skk(l,it)*pairing%de(l,it))**2)
                fac1 = 2*(uk*vl-ul*vk)**2/(Ek+El+1.d-8)
                call Jxyz_Matrix_Element(it,ib,k,l,JxME,JyME,JzME)
                Jx11MS = JxME(1,1)**2       ! |<k|Jx|l>|^2
                Jx21MS = JxME(2,1)**2       ! |<\bar{k}|Jx|l>|^2
                ! Jy11MS = abs(JyME(1,1))**2  ! |<k|Jy|l>|^2
                ! Jy21MS = abs(JyME(2,1))**2  ! |<\bar{k}|Jy|l>|^2
                ! Jz11MS = JzME(1,1)**2       ! |<k|Jz|l>|^2
                ! Jz21MS = JzME(2,1)**2       ! |<\bar{k}|Jz|l>|^2
                IMx = IMx + fac1*(Jx11MS+Jx21MS)
                ! IMy = IMy + fac*(Jy11MS+Jy21MS)
                ! IMz = IMz + fac*(Jz11MS+Jz21MS)
                IM_OA_term1_3_T = IM_OA_term1_3_T + fac1*Jx11MS
                IM_OA_term3_3_T = IM_OA_term3_3_T + fac1*Jx21MS
            enddo
            if(it .eq. 1) then
                ! if(Ek < Ek_block)
                fac2 = 2*(uk*uk_block + vk*vk_block)**2/(Ek - Ek_block + 1.d-8)
                ! fac3 = 2*(uk*uk_block - vk*vk_block)**2/(Ek - Ek_block + 0.0001)
                call Jxyz_Matrix_Element(it,ib,k,k_block,JxME,JyME,JzME)
                Jx11MS = JxME(1,1)**2       ! |<k|Jx|k_b>|^2
                Jx21MS = JxME(2,1)**2       ! |<\bar{k}|Jx|k_b>|^2
                IMx = IMx + fac2*Jx11MS + fac2*Jx21MS
                IM_OA_term2_3_T = IM_OA_term2_3_T + fac2*Jx11MS
                IM_OA_term4_3_T = IM_OA_term4_3_T + fac2*Jx21MS
            endif
        enddo
    enddo

    IM(1) = IMx
    IM(2) = IMy
    IM(3) = IMz
end subroutine

! subroutine inertia_moment_Odd_Mass_4(IM)
!     !-----------------------------------------------------------------------------
!     !   This subroutine is to calculate the moment of inertia for odd Mass nucleus
!     !
!     !   Return:
!     !       IM(1): moment of inertia in the x direction
!     !       IM(2): moment of inertia in the y direction
!     !       IM(3): moment of inertia in the z direction
!     !  Note:
!     !       set |E_k - E_kb| > 0.1
!     !-----------------------------------------------------------------------------
!     use Constants, only: zero
!     use Globals, only: dirac,pairing,OddA
!     real(r64),dimension(3) :: IM
!     real(r64), dimension(2,2) :: JxME,JzME
!     complex(16), dimension(2,2) :: JyME
!     real(r64) :: IMx,IMy,IMz,vk,vl,uk,ul,Ek,El,fac1,Jx11MS,Jy11MS,Jz11MS,Jx21MS,Jy21MS,Jz21MS,&
!                  vk_block,uk_block,Ek_block,fac2,fac3
!     integer :: it,ib,k,l,k_block
!     IMx = zero
!     IMy = zero
!     IMz = zero
!     IM_OA_term1_4 = zero
!     IM_OA_term2_4 = zero
!     IM_OA_term3_4 = zero
!     IM_OA_term4_4 = zero
!     write(444,*) 'beta2',constraint%betac(constraint%index)
!     do it=1,2
!         if(it==1) then
!             k_block = OddA%qusiparticle_state_neutron !pairing%block_level(it)
!             vk_block = dsqrt(pairing%vv(k_block,it)/2)
!             uk_block = dsqrt(1.d0 - pairing%vv(k_block,it)/2)
!             Ek_block = sqrt((dirac%ee(k_block,it)-pairing%ala(it))**2 + (pairing%skk(k_block,it)*pairing%de(k_block,it))**2)
!         else
!             k_block = 0
!         endif
!         ib = 1
!         do k = 1,dirac%nk(it)
!             if(k==k_block) cycle
!             vk = dsqrt(pairing%vv(k,it)/2)
!             uk = dsqrt(1.d0 - pairing%vv(k,it)/2)
!             Ek = sqrt((dirac%ee(k,it)-pairing%ala(it))**2 + (pairing%skk(k,it)*pairing%de(k,it))**2)
!             do l = 1,dirac%nk(it)
!                 ! if(l==k_block) cycle
!                 vl = dsqrt(pairing%vv(l,it)/2)
!                 ul = dsqrt(1.d0 - pairing%vv(l,it)/2)
!                 El = sqrt((dirac%ee(l,it)-pairing%ala(it))**2 + (pairing%skk(l,it)*pairing%de(l,it))**2)
!                 fac1 = 2*(uk*vl-ul*vk)**2/(Ek+El+1.d-8)
!                 call Jxyz_Matrix_Element(it,ib,k,l,JxME,JyME,JzME)
!                 Jx11MS = JxME(1,1)**2       ! |<k|Jx|l>|^2
!                 Jx21MS = JxME(2,1)**2       ! |<\bar{k}|Jx|l>|^2
!                 ! Jy11MS = abs(JyME(1,1))**2  ! |<k|Jy|l>|^2
!                 ! Jy21MS = abs(JyME(2,1))**2  ! |<\bar{k}|Jy|l>|^2
!                 ! Jz11MS = JzME(1,1)**2       ! |<k|Jz|l>|^2
!                 ! Jz21MS = JzME(2,1)**2       ! |<\bar{k}|Jz|l>|^2
!                 IMx = IMx + fac1*(Jx11MS+Jx21MS)
!                 ! IMy = IMy + fac*(Jy11MS+Jy21MS)
!                 ! IMz = IMz + fac*(Jz11MS+Jz21MS)
!                 IM_OA_term1_4 = IM_OA_term1_4 + fac1*Jx11MS
!                 IM_OA_term3_4 = IM_OA_term3_4 + fac1*Jx21MS
!             enddo
!             if(it .eq. 1) then
!                 if(abs(Ek - Ek_block) < 1 .and. abs(constraint%betac(constraint%index))> 0.2) cycle
!                 ! if(k==48 .or. k==49) cycle
!                 fac2 = 2*(uk*uk_block + vk*vk_block)**2/(abs(Ek - Ek_block)+1.d-8)
!                 ! fac3 = 2*(uk*uk_block - vk*vk_block)**2/(abs(Ek - Ek_block))
!                 call Jxyz_Matrix_Element(it,ib,k,k_block,JxME,JyME,JzME)
!                 Jx11MS = JxME(1,1)**2       ! |<k|Jx|k_b>|^2
!                 Jx21MS = JxME(2,1)**2       ! |<\bar{k}|Jx|k_b>|^2
!                 IMx = IMx + fac2*Jx11MS + fac2*Jx21MS
!                 IM_OA_term2_4 = IM_OA_term2_4 + fac2*Jx11MS
!                 IM_OA_term4_4 = IM_OA_term4_4 + fac2*Jx21MS
!                 write(444,"(2i8,6(f8.4))") k,k_block,fac2*Jx11MS,fac2,Jx11MS,Ek,Ek_block,IM_OA_term2_4
!             endif
!         enddo
!     enddo

!     IM(1) = IMx
!     IM(2) = IMy
!     IM(3) = IMz
! end subroutine

subroutine inertia_moment_Odd_Mass_4_add(IM)
    !-----------------------------------------------------------------------------
    !   This subroutine is to calculate the moment of inertia for odd Mass nucleus
    !
    !   Return:
    !       IM(1): moment of inertia in the x direction
    !       IM(2): moment of inertia in the y direction
    !       IM(3): moment of inertia in the z direction
    !  Note:
    !       set |E_k - E_kb| > 0.1
    !-----------------------------------------------------------------------------
    use Constants, only: zero
    use Globals, only: dirac,pairing,OddA
    real(r64),dimension(3) :: IM
    real(r64), dimension(2,2) :: JxME1,JzME1,JxME2,JzME2,JxME3,JzME3
    complex(16), dimension(2,2) :: JyME1,JyME2,JyME3
    real(r64) :: IMx,IMy,IMz,vk,uk,Ek,vk_block,uk_block,Ek_block,vk_dege,uk_dege,Ek_dege,tmp
    integer :: it,ib,k,l,k_block,k_dege
    IMx = zero
    IMy = zero
    IMz = zero

    do it=1,2
        if(it==1) then
            k_block = OddA%qusiparticle_state_neutron !pairing%block_level(it)
            vk_block = dsqrt(pairing%vv(k_block,it)/2)
            uk_block = dsqrt(1.d0 - pairing%vv(k_block,it)/2)
            Ek_block = sqrt((dirac%ee(k_block,it)-pairing%ala(it))**2 + (pairing%skk(k_block,it)*pairing%de(k_block,it))**2)
            
            k_dege = 49
            vk_dege = dsqrt(pairing%vv(k_dege,it)/2)
            uk_dege = dsqrt(1.d0 - pairing%vv(k_dege,it)/2)
            Ek_dege = sqrt((dirac%ee(k_dege,it)-pairing%ala(it))**2 + (pairing%skk(k_dege,it)*pairing%de(k_dege,it))**2)
            ib = 1
            do k = 1,dirac%nk(it)
                if(k==k_block .or. k==k_dege) cycle
                vk = dsqrt(pairing%vv(k,it)/2)
                uk = dsqrt(1.d0 - pairing%vv(k,it)/2)
                Ek = sqrt((dirac%ee(k,it)-pairing%ala(it))**2 + (pairing%skk(k,it)*pairing%de(k,it))**2)

                call Jxyz_Matrix_Element(it,ib,k_dege,k,JxME1,JyME1,JzME1)
                call Jxyz_Matrix_Element(it,ib,k,k_block,JxME2,JyME2,JzME2)
                call Jxyz_Matrix_Element(it,ib,k_block,k_dege,JxME3,JyME3,JzME3)
                tmp   = JxME1(1,1)*JxME2(1,1)*JxME3(1,1) &
                      *(uk_dege*uk+vk_dege*vk)*(uk*uk_block+vk*vk_block)*(uk_block*uk_dege+vk_block*vk_dege) &
                      + JxME1(1,2)*JxME2(2,1)*JxME3(1,1) &
                      *(uk_dege*uk-vk_dege*vk)*(uk*uk_block-vk*vk_block)*(uk_block*uk_dege+vk_block*vk_dege) &
                      + JxME1(2,1)*JxME2(1,1)*JxME3(1,2) &
                      *(uk_dege*uk-vk_dege*vk)*(uk*uk_block+vk*vk_block)*(uk_block*uk_dege-vk_block*vk_dege) &
                      + JxME1(1,1)*JxME2(2,1)*JxME3(1,2) &
                      *(uk_dege*uk+vk_dege*vk)*(uk*uk_block-vk*vk_block)*(uk_block*uk_dege-vk_block*vk_dege)
                IMx = IMx + 2*tmp/((Ek_block - Ek_dege)*(Ek_block-Ek))
            enddo
        endif
    enddo

    IM(1) = IMx
    IM(2) = IMy
    IM(3) = IMz
end subroutine

subroutine inertia_moment_Odd_Mass_4(IM)
    !-----------------------------------------------------------------------------
    !   This subroutine is to calculate the moment of inertia for odd Mass nucleus
    !
    !   Return:
    !       IM(1): moment of inertia in the x direction
    !       IM(2): moment of inertia in the y direction
    !       IM(3): moment of inertia in the z direction
    !  Note:
    !       set |E_k - E_kb| > 0.1
    !-----------------------------------------------------------------------------
    use Constants, only: zero
    use Globals, only: dirac,pairing,OddA
    real(r64),dimension(3) :: IM
    real(r64), dimension(2,2) :: JxME,JzME
    complex(16), dimension(2,2) :: JyME
    real(r64) :: IMx,IMy,IMz,vk,vl,uk,ul,Ek,El,fac1,Jx11MS,Jy11MS,Jz11MS,Jx21MS,Jy21MS,Jz21MS,&
                 vk_block,uk_block,Ek_block,fac2,fac3
    integer :: it,ib,k,l,k_block
    IMx = zero
    IMy = zero
    IMz = zero
    IM_OA_term1_4 = zero
    IM_OA_term2_4 = zero
    IM_OA_term3_4 = zero
    IM_OA_term4_4 = zero
    do it=1,2
        if(it==1) then
            k_block = OddA%qusiparticle_state_neutron !pairing%block_level(it)
            vk_block = dsqrt(pairing%vv(k_block,it)/2)
            uk_block = dsqrt(1.d0 - pairing%vv(k_block,it)/2)
            Ek_block = sqrt((dirac%ee(k_block,it)-pairing%ala(it))**2 + (pairing%skk(k_block,it)*pairing%de(k_block,it))**2)
        else
            k_block = 0
        endif
        ib = 1
        do k = 1,dirac%nk(it)
            if(k==k_block) cycle
            vk = dsqrt(pairing%vv(k,it)/2)
            uk = dsqrt(1.d0 - pairing%vv(k,it)/2)
            Ek = sqrt((dirac%ee(k,it)-pairing%ala(it))**2 + (pairing%skk(k,it)*pairing%de(k,it))**2)
            do l = 1,dirac%nk(it)
                ! if(l==k_block) cycle
                vl = dsqrt(pairing%vv(l,it)/2)
                ul = dsqrt(1.d0 - pairing%vv(l,it)/2)
                El = sqrt((dirac%ee(l,it)-pairing%ala(it))**2 + (pairing%skk(l,it)*pairing%de(l,it))**2)
                fac1 = 2*(uk*vl-ul*vk)**2/(Ek+El+1.d-8)
                call Jxyz_Matrix_Element(it,ib,k,l,JxME,JyME,JzME)
                Jx11MS = JxME(1,1)**2       ! |<k|Jx|l>|^2
                Jx21MS = JxME(2,1)**2       ! |<\bar{k}|Jx|l>|^2
                ! Jy11MS = abs(JyME(1,1))**2  ! |<k|Jy|l>|^2
                ! Jy21MS = abs(JyME(2,1))**2  ! |<\bar{k}|Jy|l>|^2
                ! Jz11MS = JzME(1,1)**2       ! |<k|Jz|l>|^2
                ! Jz21MS = JzME(2,1)**2       ! |<\bar{k}|Jz|l>|^2
                IMx = IMx + fac1*(Jx11MS+Jx21MS)
                ! IMy = IMy + fac*(Jy11MS+Jy21MS)
                ! IMz = IMz + fac*(Jz11MS+Jz21MS)
                IM_OA_term1_4 = IM_OA_term1_4 + fac1*Jx11MS
                IM_OA_term3_4 = IM_OA_term3_4 + fac1*Jx21MS
            enddo
            if(it .eq. 1) then
                fac2 = 2*(uk*uk_block + vk*vk_block)**2/(abs(Ek - Ek_block)+1.d-8)
                ! if(abs(constraint%betac(constraint%index))> 0.2) then
                !     ! fac2 = fac2 * 1/(1+2.71828**(-20*(Ek - Ek_block-1)))
                !     fac2 = fac2 * (1-2.71828**(-0.05*(Ek-Ek_block)**4))
                ! end if 
                fac2 = fac2 * regulator(Ek-Ek_block,constraint%betac(constraint%index))
                ! fac3 = 2*(uk*uk_block - vk*vk_block)**2/(abs(Ek - Ek_block))
                call Jxyz_Matrix_Element(it,ib,k,k_block,JxME,JyME,JzME)
                Jx11MS = JxME(1,1)**2       ! |<k|Jx|k_b>|^2
                Jx21MS = JxME(2,1)**2       ! |<\bar{k}|Jx|k_b>|^2
                IMx = IMx + fac2*Jx11MS + fac2*Jx21MS
                IM_OA_term2_4 = IM_OA_term2_4 + fac2*Jx11MS
                IM_OA_term4_4 = IM_OA_term4_4 + fac2*Jx21MS
            endif
        enddo
    enddo

    IM(1) = IMx
    IM(2) = IMy
    IM(3) = IMz
end subroutine

function regulator(x,beta2)
    real(r64) :: regulator
    real(r64) :: x ,beta2
    real(r64) :: k1,k2,beta2_0 
    real(r64) :: e = 2.71828
    k1 = 0.05
    k2 = 100
    beta2_0 = 0.2
    regulator = 1 - e**(-k1*x**4)/(1+e**(-k2*(abs(beta2) - beta2_0)))
end function

END Module ExpectationRotation