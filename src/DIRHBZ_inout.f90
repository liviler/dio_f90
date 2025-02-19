!==============================================================================!
! MODULE Inout                                                                 !
!                                                                              !
! This module contains functions and routines for reading and writing files.   !
!                                                                              !
! List of routines and functions:                                              !
! - subroutine                                                                 !
!==============================================================================!
MODULE Inout
use Globals, only: outputfile
use Constants, only: i8,i16,i64,r64,u_start,pi,ngh,ngl,nghl,itx,OUTPUT_PATH
implicit none
integer, private :: u_b23 = u_start + 1  ! the unit of b23.dat
integer, private :: u_dio = u_start + 2  ! the unit of dio.dat
integer, private :: u_wfs = u_start + 3  ! the unit of  ???
contains

subroutine read_file_b23
    use Globals, only: constraint
    character(len=*), parameter :: file_path_b23 = './input/b23.dat' ! the path of b23.dat
    integer(i64) :: fileb23_lines_max = 1000   ! maximum number of lines in file b23.dat

    integer :: length_count = 0 ! count the exact number of data rows in the 'b23.dat' file
    integer :: iostat ! store the io status of the 'read' function 
    integer(i64) :: index

    allocate(constraint%betac(fileb23_lines_max))
    allocate(constraint%bet3c(fileb23_lines_max))
    allocate(constraint%clam2(fileb23_lines_max))

    open(u_b23, file=file_path_b23, status='old')
    do index = 1, fileb23_lines_max+1
        if(length_count >= fileb23_lines_max) then
            write(*,*) "[Initialization]: 'fileb23_lines_max' too small"
            stop
        end if 
        read(u_b23,*,iostat=iostat) constraint%betac(index),constraint%bet3c(index),constraint%clam2(index)
        if (iostat /= 0) then
            if (iostat < 0) then
                ! read the end
                exit
            else
                ! error
                write(*, *) "Error reading from b23.dat!","iostat:",iostat
                stop
            end if
        end if
        length_count = length_count + 1
    end do
    constraint%length = length_count
    close(u_b23)
end subroutine read_file_b23

subroutine read_file_dio(ifPrint)
    use Globals, only: input_par,iteration,option,OddA
    logical,intent(in),optional :: ifPrint
    character(len=*), parameter :: file_path_dio = './input/dio.dat'

    character :: first_character
    real(r64) :: tmp
    character(len=*), parameter ::  format1= "(10x, i5)", &
                                    format2= "(10x, 2i3)", &
                                    format3= "(10x, f10.6)", &
                                    format4= "(10x, 2f10.6)", &
                                    format5= "(12x, a10)", &
                                    format6= "(a1, 9x, 2f10.3)", &
                                    format7= "(a2, i5)",&
                                    format8= "(10x, i8)"

    open(u_dio, file=file_path_dio, status='old')
    read(u_dio, format2) input_par%basis_n0f, input_par%basis_n0b    
    read(u_dio, format3) input_par%basis_b0
    read(u_dio, format6) first_character, tmp
    if(first_character == 'q') then
        input_par%basis_q = tmp
        input_par%basis_beta0 = dlog(input_par%basis_q)/(3*sqrt(5/(16*pi)))
    else
        input_par%basis_beta0  = tmp
        input_par%basis_q = exp(input_par%basis_beta0* 3*sqrt(5/(16*pi)))
    endif
    read(u_dio, format6) first_character, tmp
    if(first_character == 'q') then
        input_par%woodssaxon_qs = tmp
        input_par%woodssaxon_beta2 = dlog(input_par%woodssaxon_qs)/(3*sqrt(5/(16*pi)))
    else
        input_par%woodssaxon_beta2 = tmp
        input_par%woodssaxon_qs = exp(input_par%woodssaxon_beta2 * 3*sqrt(5/(16*pi)))
    endif
    read(u_dio, format3) input_par%woodssaxon_beta3
    read(u_dio, format1) input_par%iteration_max
    read(u_dio, format3) input_par%potential_mix
    read(u_dio, format1) input_par%inin
    read(u_dio, format7) input_par%nucleus_name, input_par%nucleus_mass_number
    read(u_dio, format2) input_par%pairing_ide
    read(u_dio, format4) input_par%pairing_dec
    read(u_dio, format4) input_par%pairing_ga
    read(u_dio, format4) input_par%pairing_del
    read(u_dio, format4) input_par%pairing_vpair
    read(u_dio, format5) input_par%force_name
    read(u_dio, format1) input_par%constraint_icstr
    read(u_dio, format3) input_par%constraint_cspr
    read(u_dio, format3) input_par%constraint_cmax
    read(u_dio, format1) input_par%option_iRHB
    read(u_dio, format1) input_par%option_iBlock
    read(u_dio, format8) input_par%block_level(1)
    read(u_dio, format8) input_par%block_level(2)
    read(u_dio, format8) input_par%K(1)
    read(u_dio, format1) input_par%Pi(1)
    read(u_dio, format8) input_par%K(2)
    read(u_dio, format1) input_par%Pi(2)
    read(u_dio, format1) input_par%option_blockMethod
    read(u_dio, format1) input_par%option_Erot
    close(u_dio)
    ! Iteration parameters is set here, while other parameters are set in their respective modules.
    iteration%iteration_max = input_par%iteration_max
    iteration%xmix = input_par%potential_mix
    !Option set here
    option%eqType = input_par%option_iRHB
    option%block_type =  input_par%option_iBlock
    option%block_method = input_par%option_blockMethod
    option%Erot_type = input_par%option_Erot
    if(ifPrint) call printInputConfig
    contains
    subroutine printInputConfig
        ! use Globals, only: outputfile
        ! character(len=*), parameter :: format11 = "(a,2i5)", &
        !                                format12 = "(a,2f12.6)", &
        !                                format13 = "(a,f8.5,5x,a,f8.5)", &
        !                                format14 = "(a,i6)", &
        !                                format15 = "(a,20x,a,a2,i4)"

        ! write(outputfile%u_config,*)'*************************BEGIN read_file_dio*************************'
        ! write(outputfile%u_config,format11) ' Number of oscillator shells : ',HO%n0f,HO%n0b
        ! write(outputfile%u_config,format12) ' Oscillator length b0 (fm)   : ',HO%b0  
        ! write(outputfile%u_config,format13) ' Basis deformation           :  beta0 =',HO%beta0,'q =',HO%q  
        ! write(outputfile%u_config,format13) ' Initial deformation         :  betas =',woodssaxon%beta2,'qs =',woodssaxon%qs       
        ! ! write(u_config,format12) ' Q1/Q2=                      : ',delta3
        ! write(outputfile%u_config,format14) ' Maximal number of iterations: ',iteration%iteration_max         
        ! write(outputfile%u_config,format12) ' Mixing parameter            : ',iteration%xmix             
        ! if (inin.eq.0) then                                                    
        !     write(outputfile%u_config,*) 'Initial wavefunctions       :  will read from tape'  
        ! else if (inin.eq.1) then                                                    
        !     write(outputfile%u_config,*) 'Initial wavefunctions       :  Saxon-Woods'  
        ! endif        
        ! write(outputfile%u_config,format15) ' Nucleus ',':     ',nucleus_attributes%name,nucleus_attributes%mass_number_int      
        ! write(outputfile%u_config,format11) ' Pairing control number      : ',pairing%ide
        ! write(outputfile%u_config,format12) ' Frozen Gap Parameters       : ',pairing%dec               
        ! write(outputfile%u_config,format12) ' Pairing-Constants           : ',pairing%ga                
        ! write(outputfile%u_config,format12) ' Initial values for Gap      : ',pairing%del               
        ! write(outputfile%u_config,format12) ' Vpair for delta force       : ',pairing%vpair
        ! write(outputfile%u_config,format14) ' Calculation with constraint : ',constraint%icstr
        ! ! write(u_config,format12) ' Constraint beta2-values     : ',constraint%betac
        ! ! write(u_config,format12) ' Constraint beta3-values     : ',constraint%bet3c
        ! write(outputfile%u_config,format12) ' Spring constant cspr        : ',constraint%cspr
        ! write(outputfile%u_config,format12) ' cutoff for dE/db            : ',constraint%cmax
        ! write(outputfile%u_config,"(a,/)")'*************************END read_file_dio****************************'
    end subroutine printInputConfig
end subroutine read_file_dio

subroutine set_output_filename(constraint_beta2,constraint_beta3)
    use Globals, only: BS,nucleus_attributes
    real(r64),intent(in):: constraint_beta2,constraint_beta3
    character :: sign_beta2, sign_beta3
    real(r64) :: abs2c, abs3c
    integer(i16), dimension(6) :: name
    integer(i16) :: name_nf1,name_nf2
    if (constraint_beta2 >= 0.d0) then
        sign_beta2 = '+'
    else
        sign_beta2 = '-'
    end if
    if (constraint_beta3 >= 0.d0) then
        sign_beta3 = '+'
    else
        sign_beta3 = '-'
    end if
    abs2c = abs(constraint_beta2)
    abs3c = abs(constraint_beta3)
    name(1) = abs2c + 48 !In ASCII, character '0' start from 48. 
    name(2) = mod(abs2c*10,10.d0)+48
    name(3) = mod(abs2c*100,10.d0)+48
    name(4) = abs3c+48
    name(5) = mod(abs3c*10,10.d0)+48
    name(6) = mod(abs3c*100,10.d0)+48
    name_nf1 = mod(BS%HO_cyl%n0f/10,10) + 48
    name_nf2 = mod(BS%HO_cyl%n0f,10) + 48
    ! the structure of ouput filenames are `dio_eMax`//HO%n0f//constraint_beta2*100//constraint_beta3*100//type
    ! like 'dio_eMax08+140+080.out' means HO%n0f=08, constraint_beta2= +1.40, constraint_beta3= +0.80, type is '.out'
    outputfile%outputf = OUTPUT_PATH//'dio'//'_eMax'//char(name_nf1)//char(name_nf2) &
                        //sign_beta2//char(name(1))//char(name(2))//char(name(3)) &
                        //'.'//sign_beta3//char(name(4))//char(name(5))//char(name(6))//'.out'
    outputfile%outputw = OUTPUT_PATH//'dio'//'_eMax'//char(name_nf1)//char(name_nf2) &
                        //sign_beta2//char(name(1))//char(name(2))//char(name(3)) &
                        //'.'//sign_beta3//char(name(4))//char(name(5))//char(name(6))//'.wel'
     outputfile%outdel = OUTPUT_PATH//'dio'//'_eMax'//char(name_nf1)//char(name_nf2) &
                        //sign_beta2//char(name(1))//char(name(2))//char(name(3)) &
                        //'.'//sign_beta3//char(name(4))//char(name(5))//char(name(6))//'.del'
     outputfile%outputwf=OUTPUT_PATH//nucleus_attributes%name//'_eMax'//char(name_nf1)//char(name_nf2) &
                        //'.dio'//sign_beta2//char(name(1))//char(name(2))//char(name(3)) &
                        //'.'//sign_beta3//char(name(4))//char(name(5))//char(name(6))//'.wf'
     outputfile%outputd= OUTPUT_PATH//'dio'//'_eMax'//char(name_nf1)//char(name_nf2) &
                        //sign_beta2//char(name(1))//char(name(2))//char(name(3)) &
                        //'.'//sign_beta3//char(name(4))//char(name(5))//char(name(6))//'.dens'
end subroutine set_output_filename

!--------------------------------------------------------------------------------------!
!subroutine read_fields:                                                               ! 
!reading of meson fields from the specified file                                       !
!--------------------------------------------------------------------------------------!
subroutine read_fields
    use Globals, only: expectations,pairing,constraint,fields
    character(len=500) :: WFS_DIR
    character(len=*),parameter :: format1 = "(1x,a2,8i4)", &
                                  format2 = "(2(6x,f12.6),6x,f12.8)",&
                                  format3 = "(2(6x,f12.6),6x,f20.12)",&
                                  format4 = "(10x,5f10.4)", &
                                  format5 = "(10x,5f12.6)", &
                                  format6 = "(4e20.12)"
    character(len=2) :: nucleus_name
    real(r64):: nucleus_mass_number,nucleus_proton_number,gauss_nh,basis_b0,basis_beta0,&
                iteration_si, force_amsig, force_amome,force_amrho,force_gsig,force_gome,&
                force_grho,force_g2,force_g3,tmp_vcn
    integer :: basis_NF,basis_NB, basis_block_number, basis_levels_number_fermions,&
                   index,it
    real(r64),dimension(nghl,itx) :: tmp_vps,tmp_vms,tmp_vpstot,tmp_vmstot


    WFS_DIR = find_file("GCM_FILES_DIR",'HFB.wfs')
    open(u_wfs,file=trim(WFS_DIR)//outputfile%outputw,status='old')
    read(u_wfs,format1)  nucleus_name,nucleus_mass_number,nucleus_proton_number,gauss_nh,&
                         basis_NF,basis_NB,basis_block_number,basis_levels_number_fermions
    read(u_wfs,format2) basis_b0,basis_beta0,iteration_si
    read(u_wfs,format3) expectations%ea, expectations%rms, expectations%betg
    read(u_wfs,format4) force_amsig, force_amome, force_amrho
    read(u_wfs,format4) force_gsig, force_gome, force_grho, force_g2, force_g3
    read(u_wfs,format5) pairing%ga,pairing%gg,pairing%pwi
    read(u_wfs,format5) pairing%del,pairing%dec
    read(u_wfs,format5) pairing%spk
    read(u_wfs,format5) pairing%ala!,fermi%tz
    read(u_wfs,format5) constraint%c1x, constraint%c2x, constraint%c3x
    read(u_wfs,format5) constraint%c1xold, constraint%c2xold, constraint%c3xold
    read(u_wfs,format5) constraint%calq1, constraint%calq2, constraint%calq3
    read(u_wfs,format6) fields%vps
    read(u_wfs,format6) fields%vms
    if (pairing%ide(1)==4) then
        read(u_wfs,*)
        read(u_wfs,format6) pairing%delq
    endif
    close(u_wfs)
    tmp_vps = reshape(fields%vps,[nghl,itx])
    tmp_vms = reshape(fields%vms,[nghl,itx])
    do index=1,nghl
        tmp_vcn = constraint%c1x * constraint%vc(index,1) +&
                  constraint%c2x * constraint%vc(index,2) + &
                  constraint%c3x * constraint%vc(index,3)
        do it=1,itx
            tmp_vpstot(index,it) = tmp_vps(index,it) + tmp_vcn
            tmp_vmstot(index,it) = tmp_vms(index,it) + tmp_vcn
        enddo
    enddo
    fields%vpstot = reshape(tmp_vpstot,[ngh,ngl,itx])
    fields%vmstot = reshape(tmp_vmstot,[ngh,ngl,itx])
end subroutine 


subroutine write_fields
    use Globals, only: nucleus_attributes,gauss,BS,iteration,force,expectations,pairing,constraint,fields
    character(len=500) :: WFS_DIR
    character(len=*),parameter :: format1 = "(1x,a2,8i4)", &
                                  format2 = "(3(a,f12.6),a,f12.8)",&
                                  format3 = "(2(a,f12.6),a,f20.12)",&
                                  format4 = "(10x,5f10.4)", &
                                  format5 = "(a,5f12.6)", &
                                  format6 = "(4e20.12)"
    character(len=2) :: nucleus_name
    real(r64):: nucleus_mass_number,nucleus_proton_number,gauss_nh,basis_b0,basis_beta0,&
                iteration_si, force_amsig, force_amome,force_amrho,force_gsig,force_gome,&
                force_grho,force_g2,force_g3,tmp_vcn
    integer :: basis_NF,basis_NB, basis_block_number, basis_levels_number_fermions,&
                   index,it
    real(r64),dimension(nghl,itx) :: tmp_vps,tmp_vms,tmp_vpstot,tmp_vmstot

    open(outputfile%u_outputw,file=outputfile%outputw,status='unknown')
    write(outputfile%u_outputw,format1) nucleus_attributes%name,nucleus_attributes%mass_number_int,&
                    int(nucleus_attributes%neutron_number),int(nucleus_attributes%proton_number),&
                    gauss%nh,BS%HO_cyl%n0f,BS%HO_cyl%n0b,BS%HO_cyl%nb,BS%HO_cyl%nt
    write(outputfile%u_outputw,format2) 'b0 =  ',BS%HO_cyl%b0,' bet0=',BS%HO_cyl%beta0,' si = ',iteration%si
    write(outputfile%u_outputw,format3) 'E/A = ',expectations%ea, ' RMS =',expectations%rms, '  Q = ',expectations%betg
    write(outputfile%u_outputw,format4) force%masses%amsig, force%masses%amome, force%masses%amrho
    write(outputfile%u_outputw,format4) force%couplm%gsig, force%couplm%gome, force%couplm%grho, force%nonlin%g2, force%nonlin%g3
    write(outputfile%u_outputw,format5) 'Pairing:  ',pairing%ga,pairing%gg,pairing%pwi
    write(outputfile%u_outputw,format5) 'Delta:    ',pairing%del,pairing%dec
    write(outputfile%u_outputw,format5) 'Spk:      ',pairing%spk
    write(outputfile%u_outputw,format5) 'Lambda:   ',pairing%ala,nucleus_attributes%neutron_number,nucleus_attributes%proton_number
    write(outputfile%u_outputw,format5) 'dE/db:    ',constraint%c1x, constraint%c2x, constraint%c3x
    write(outputfile%u_outputw,format5) 'cxold:    ',constraint%c1xold, constraint%c2xold, constraint%c3xold
    write(outputfile%u_outputw,format5) '<Q>:      ',constraint%calq1, constraint%calq2, constraint%calq3
    write(outputfile%u_outputw,format6) fields%vps
    write(outputfile%u_outputw,format6) fields%vms
    if (pairing%ide(1)==4) then
        write(outputfile%u_outputw,*) 'Delq:'
        write(outputfile%u_outputw,format6) pairing%delq
    endif
    close(outputfile%u_outputw)
end subroutine 

!--------------------------------------------------------------------------------------!
!function find_file:                                                                   ! 
!returns directory that find file_name from input environment variable(directory_name).!
!--------------------------------------------------------------------------------------!
function find_file(environment_name,file_name)
    character(500) :: find_file
    character(*) :: file_name,environment_name
    character(1000) :: path
    integer :: istart,iend,i
    logical :: isthere,last_chance
    call GETENV(trim(environment_name),path)
        path=adjustl(path)
        istart = 1
        do while (.true.)
            i = istart
            do while (path(i:i).ne.':')
                iend = i
                i = i+ 1
                if (path(i:i) == ' ') then
                    i = 0
                exit
                end if
            end do
            inquire(file=path(istart:iend)//'/'//trim(file_name),exist=isthere)
            if (isthere) then
                find_file = path(istart:iend)//'/'
                return
            else if (i == 0) then
                inquire(file='./'//file_name,exist=last_chance)
                if (last_chance) then
                    find_file = './'
                    return
                else
                    print*,'FILE NOT FOUND: ',trim(adjustl(file_name))
                    print*,'CHECK ENVIRONMENT VARIABLE: ',trim(adjustl(environment_name))
                    stop
                end if
            end if
            istart = iend+2
        end do
end function find_file


subroutine write_result_DIR
    use Expectation, only: calculate_expectation_DIR
    call calculate_expectation_DIR(.True.)
    call print_single_particle_energies
    call write_expectation
    call write_wavefuntion
    call write_fields
    call write_densit
end subroutine

subroutine print_single_particle_energies
    !-----------------------------------------
    !  write single particle energies 
    !----------------------------------------
    use Globals, only: outputfile,dirac,BS,pairing,constraint,expectations
    use Constants, only: nkx,tit,tp
    integer :: it,ib,kk,ip,nf,i0f,k1,k2,k,n,imax,min_qusiparticle_level
    real(r64) :: e0,smax,s,min_qusiparticle_Energy,block_qusiparticle_Energy
    real(r64),dimension(nkx) :: sln
    character(len=*), parameter ::  format1 = "(//,' Single-particle Energies ',a,/,1x,33(1h-))", &
                                    format2 = "(i3,' K =',i2,'/2',a1,'  ',a8,f5.1,2f10.3,f8.3,f8.3)",&
                                    format3 = "(2(a16,2x), 4(a17,2x))",&
                                    format4 = "(2(f5.2,14x), 2(3x,i4,8x,8x,f8.2,8x))"

    do it = 1,2
        write(outputfile%u_outputf,format1) tit(it)
        min_qusiparticle_Energy = sqrt((dirac%ee(1,it)-pairing%ala(it))**2+(pairing%skk(1,it)*pairing%de(1,it))**2)
        min_qusiparticle_level = 0
        e0 = dirac%ee(1,it)
        do ib = 1,BS%HO_cyl%nb
            kk = BS%HO_cyl%ikb(ib) ! K+1/2
            ip = BS%HO_cyl%ipb(ib)
            nf = BS%HO_cyl%id(ib,1)
            ! ng = BS%HO_cyl%id(ib,2)
            i0f = BS%HO_cyl%ia(ib,1)
            ! nh = nf + ng
            k1 = dirac%ka(ib,it) + 1
            k2 = dirac%ka(ib,it) + dirac%kd(ib,it)
            ! s0=0.d0
            do k = k1,k2
                if (dirac%ee(k,it)-pairing%ala(it).gt.30.0) cycle
                ! search for main oscillator component
                smax = 0.d0
                do n = 1,nf
                    s = abs(dirac%fg(n,k,it))
                    if (s.gt.smax) then
                        smax = s
                        imax = n
                    endif
                enddo
                ! printing
	            sln(k)=sqrt((dirac%ee(k,it)-pairing%ala(it))**2+(pairing%skk(k,it)*pairing%de(k,it))**2) !+ pairing%ala(it)
                if(sln(k) < min_qusiparticle_Energy ) then
                    min_qusiparticle_Energy = sln(k)
                    min_qusiparticle_level = k
                endif
                write(outputfile%u_outputf,format2) k, 2*kk-1,' ',BS%HO_cyl%tb(i0f+imax)(:2)//BS%HO_cyl%tb(i0f+imax)(4:),smax,&
                                                    dirac%ee(k,it),dirac%ee(k,it)-e0,pairing%vv(k,it)/2,sln(k)
	            ! if(it.eq.1.and.sln(k).lt.s0) then
                !     ibk(kk) = k
                !     s0 = sln(k)
	            ! end if
                if(it.eq.1.and.pairing%vv(k,it)/2.eq.0.5d0) then
                    pairing%ibk(kk) =  k  
                end if  
            enddo
	    enddo   
    enddo 

end subroutine

subroutine write_expectation
    use Globals,only: outputfile,constraint,iteration,expectations
    character(len=*), parameter ::  format1 = "(a,2x,a5,4x,2(a5,2x),2(a15,2x),3(a9,5x))", &
                                    format2 = "(i3,3x,f12.6,2x,2(f5.2,2x),2(6x,f5.2,6x),3(f12.6,2x))"
    if(constraint%index==1) then
        write(outputfile%u_outExpectation,format1) "iteration","  epsi","beta2 ","beta3 ",&
                                                    "beta2_calculate","beta3_calculate","Etot","Erot","Etot-Erot"
    endif
    write(outputfile%u_outExpectation,format2) iteration%ii,iteration%si,&
            constraint%betac(constraint%index),constraint%bet3c(constraint%index), &
            expectations%beta2,expectations%beta3,&
            expectations%etot,expectations%Erot,expectations%etot-expectations%Erot
end subroutine

subroutine write_wavefuntion
    use Globals, only: outputfile,dirac,pairing,BS,expectations
    open(outputfile%u_outputwf,file=outputfile%outputwf,form='unformatted',status='unknown')      
    write(outputfile%u_outputwf) dirac%ka
    write(outputfile%u_outputwf) dirac%kd
    write(outputfile%u_outputwf) pairing%vv*0.5d0
    write(outputfile%u_outputwf) pairing%vv 
    write(outputfile%u_outputwf) pairing%skk
    write(outputfile%u_outputwf) BS%HO_sph%fg  
    write(outputfile%u_outputwf) BS%HO_cyl%nb
    write(outputfile%u_outputwf) expectations%ecm 
    write(outputfile%u_outputwf) pairing%ibk
    close(outputfile%u_outputwf)
end subroutine

subroutine write_densit_old
    use Globals, only: gauss,densities,BS
    integer :: ih,il,ihl
    character(len=*), parameter ::  format1 = "(4(a12,2x))", &
                                    format2 = "(4(f12.6,2x))"

    open(unit=23,file='dens_neut.out',status='unknown',form='formatted')
    write(23,format1) 'zb', 'rb', 'rhov', 'wdcor'
    do ih=1,ngh
        do il=ngl,1,-1
            ihl = ih + (il-1)*ngh
            write(23,format2) gauss%zb(ih),-gauss%rb(il),densities%rv(ihl,1)/gauss%wdcor(ih,il), gauss%wdcor(ih,il)
        enddo
        do il=1,ngl
            ihl = ih + (il-1)*ngh
            write(23,format2) gauss%zb(ih), gauss%rb(il), densities%rv(ihl,1)/ gauss%wdcor(ih,il),gauss%wdcor(ih,il)
        enddo
    enddo      
    close(23)

    open(unit=23,file='dens_prot.out',status='unknown',form='formatted')
    write(23,format1) 'zb', 'rb', 'rhov', 'wdcor'
    do ih=1,ngh
        do il=ngl,1,-1
            ihl = ih + (il-1)*ngh
            write(23,format2)  gauss%zb(ih),- gauss%rb(il),densities%rv(ihl,2)/ gauss%wdcor(ih,il), gauss%wdcor(ih,il)
        enddo
        do il=1,ngl
            ihl = ih + (il-1)*ngh
            write(23,format2)  gauss%zb(ih), gauss%rb(il),densities%rv(ihl,2)/ gauss%wdcor(ih,il), gauss%wdcor(ih,il)
        enddo
    enddo      
    close(23)

    open(unit=23,file='dens_tot.out',status='unknown', form='formatted')
    write(23,format1) 'zb', 'rb', 'rhov', 'wdcor'
    do ih=1,ngh
        do il=ngl,1,-1
            ihl = ih + (il-1)*ngh
            write(23,format2)  gauss%zb(ih),-gauss%rb(il), (densities%rv(ihl,1)+densities%rv(ihl,2))/ gauss%wdcor(ih,il),&
            gauss%wdcor(ih,il)
        enddo
        do il=1,ngl
            ihl = ih + (il-1)*ngh
            write(23,format2)  gauss%zb(ih), gauss%rb(il), (densities%rv(ihl,1)+densities%rv(ihl,2))/ gauss%wdcor(ih,il),&
            gauss%wdcor(ih,il)
        enddo
    enddo      
    close(23)

end subroutine 

subroutine write_densit
    use Globals, only: gauss,densities,outputfile
    integer :: ih,il,ihl
    character(len=*), parameter ::  format1 = "(4(a12,2x))", &
                                    format2 = "(4(f12.6,2x))"
    open(outputfile%u_outputd,file=outputfile%outputd,status='unknown')
    write(outputfile%u_outputd,format1) 'zb', 'rb', 'rhov_n', 'rhov_p'
    do ih=1,ngh
        do il=ngl,1,-1
            ihl = ih + (il-1)*ngh
            write(outputfile%u_outputd,format2) gauss%zb(ih), -gauss%rb(il), densities%rv(ihl,1)/gauss%wdcor(ih,il), &
                densities%rv(ihl,2)/gauss%wdcor(ih,il)
        enddo
        do il=1,ngl
            ihl = ih + (il-1)*ngh
            write(outputfile%u_outputd,format2) gauss%zb(ih), gauss%rb(il), densities%rv(ihl,1)/gauss%wdcor(ih,il), &
                densities%rv(ihl,2)/gauss%wdcor(ih,il)
        enddo
    enddo  
    close(outputfile%u_outputd)
end subroutine



END MODULE Inout