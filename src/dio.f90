!==============================================================================!
! PROGRAM DIO                                                                  !
!                                                                              !
! This program ...                                                             !
!                                                                              !
! Licence:                                                                     !
! Article describing the code:                                                 !
! GitHub repository:                                                           !
!==============================================================================!
PROGRAM DIO
use Globals, only: constraint, iteration, outputfile,expectations,option,pairing
use MathMethods, only: math_gfv
use Inout, only: read_file_b23,read_file_dio,set_output_filename,write_result_DIR
use Field, only: set_woodssaxon_parameters,calculate_meson_propagators,initial_potential_fields,calculate_fields
use Forces, only : set_force_parameters,calculate_density_dependence_of_coupling_constants
use Nucleus, only: set_nucleus_attributes
use Basis, only : set_basis,transform_coefficients_form_cylindrical_to_spherical
use Constraints, only: set_constraint_parameters,calculate_constraint_potential_coefficients
use BCS, only: initial_pairing_field,calculate_occupation,set_block_level_of_KPi
use DiracEquation, only: calculate_sigma_nabla,solve_dirac_equation
use Density, only: calculate_densities_DIR,calculate_densities_RHB
use Expectation, only: calculate_expectation_DIR,calculate_expectation_RHB
use Broyden, only: set_broyden_parameters ,store_initial_matrix_elements_RHB ,mix_potentials_DIR,mix_matrix_elements_RHB
use DeltaField, only: calculate_pairing_matrix_element_W,set_separable_pairing_parameters,initial_delta_field
use RHBEquation, only: solve_RHB_equation
use ExpectationRotation, only: calculate_rotational_correction_energy_DIR
implicit none
logical :: ifPrint=.False.
integer :: constraint_index,iteration_index
character(len=150) :: format1,format2,format3
character(len=20) :: length
format1 = "(110(1h=),/, 48x,'Iteration: ',i4,48x,/,(110(1h-)))"
format2 = "(i3,'. It. si =',f15.10,'  E =',f15.10,'  R =',f15.10,'  b2 =',f15.10,'  b3 =',f15.10)"
! format2 = "(i3,'. It. si =',f30.25,'  E =',f30.25,'  R =',f30.25,'  b2 =',f30.25,'  b3 =',f30.25)"
format3 = "(/,110(1h#),/,'#',27x ,'Iteration converged after',i4,' steps, si =',f13.10,27x,'#',/,110(1h#))"

write(*,*)"PROGRAM Start!"
open(outputfile%u_config, file=outputfile%config, status='unknown')
call math_gfv
call read_file_b23
call read_file_dio(ifPrint .and. .True.)
call set_woodssaxon_parameters(ifPrint .and. .True.)
call set_force_parameters(ifPrint .and. .True.)
call set_nucleus_attributes(ifPrint .and. .True.)
call set_basis(ifPrint .and. .True.)
call set_broyden_parameters

call calculate_sigma_nabla(ifPrint .and. .True.)
call calculate_meson_propagators(ifPrint .and. .True.)
call set_separable_pairing_parameters
call calculate_pairing_matrix_element_W 
close(outputfile%u_config)

open(outputfile%u_rotationalE, file=outputfile%rotationalE, status='unknown')
open(outputfile%u_outExpectation, file=outputfile%outExpectation, status='unknown')

do constraint_index = 1, constraint%length ! loop for different deformation parameters 
    constraint%index = constraint_index  
    call set_constraint_parameters
    call set_output_filename(constraint%betac(constraint%index),constraint%bet3c(constraint%index))
    open(outputfile%u_outputf, file=outputfile%outputf, status='unknown')
    call initial_pairing_field(ifPrint .and. .True.)
    call initial_delta_field
    call initial_potential_fields(ifPrint .and. .True.)

    write(length ,*) constraint%length
    write(*,"('Constraint: ',i4,'/',a, 'beta2=',f6.2,'   beat3=',f6.2)")  &
    constraint%index, adjustl(length),constraint%betac(constraint%index),constraint%bet3c(constraint%index)
    do iteration_index = 1, iteration%iteration_max ! iteration loop
        iteration%ii = iteration_index
        write(outputfile%u_outputf,format1) iteration_index
        if(option%eqType .eq. 0) then
            call solve_dirac_equation(ifPrint .and. .True.)
            call calculate_occupation(ifPrint .and. .True.)
            call calculate_densities_DIR(ifPrint .and. .True.)
        else
            if(iteration%ii .eq. 1) call store_initial_matrix_elements_RHB
            call solve_RHB_equation(ifPrint .and. .True.)
            call calculate_densities_RHB(ifPrint .and. .True.)
        endif
        call calculate_density_dependence_of_coupling_constants(ifPrint .and. .True.)
        call calculate_fields
        if(option%eqType .eq. 0) then
            call calculate_expectation_DIR(ifprint .and. .True.)
        else
            call calculate_expectation_RHB(ifprint .and. .True.)
        endif 
        call calculate_constraint_potential_coefficients(ifPrint .and. .True.)
        if(option%eqType .eq. 0) then
            call mix_potentials_DIR(ifPrint .and. .True.)
        else
            call mix_matrix_elements_RHB
        endif
        write(outputfile%u_outputf,format2)iteration%ii,iteration%si,expectations%ea, &
              expectations%rms,expectations%betg,expectations%beto
        
        if(option%block_method == 2 .and. pairing%allow_block) exit
        if(iteration%ii.ge.2 .and. abs(iteration%si).lt.iteration%epsi) then
            if(option%block_type==0 .or. pairing%allow_block .or. option%block_method==1) exit ! exit iteration
            if(option%block_type==2) call set_block_level_of_KPi(.True.)
            pairing%allow_block = .True.
        endif

    end do
    write(outputfile%u_outputf,format3) iteration%ii,iteration%si
    write(*,"('Iteration converged after',i4,' steps, si =',f13.8, '  ,beta2 =',f6.3,'  ,beta3 =',f6.3)") &
            iteration%ii,iteration%si,expectations%betg,expectations%beto
    if(option%eqType .eq. 0) then
        call transform_coefficients_form_cylindrical_to_spherical(ifPrint .and. .True.)
        call calculate_rotational_correction_energy_DIR
        call write_result_DIR
    endif

    close(outputfile%u_outputf)
end do
close(outputfile%u_rotationalE)
close(outputfile%u_outExpectation)
write(*,*)"PROGRAM END!"
END PROGRAM DIO

!==============================================================================!
! End of file                                                                  !
!==============================================================================!