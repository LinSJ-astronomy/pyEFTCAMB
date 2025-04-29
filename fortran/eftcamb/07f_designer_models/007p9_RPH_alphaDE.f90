!----------------------------------------------------------------------------------------
!
! This file is part of EFTCAMB.
!
! Copyright (C) 2013-2020 by the EFTCAMB authors
!
! The EFTCAMB code is free software;
! You can use it, redistribute it, and/or modify it under the terms
! of the GNU General Public License as published by the Free Software Foundation;
! either version 3 of the License, or (at your option) any later version.
! The full text of the license can be found in the file eftcamb/LICENSE at
! the top level of the EFTCAMB distribution.
!
!----------------------------------------------------------------------------------------

!> @file 007p9_RPH_Omegade.f90
!! This file contains the definition of an EFT reparametrization of Horndeski models.
!! This is described by w_DE and four functions alpha_i = c_i * Omega_DE.

!----------------------------------------------------------------------------------------
!> This module contains the definition of an EFT reparametrization of Horndeski models.
!! This is described by four functions of time and w_DE.
!! Please refer to the numerical notes for details.

!> @author Gen Ye

module EFTCAMB_Reparametrized_Horndeski_alphaDE

   use precision
   use IniObjects
   use EFTCAMB_cache
   use MpiUtils
   use EFT_def
   use EFTCAMB_abstract_parametrizations_1D
   use EFTCAMB_parametrizations_1D
   use EFTCAMB_abstract_model_designer
   use equispaced_linear_interpolation_1D
   use MassiveNu

   implicit none

   private

   public EFTCAMB_RPH_alphaDE

   !----------------------------------------------------------------------------------------
   !> This is the type that contains the definition of RPH.
   type, extends ( EFTCAMB_designer_model ) :: EFTCAMB_RPH_alphaDE

      ! the background DE model selection flags:
      integer  :: RPHalphaDEwDEmodel              !< Model selection flag for RPH w DE.

      ! the c_i coefficients of alpha_i = c_i * Omega_DE
      real(dl)  :: RPHalphaDEc_M      !< coefficient for alpha_M.
      real(dl)  :: RPHalphaDEc_K      !< coefficient for alpha_K.
      real(dl)  :: RPHalphaDEc_B      !< coefficient for alpha_B.
      real(dl)  :: RPHalphaDEc_T      !< coefficient for alpha_T.

      ! the RPH functions:
      class( parametrized_function_1D ), allocatable    :: RPHalphaDE_wDE         !< The RPH function w_DE.

      type( equispaced_linear_interpolate_function_1D ) :: RPHalphaDE_c
      type( equispaced_linear_interpolate_function_1D ) :: RPHalphaDE_dc
      type( equispaced_linear_interpolate_function_1D ) :: RPHalphaDE_Lmd
      type( equispaced_linear_interpolate_function_1D ) :: RPHalphaDE_dLmd
      type( equispaced_linear_interpolate_function_1D ) :: RPHalphaDE_d2Lmd
      type( equispaced_linear_interpolate_function_1D ) :: RPHalphaDE_Omg
      type( equispaced_linear_interpolate_function_1D ) :: RPHalphaDE_dOmg
      type( equispaced_linear_interpolate_function_1D ) :: RPHalphaDE_d2Omg
      type( equispaced_linear_interpolate_function_1D ) :: RPHalphaDE_d3Omg
      type( equispaced_linear_interpolate_function_1D ) :: RPHalphaDE_d4Omg
      ! type( equispaced_linear_interpolate_function_1D ) :: RPHalphaDE_gm1
      ! type( equispaced_linear_interpolate_function_1D ) :: RPHalphaDE_dgm1
      type( equispaced_linear_interpolate_function_1D ) :: RPHalphaDE_aK
      type( equispaced_linear_interpolate_function_1D ) :: RPHalphaDE_daK
      type( equispaced_linear_interpolate_function_1D ) :: RPHalphaDE_gm2
      type( equispaced_linear_interpolate_function_1D ) :: RPHalphaDE_dgm2
      type( equispaced_linear_interpolate_function_1D ) :: RPHalphaDE_d2gm2
      type( equispaced_linear_interpolate_function_1D ) :: RPHalphaDE_d3gm2
      type( equispaced_linear_interpolate_function_1D ) :: RPHalphaDE_gm3
      type( equispaced_linear_interpolate_function_1D ) :: RPHalphaDE_dgm3
      type( equispaced_linear_interpolate_function_1D ) :: RPHalphaDE_d2gm3
      type( equispaced_linear_interpolate_function_1D ) :: RPHalphaDE_d3gm3
      type( equispaced_linear_interpolate_function_1D ) :: RPHalphaDE_d4gm3


      type( equispaced_linear_interpolate_function_1D ) :: RPHalphaDE_Ode        !< The interpolated Omega_DE.
      type( equispaced_linear_interpolate_function_1D ) :: RPHalphaDE_Odedot     !< The interpolated d Omega_DE/ dt, conformal time.
      ! type( equispaced_linear_interpolate_function_1D ) :: RPHalphaDE_d2Ode        !< The interpolated d^2 Omega_DE/ da^2.
      ! type( equispaced_linear_interpolate_function_1D ) :: RPHalphaDE_d3Ode        !< The interpolated d^3 Omega_DE/ da^3.
      ! type( equispaced_linear_interpolate_function_1D ) :: RPHalphaDE_d4Ode        !< The interpolated d^4 Omega_DE/ da^4.
      type( equispaced_linear_interpolate_function_1D ) :: RPHalphaDE_M        !< The interpolated conformal Hubble parameter.
      type( equispaced_linear_interpolate_function_1D ) :: RPHalphaDE_dM        !< The interpolated d M/ da.
      ! type( equispaced_linear_interpolate_function_1D ) :: RPHalphaDE_d2M        !< The interpolated d^2 M/ da^2.
      ! type( equispaced_linear_interpolate_function_1D ) :: RPHalphaDE_d3M        !< The interpolated d^3 M/ da^3.
      ! type( equispaced_linear_interpolate_function_1D ) :: RPHalphaDE_d4M        !< The interpolated d^4 M/ da^4.

      ! background solver parameters:
      integer  :: designer_num_points = 1000                            !< Number of points sampled by the designer code.
      real(dl) :: x_initial           = log(10._dl**(-8))               !< log(a start)
      real(dl) :: x_final             = 0.0_dl                          !< log(a final)

   contains

      ! initialization of the model:
      procedure :: read_model_selection            => EFTCAMBRPHalphaDEReadModelSelectionFromFile  !< subroutine that reads the parameters of the model from file
      procedure :: allocate_model_selection        => EFTCAMBRPHalphaDEAllocateModelSelection      !< subroutine that allocates the model selection.
      procedure :: init_model_parameters           => EFTCAMBRPHalphaDEInitModelParameters         !< subroutine that initializes the model parameters based on the values found in an input array.
      procedure :: init_model_parameters_from_file => EFTCAMBRPHalphaDEInitModelParametersFromFile !< subroutine that reads the parameters of the model from file.

      ! background solver:
      procedure :: initialize_background           => EFTCAMBRPHalphaDEInitBackground               !< subroutine that initializes the background of designer DGP.
      procedure :: solve_designer_equations        => EFTCAMBRPHalphaDESolveDesignerEquations       !< subroutine that solves the designer DGP background equations.

      ! utility functions:
      procedure :: compute_param_number  => EFTCAMBRPHalphaDEComputeParametersNumber    !< subroutine that computes the number of parameters of the model.
      procedure :: feedback              => EFTCAMBRPHalphaDEFeedback                   !< subroutine that prints on the screen feedback information about the model.
      procedure :: parameter_names       => EFTCAMBRPHalphaDEParameterNames             !< subroutine that returns the i-th parameter name of the model.
      procedure :: parameter_names_latex => EFTCAMBRPHalphaDEParameterNamesLatex        !< subroutine that returns the i-th parameter name of the model.
      procedure :: parameter_values      => EFTCAMBRPHalphaDEParameterValues            !< subroutine that returns the i-th parameter value.

      ! CAMB related procedures:
      procedure :: compute_background_EFT_functions  => EFTCAMBRPHalphaDEBackgroundEFTFunctions   !< subroutine that computes the value of the background EFT functions at a given time.
      procedure :: compute_secondorder_EFT_functions => EFTCAMBRPHalphaDESecondOrderEFTFunctions  !< subroutine that computes the value of the second order EFT functions at a given time.
      procedure :: compute_dtauda                    => EFTCAMBRPHalphaDEComputeDtauda            !< function that computes dtauda = 1/sqrt(a^2H^2).
      procedure :: compute_adotoa                    => EFTCAMBRPHalphaDEComputeAdotoa            !< subroutine that computes adotoa = H and its two derivatives wrt conformal time.
      procedure :: compute_H_derivs                  => EFTCAMBRPHalphaDEComputeHubbleDer         !< subroutine that computes the two derivatives wrt conformal time of H.
      ! stability procedures:
      procedure :: additional_model_stability        => EFTCAMBRPHalphaDEAdditionalModelStability !< function that computes model specific stability requirements.

   end type EFTCAMB_RPH_alphaDE

   ! ---------------------------------------------------------------------------------------------

contains

   ! ---------------------------------------------------------------------------------------------
   !> Subroutine that reads the parameters of the model from file.
   subroutine EFTCAMBRPHalphaDEReadModelSelectionFromFile( self, Ini, eft_error )

      implicit none

      class(EFTCAMB_RPH_alphaDE) :: self      !< the base class
      type(TIniFile)     :: Ini       !< Input ini file
      integer            :: eft_error !< error code: 0 all fine, 1 initialization failed

      ! read model selection flags:
      self%RPHalphaDEwDEmodel             = Ini%Read_Int( 'RPHalphaDEwDEmodel'            , 0 )
      self%x_initial                      = Log(Ini%Read_Double( 'RPHalphaDEa_ini'            , 10._dl**(-8) ))

   end subroutine EFTCAMBRPHalphaDEReadModelSelectionFromFile

   ! ---------------------------------------------------------------------------------------------
   !> Subroutine that allocates the model selection.
   subroutine EFTCAMBRPHalphaDEAllocateModelSelection( self, Ini, eft_error )

      implicit none

      class(EFTCAMB_RPH_alphaDE)  :: self      !< the base class
      type(TIniFile)      :: Ini       !< Input ini file
      integer             :: eft_error !< error code: 0 all fine, 1 initialization failed

      integer                     :: temp_feedback 

      ! get feedback flag:
      temp_feedback = Ini%Read_Int('feedback_level', 0)
      
      ! allocate wDE:
      if ( allocated(self%RPHalphaDE_wDE) ) deallocate(self%RPHalphaDE_wDE)
      select case ( self%RPHalphaDEwDEmodel )
       case(0)
         allocate( wDE_LCDM_parametrization_1D::self%RPHalphaDE_wDE )
         call self%RPHalphaDE_wDE%set_name( 'RPHalphaDEw', 'w' )
       case(1)
         allocate( constant_parametrization_1D::self%RPHalphaDE_wDE )
         call self%RPHalphaDE_wDE%set_name( 'RPHalphaDEw', 'w' )
       case(2)
         allocate( CPL_parametrization_1D::self%RPHalphaDE_wDE )
         call self%RPHalphaDE_wDE%set_param_names( ['RPHalphaDEw0', 'RPHalphaDEwa'], ['w_0', 'w_a'] )
         call self%RPHalphaDE_wDE%set_name( 'RPHalphaDEw', 'w' )
       case(3)
         allocate( JBP_parametrization_1D::self%RPHalphaDE_wDE )
         call self%RPHalphaDE_wDE%set_param_names( ['RPHalphaDEw0', 'RPHalphaDEwa', 'RPHalphaDEwn'], [ 'w_0', 'w_a', 'n  ' ] )
         call self%RPHalphaDE_wDE%set_name( 'RPHalphaDEw', 'w' )
       case(4)
         allocate( turning_point_parametrization_1D::self%RPHalphaDE_wDE )
         call self%RPHalphaDE_wDE%set_param_names( ['RPHalphaDEw0 ', 'RPHalphaDEwa ', 'RPHalphaDEwat'], ['w_0', 'w_a', 'a_t'] )
         call self%RPHalphaDE_wDE%set_name( 'RPHalphaDEw', 'w' )
       case(5)
         allocate( taylor_parametrization_1D::self%RPHalphaDE_wDE )
         call self%RPHalphaDE_wDE%set_param_names( ['RPHalphaDEw0', 'RPHalphaDEwa', 'RPHalphaDEw2', 'RPHalphaDEw3'], ['w_0', 'w_a', 'w_2', 'w_3'] )
         call self%RPHalphaDE_wDE%set_name( 'RPHalphaDEw', 'w' )
       case default
         call allocate_parametrized_1D_function( self%RPHalphaDE_wDE, self%RPHalphaDEwDEmodel, 'RPHalphaDEw', 'w', eft_error, temp_feedback )
         if ( eft_error == 1 ) return
      end select

      ! additional initialization of the function:
      call self%RPHalphaDE_wDE%init_func_from_file       ( Ini, eft_error )

   end subroutine EFTCAMBRPHalphaDEAllocateModelSelection

   ! ---------------------------------------------------------------------------------------------
   !> Subroutine that initializes the model parameters based on the values found in an input array.
   subroutine EFTCAMBRPHalphaDEInitModelParameters( self, array )

      implicit none

      class(EFTCAMB_RPH_alphaDE)                            :: self   !< the base class
      real(dl), dimension(self%parameter_number), intent(in) :: array  !< input array with the values of the parameters of the model.

      real(dl), allocatable, dimension(:)                    :: temp
      integer :: num_params_wDE, i

      ! first four are the c_i's, in order M K B T
      self%RPHalphaDEc_M = array(1)
      self%RPHalphaDEc_K = array(2)
      self%RPHalphaDEc_B = array(3)
      self%RPHalphaDEc_T = array(4)
      ! then w_DE parameters:
      num_params_wDE = self%RPHalphaDE_wDE%parameter_number
      allocate( temp(num_params_wDE) )
      do i = 1, num_params_wDE
         temp(i)         = array(i + 4)
      end do
      call self%RPHalphaDE_wDE%init_parameters(temp)
      deallocate( temp )

      ! now check the length of the parameters:
      if ( num_params_wDE + 4 /= self%parameter_number ) then
         write(*,*) 'In EFTCAMBRPHalphaDEInitModelParameters:'
         write(*,*) 'Length of self%parameter_number does not coincide with the number of required model parameters.'
         write(*,*) 'number of required model parameters:', num_params_wDE + 4
         write(*,*) 'self%parameter_number:', self%parameter_number
         call MpiStop('EFTCAMB error')
      end if

   end subroutine EFTCAMBRPHalphaDEInitModelParameters

   ! ---------------------------------------------------------------------------------------------
   !> Subroutine that reads the parameters of the model from file.
   subroutine EFTCAMBRPHalphaDEInitModelParametersFromFile( self, Ini, eft_error )

      implicit none

      class(EFTCAMB_RPH_alphaDE) :: self      !< the base class
      type(TIniFile)     :: Ini       !< Input ini file
      integer            :: eft_error !< error code: 0 all fine, 1 initialization failed

      ! read c_i's
      self%RPHalphaDEc_M = Ini%Read_Double( 'EFTRPHalphaDEc_M', 0._dl )
      self%RPHalphaDEc_K = Ini%Read_Double( 'EFTRPHalphaDEc_K', 0._dl )
      self%RPHalphaDEc_B = Ini%Read_Double( 'EFTRPHalphaDEc_B', 0._dl )
      self%RPHalphaDEc_T = Ini%Read_Double( 'EFTRPHalphaDEc_T', 0._dl )

      ! read wDE pars
      call self%RPHalphaDE_wDE%init_from_file( Ini, eft_error )

   end subroutine EFTCAMBRPHalphaDEInitModelParametersFromFile

   ! ---------------------------------------------------------------------------------------------
   !> Subroutine that computes the number of parameters of the model.
   subroutine EFTCAMBRPHalphaDEComputeParametersNumber( self )

      implicit none

      class(EFTCAMB_RPH_alphaDE)  :: self   !< the base class

      self%parameter_number = self%RPHalphaDE_wDE%parameter_number + 4

   end subroutine EFTCAMBRPHalphaDEComputeParametersNumber

   ! ---------------------------------------------------------------------------------------------
   !> Subroutine that prints on the screen feedback information about the model.
   subroutine EFTCAMBRPHalphaDEFeedback( self, print_params )

      implicit none

      class(EFTCAMB_RPH_alphaDE)  :: self         !< the base class
      logical, optional   :: print_params !< optional flag that decised whether to print numerical values
      !! of the parameters.

      ! print general model informations:
      write(*,*)
      write(*,'(a,a)')    '   Model               =  ', self%name
      write(*,'(a,I3)')   '   Number of params    ='  , self%parameter_number

      ! print model functions informations:
      write(*,*)
      if ( self%RPHalphaDEwDEmodel             /= 0 ) write(*,'(a,I3)') '   RPHalphaDEwDEmodel             =', self%RPHalphaDEwDEmodel

      write(*,*)
      ! print functions informations:
      call self%RPHalphaDE_wDE%feedback        ( print_params )

   end subroutine EFTCAMBRPHalphaDEFeedback

   ! ---------------------------------------------------------------------------------------------
   !> Subroutine that returns the i-th parameter name of the model
   subroutine EFTCAMBRPHalphaDEParameterNames( self, i, name )

      implicit none

      class(EFTCAMB_RPH_alphaDE) :: self   !< the base class
      integer     , intent(in)    :: i      !< the index of the parameter
      character(*), intent(out)   :: name   !< the output name of the i-th parameter

      integer  :: j

      ! check validity of input:
      if ( i > self%parameter_number .or. i <= 0 ) then
         write(*,'(a,I3)') 'No parameter corresponding to: ', i
         write(*,'(a,I3)') 'Total number of parameters is: ', self%parameter_number
         call MpiStop('EFTCAMB error')

         ! parameter from wDE function
      else if ( i <= 4 ) then
         select case ( i )
          case(1)
            name = TRIM('EFTRPHalphaDEc_M')
          case(2)
            name = TRIM('EFTRPHalphaDEc_K')
          case(3)
            name = TRIM('EFTRPHalphaDEc_B')
          case(4)
            name = TRIM('EFTRPHalphaDEc_T')
         end select
         return

         ! parameter from wDE function
      else
         do j = 1, self%RPHalphaDE_wDE%parameter_number
            if ( i - 4 == j ) call self%RPHalphaDE_wDE%parameter_names( j, name )
         end do
         return

      end if

   end subroutine EFTCAMBRPHalphaDEParameterNames

   ! ---------------------------------------------------------------------------------------------
   !> Subroutine that returns the i-th parameter name of the model
   subroutine EFTCAMBRPHalphaDEParameterNamesLatex( self, i, latexname )

      implicit none

      class(EFTCAMB_RPH_alphaDE) :: self       !< the base class
      integer     , intent(in)    :: i          !< The index of the parameter
      character(*), intent(out)   :: latexname  !< the output latex name of the i-th parameter

      integer  :: j

      ! check validity of input:
      if ( i > self%parameter_number .or. i <= 0 ) then
         write(*,'(a,I3)') 'No parameter corresponding to: ', i
         write(*,'(a,I3)') 'Total number of parameters is: ', self%parameter_number
         call MpiStop('EFTCAMB error')

         ! parameter from wDE function
      else if ( i <= 4 ) then
         select case ( i )
          case(1)
            latexname = TRIM('c_M')
          case(2)
            latexname = TRIM('c_K')
          case(3)
            latexname = TRIM('c_B')
          case(4)
            latexname = TRIM('c_T')
         end select
         return

         ! parameter from wDE function
      else
         do j = 1, self%RPHalphaDE_wDE%parameter_number
            if ( i - 4 == j ) call self%RPHalphaDE_wDE%parameter_names_latex( j, latexname )
         end do
         return

      end if

   end subroutine EFTCAMBRPHalphaDEParameterNamesLatex

   ! ---------------------------------------------------------------------------------------------
   !> Subroutine that returns the i-th parameter name of the model
   subroutine EFTCAMBRPHalphaDEParameterValues( self, i, value )

      implicit none

      class(EFTCAMB_RPH_alphaDE) :: self   !< the base class
      integer , intent(in)        :: i      !< The index of the parameter
      real(dl), intent(out)       :: value  !< the output value of the i-th parameter

      integer  :: j

      ! check validity of input:
      if ( i > self%parameter_number .or. i <= 0 ) then
         write(*,'(a,I3)') 'No parameter corresponding to: ', i
         write(*,'(a,I3)') 'Total number of parameters is: ', self%parameter_number
         call MpiStop('EFTCAMB error')

         ! parameter from wDE function
      else if ( i <= 4 ) then
         select case ( i )
          case(1)
            value = self%RPHalphaDEc_M
          case(2)
            value = self%RPHalphaDEc_K
          case(3)
            value = self%RPHalphaDEc_B
          case(4)
            value = self%RPHalphaDEc_T
         end select
         return

         ! parameter from wDE function
      else
         do j = 1, self%RPHalphaDE_wDE%parameter_number
            if ( i - 4 == j ) call self%RPHalphaDE_wDE%parameter_value( j, value )
         end do
         return

      end if

   end subroutine EFTCAMBRPHalphaDEParameterValues

   ! ---------------------------------------------------------------------------------------------
   !> Subroutine that initializes the background.
   subroutine EFTCAMBRPHalphaDEInitBackground( self, params_cache, feedback_level, success, outroot )

      implicit none

      class(EFTCAMB_RPH_alphaDE)                    :: self           !< the base class
      type(TEFTCAMB_parameter_cache), intent(inout) :: params_cache   !< a EFTCAMB parameter cache containing cosmological parameters
      integer                       , intent(in)    :: feedback_level !< level of feedback from the background code. 0=none; 1=some; 2=chatty.
      logical                       , intent(out)   :: success        !< wether the background initialization succeded or not
      character(LEN=*), optional    , intent(in)    :: outroot        !< the output root for the debug files

      ! some feedback:
      if ( feedback_level>0 ) then
         write(*,'(a)') "***************************************************************"
         write(*,'(a)') ' EFTCAMB designer RPH background solver'
         write(*,'(a)')
      end if

      ! initialize interpolating functions:
      call self%RPHalphaDE_c%initialize  ( self%designer_num_points, self%x_initial, self%x_final )
      call self%RPHalphaDE_dc%initialize  ( self%designer_num_points, self%x_initial, self%x_final )
      call self%RPHalphaDE_Lmd%initialize  ( self%designer_num_points, self%x_initial, self%x_final )
      call self%RPHalphaDE_dLmd%initialize  ( self%designer_num_points, self%x_initial, self%x_final )
      call self%RPHalphaDE_d2Lmd%initialize  ( self%designer_num_points, self%x_initial, self%x_final )
      call self%RPHalphaDE_Omg%initialize  ( self%designer_num_points, self%x_initial, self%x_final )
      call self%RPHalphaDE_dOmg%initialize  ( self%designer_num_points, self%x_initial, self%x_final )
      call self%RPHalphaDE_d2Omg%initialize  ( self%designer_num_points, self%x_initial, self%x_final )
      call self%RPHalphaDE_d3Omg%initialize  ( self%designer_num_points, self%x_initial, self%x_final )
      call self%RPHalphaDE_d4Omg%initialize  ( self%designer_num_points, self%x_initial, self%x_final )
      ! call self%RPHalphaDE_gm1%initialize  ( self%designer_num_points, self%x_initial, self%x_final )
      ! call self%RPHalphaDE_dgm1%initialize  ( self%designer_num_points, self%x_initial, self%x_final )
      call self%RPHalphaDE_aK%initialize  ( self%designer_num_points, self%x_initial, self%x_final )
      call self%RPHalphaDE_daK%initialize  ( self%designer_num_points, self%x_initial, self%x_final )
      call self%RPHalphaDE_gm2%initialize  ( self%designer_num_points, self%x_initial, self%x_final )
      call self%RPHalphaDE_dgm2%initialize  ( self%designer_num_points, self%x_initial, self%x_final )
      call self%RPHalphaDE_d2gm2%initialize  ( self%designer_num_points, self%x_initial, self%x_final )
      call self%RPHalphaDE_d3gm2%initialize  ( self%designer_num_points, self%x_initial, self%x_final )
      call self%RPHalphaDE_gm3%initialize  ( self%designer_num_points, self%x_initial, self%x_final )
      call self%RPHalphaDE_dgm3%initialize  ( self%designer_num_points, self%x_initial, self%x_final )
      call self%RPHalphaDE_d2gm3%initialize  ( self%designer_num_points, self%x_initial, self%x_final )
      call self%RPHalphaDE_d3gm3%initialize  ( self%designer_num_points, self%x_initial, self%x_final )
      call self%RPHalphaDE_d4gm3%initialize  ( self%designer_num_points, self%x_initial, self%x_final )


      call self%RPHalphaDE_M%initialize  ( self%designer_num_points, self%x_initial, self%x_final )
      call self%RPHalphaDE_Ode%initialize  ( self%designer_num_points, self%x_initial, self%x_final )
      call self%RPHalphaDE_Odedot%initialize  ( self%designer_num_points, self%x_initial, self%x_final )
      ! call self%RPHalphaDE_d2Ode%initialize  ( self%designer_num_points, self%x_initial, self%x_final )
      ! call self%RPHalphaDE_d3Ode%initialize  ( self%designer_num_points, self%x_initial, self%x_final )
      ! call self%RPHalphaDE_d4Ode%initialize  ( self%designer_num_points, self%x_initial, self%x_final )
      call self%RPHalphaDE_dM%initialize  ( self%designer_num_points, self%x_initial, self%x_final )
      ! call self%RPHalphaDE_d2M%initialize  ( self%designer_num_points, self%x_initial, self%x_final )
      ! call self%RPHalphaDE_d3M%initialize  ( self%designer_num_points, self%x_initial, self%x_final )
      ! call self%RPHalphaDE_d4M%initialize  ( self%designer_num_points, self%x_initial, self%x_final )

      success = .True.
      ! solve the background equations and store the solution:
      call self%solve_designer_equations( params_cache, success=success, feedback_level=feedback_level )

   end subroutine EFTCAMBRPHalphaDEInitBackground

   ! ---------------------------------------------------------------------------------------------
   !> Subroutine that solves for M by integrating apha_M
   subroutine EFTCAMBRPHalphaDESolveDesignerEquations( self, params_cache, success, feedback_level )

      implicit none

      class(EFTCAMB_RPH_alphaDE)                   :: self          !< the base class.
      type(TEFTCAMB_parameter_cache), intent(in)   :: params_cache  !< a EFTCAMB parameter cache containing cosmological parameters.
      logical , intent(out)                        :: success       !< whether the calculation ended correctly or not
      integer , intent(in)                         :: feedback_level   !< whether be noisy

      integer, parameter :: num_eq = 1   !<  Number of equations
      real(dl) :: y(num_eq), ydot(num_eq)

      ! real(dl) :: Omegam_EFT, Omegavac_EFT, Omegar_EFT
      ! real(dl) :: EFT_E2_DE, OmegaDE_t
      ! real(dl) :: EFT_E2_nu, EFT_E2P_nu, rhonu_tot, presnu_tot, rhonu, presnu, grhormass_t
      ! integer :: nu_i
      real(dl) :: grhom0, grhor0, grhov0

      integer  :: itol, itask, istate, iopt, LRN, LRS, LRW, LIS, LIN, LIW, JacobianMode, i
      real(dl) :: rtol, atol, t1, t2
      real(dl), allocatable :: rwork(:)
      integer,  allocatable :: iwork(:)

      ! This routine integrates the simple equation \log M = \int d\log a \alpha_M

      ! 1) Cosmological parameters:
      grhom0 = params_cache%grhob + params_cache%grhoc
      grhor0 = 3._dl * params_cache%h0_Mpc**2 * (params_cache%omegag + params_cache%omegar)
      grhov0 = params_cache%grhov
      ! Omegam_EFT         = params_cache%omegab + params_cache%omegac
      ! Omegavac_EFT       = params_cache%omegav
      ! Omegar_EFT         = params_cache%omegag + params_cache%omegar

      ! 2) Initial values
      y(1) = 0._dl

      ! 3) Initialize DLSODA:
      ! set-up the relative and absolute tollerances:
      itol = 1
      rtol = 1.d-12
      atol = 1.d-16
      ! initialize task to do:
      itask  = 1
      istate = 1
      iopt   = 1
      ! initialize the work space:
      LRN = 20 + 16*num_eq
      LRS = 22 + 9*num_eq + num_eq**2
      LRW = max(LRN,LRS)
      LIS = 20 + num_eq
      LIN = 20
      LIW = max(LIS,LIN)
      ! allocate the arrays:
      allocate(rwork(LRW))
      allocate(iwork(LIW))
      ! optional lsoda input:
      RWORK(5) = 0._dl  ! the step size to be attempted on the first step. The default value is determined by the solver.
      RWORK(6) = 0._dl  ! the maximum absolute step size allowed. The default value is infinite.
      RWORK(7) = 0._dl  ! the minimum absolute step size allowed. The default value is 0.
      IWORK(5) = 0      ! flag to generate extra printing at method switches. IXPR = 0 means no extra printing (the default). IXPR = 1 means print data on each switch.
      IWORK(6) = 100    ! maximum number of (internally defined) steps allowed during one call to the solver. The default value is 500.
      IWORK(7) = 0      ! maximum number of messages printed (per problem) warning that T + H = T on a step (H = step size). This must be positive to result in a non-default value.  The default value is 10.
      IWORK(8) = 0      ! the maximum order to be allowed for the nonstiff (Adams) method.  the default value is 12.
      IWORK(9) = 0      ! the maximum order to be allowed for the stiff (BDF) method.  The default value is 5.
      ! additional lsoda stuff:
      CALL XSETF(0) ! suppress odepack printing
      ! Jacobian mode: 1=fullJacobian, 2=not provided
      JacobianMode = 2

      ! initial time:
      t1  = self%RPHalphaDE_gm3%x(1)
      ! store initial step:
      call output( num_eq, 1, t1, y)

      ! 3) solve for M:
      do i=1, self%RPHalphaDE_gm3%num_points-1

         ! set the time step:
         t1 = self%RPHalphaDE_gm3%x(i)
         t2 = self%RPHalphaDE_gm3%x(i+1)
         ! solve the system:
         call DLSODA ( derivs, num_eq, y, t1, t2, itol, rtol, atol, itask, istate, iopt, RWORK, LRW, IWORK, LIW, jacobian, JacobianMode)
         ! check istate for LSODA good completion:
         if ( istate < 0 ) then
            if ( istate == -1 ) then
               if ( feedback_level>1 ) write(*,*) 'DLSODA excessive work'
               istate = 1
            else
               success = .False.
               if ( feedback_level>1 ) write(*,*) ' EFTCAMB RPHalphaDE failed to integrate M_eff in segment a = ', Exp(t2), Exp(t1), istate
               return
            end if
         end if

         ! output one step:
         call output( num_eq, i+1, t2, y)

      end do

      return

   contains
      ! ---------------------------------------------------------------------------------------------
      !> Subroutine that computes y' given y
      subroutine derivs( num_eq, x, y, ydot )

         implicit none

         integer , intent(in)                     :: num_eq !< number of equations in the ODE system
         real(dl), intent(in)                     :: x      !< time at which the derivatives of the system are computed
         real(dl), intent(in) , dimension(num_eq) :: y      !< input status of the system
         real(dl), intent(out), dimension(num_eq) :: ydot   !< value of the derivative at time x

         real(dl) :: Ode
         real(dl) :: a
         real(dl) :: grhor, grhom, grhov, grhon
         real(dl) :: grhonu, gpinu, grhormass_t
         integer  :: nu_i
         
         a = Exp(x)

         grhom = grhom0 / a
         grhor = grhor0 / a**2
         grhov = grhov0 * self%RPHalphaDE_wDE%integral(a)
         grhon = 0._dl
         if ( params_cache%Num_Nu_Massive /= 0 ) then
            do nu_i = 1, params_cache%Nu_mass_eigenstates
               grhonu      = 0._dl
               gpinu       = 0._dl
               grhormass_t = params_cache%grhormass(nu_i)/a**2
               call ThermalNuBack%rho_P(a*params_cache%nu_masses(nu_i), grhonu, gpinu )
               grhon = grhon +grhormass_t*grhonu
           end do
         end if
         Ode = grhov / ( grhom + grhor + grhov + grhon )

         ydot(1) = self%RPHalphaDEc_M * Ode

      end subroutine derivs

      ! ---------------------------------------------------------------------------------------------
      !> Subroutine that computes the Jacobian of the system. Now a dummy function.
      !! Implementing it might increase performances.
      subroutine jacobian( num_eq, x, y, ml, mu, pd, nrowpd )

         implicit none

         integer                            :: num_eq !< number of components of the Jacobian
         integer                            :: ml     !< ignored
         integer                            :: mu     !< ignored
         integer                            :: nrowpd !< ignored
         real(dl)                           :: x      !< time at which the Jacobian is computed
         real(dl), dimension(num_eq)        :: y      !< input status of the system
         real(dl), dimension(nrowpd,num_eq) :: pd     !< output Jacobian

      end subroutine jacobian

      ! ---------------------------------------------------------------------------------------------
      !> Subroutine that takes the solution of the background DGP equations and stores the values of the EFT functions.
      subroutine output( num_eq, ind, x, y)

         implicit none

         integer , intent(in)                     :: num_eq !< number of equations in the ODE system.
         integer , intent(in)                     :: ind    !< index of the EFT functions interpolation tables to fill.
         real(dl), intent(in)                     :: x      !< time at which the derivatives of the system are computed.
         real(dl), intent(in) , dimension(num_eq) :: y      !< input status of the system.

         real(dl) :: a
         real(dl) :: wDE, dwDE, d2wDE, d3wDE, d4wDE
         real(dl) :: H, dH, d2H, d3H, d4H, ca, Ht, Ht2, Ht3, Ht4
         real(dl) :: grhor, grhom, grhov, grhon, gpn, gdpn, gddpn, gdddpn
         real(dl) :: grhonu, gpinu, gpinudot, gpinudotdot, gpinudotdotdot, grhormass_t
         real(dl) :: Ode, dOde, d2Ode, d3Ode, d4Ode
         integer  :: nu_i
         real(dl) :: M, Mp1, dM, d2M, d3M, d4M, aT, daT, d2aT, d3aT, d4aT, aK, daK, aB, daB, d2aB, d3aB, aM, daM, d2aM, d3aM
         real(dl) :: Omg, dOmg, d2Omg, d3Omg, d4Omg
         real(dl) :: c, dc

         a = Exp(x)
         ! short vars for wDE
         wDE   = self%RPHalphaDE_wDE%value(a)
         dwDE  = a * self%RPHalphaDE_wDE%first_derivative(a)
         d2wDE = a**2 * self%RPHalphaDE_wDE%second_derivative(a)
         d3wDE = a**3 * self%RPHalphaDE_wDE%third_derivative(a)

         !1) Hubble and derivatives
         ! H, Ht
         grhom = grhom0 / a
         grhor = grhor0 / a**2
         grhov = grhov0 * self%RPHalphaDE_wDE%integral(a)
         grhon = 0._dl
         gpn = 0._dl
         gdpn = 0._dl
         gddpn = 0._dl
         gdddpn = 0._dl
         if ( params_cache%Num_Nu_Massive /= 0 ) then
            do nu_i = 1, params_cache%Nu_mass_eigenstates
               grhonu      = 0._dl
               gpinu       = 0._dl
               grhormass_t = params_cache%grhormass(nu_i)/a**2
               call ThermalNuBack%rho_P(a*params_cache%nu_masses(nu_i), grhonu, gpinu )
               grhon = grhon +grhormass_t*grhonu
               gpn   = gpn  +grhormass_t*gpinu
           end do
         end if
         H = sqrt( ( grhom + grhor + grhov + grhon )/3._dl )
         Ht = -0.5_dl*( H**2 + grhor/3._dl + grhov*wDE + gpn )
         ! Ht2
         grhon = 0._dl
         gpn = 0._dl
         gdpn = 0._dl
         gddpn = 0._dl
         gdddpn = 0._dl
         if ( params_cache%Num_Nu_Massive /= 0 ) then
            do nu_i = 1, params_cache%Nu_mass_eigenstates
               grhonu      = 0._dl
               gpinu       = 0._dl
               gpinudot    = 0._dl
               grhormass_t = params_cache%grhormass(nu_i)/a**2
               call ThermalNuBack%rho_P(a*params_cache%nu_masses(nu_i), grhonu, gpinu )
               grhon    = grhon +grhormass_t*grhonu
               gpn      = gpn  +grhormass_t*gpinu
               gpinudot = ThermalNuBack%pidot( a*params_cache%nu_masses(nu_i), H, gpinu )
               gdpn     = gdpn  + grhormass_t*( gpinudot - 4._dl*H*gpinu)
           end do
         end if
         Ht2 = H*( grhom/6._dl + 2._dl*grhor/3._dl ) &
         & + H*grhov*( 1._dl/6._dl +wDE +1.5_dl*wDE**2 -0.5_dl*dwDE ) &
         & + H*grhon/6._dl -0.5_dl*H*gpn -0.5_dl*gdpn
         ! Ht3
         grhon = 0._dl
         gpn = 0._dl
         gdpn = 0._dl
         gddpn = 0._dl
         gdddpn = 0._dl
         if ( params_cache%Num_Nu_Massive /= 0 ) then
            do nu_i = 1, params_cache%Nu_mass_eigenstates
               grhonu      = 0._dl
               gpinu       = 0._dl
               gpinudot    = 0._dl
               gpinudotdot = 0._dl
               grhormass_t = params_cache%grhormass(nu_i)/a**2
               call ThermalNuBack%rho_P(a*params_cache%nu_masses(nu_i), grhonu, gpinu )
               grhon       = grhon +grhormass_t*grhonu
               gpn         = gpn  +grhormass_t*gpinu
               gpinudot    = ThermalNuBack%pidot( a*params_cache%nu_masses(nu_i), H, gpinu )
               gdpn        = gdpn  + grhormass_t*( gpinudot - 4._dl*H*gpinu)
               gpinudotdot = ThermalNuBack%pidotdot( a*params_cache%nu_masses(nu_i), H, Ht, gpinu, gpinudot)
               gddpn       = gddpn - 4._dl*H*grhormass_t*( gpinudot - 4._dl*H*gpinu )+ grhormass_t*( gpinudotdot - 4._dl*Ht*gpinu - 4._dl*H*gpinudot )
           end do
         end if
         Ht3 = H**2*( -grhom/6._dl - 4._dl*grhor/3._dl &
            & + grhov*(-1._dl/6._dl - 1.5_dl*wDE - 4.5_dl*wDE**2 - 4.5_dl*wDE**3 + dwDE + 4.5_dl*wDE*dwDE - 0.5_dl*d2wDE) &
            & -1._dl/6._dl*grhon - 1.5_dl*gpn ) &
            & + Ht*( grhom/6._dl + 2._dl/3._dl*grhor &
            & + grhov*(1._dl/6._dl + wDE +1.5_dl*wDE**2 - 0.5_dl*dwDE) &
            & + grhon/6._dl - 0.5_dl*gpn ) &
            & -1.5_dl*H*gdpn - 0.5_dl*gddpn
         ! Ht4
         grhon = 0._dl
         gpn = 0._dl
         gdpn = 0._dl
         gddpn = 0._dl
         gdddpn = 0._dl
         if ( params_cache%Num_Nu_Massive /= 0 ) then
            do nu_i = 1, params_cache%Nu_mass_eigenstates
               grhonu         = 0._dl
               gpinu          = 0._dl
               gpinudot       = 0._dl
               gpinudotdot    = 0._dl
               gpinudotdotdot = 0._dl
               grhormass_t = params_cache%grhormass(nu_i)/a**2
               call ThermalNuBack%rho_P(a*params_cache%nu_masses(nu_i), grhonu, gpinu )
               grhon          = grhon +grhormass_t*grhonu
               gpn            = gpn  +grhormass_t*gpinu
               gpinudot       = ThermalNuBack%pidot( a*params_cache%nu_masses(nu_i), H, gpinu )
               gdpn           = gdpn  + grhormass_t*( gpinudot - 4._dl*H*gpinu)
               gpinudotdot    = ThermalNuBack%pidotdot( a*params_cache%nu_masses(nu_i), H, Ht, gpinu, gpinudot)
               gddpn          = gddpn - 4._dl*H*grhormass_t*( gpinudot - 4._dl*H*gpinu )+ grhormass_t*( gpinudotdot - 4._dl*Ht*gpinu - 4._dl*H*gpinudot )
               gpinudotdotdot = ThermalNuBack%pidotdotdot(a*params_cache%nu_masses(nu_i),H,Ht,Ht2,gpinu,gpinudot,gpinudotdot)
               gdddpn = gdddpn + grhormass_t*(gpinudotdotdot &
                  & - 12._dl*H*gpinudotdot &
                  & + (48._dl*H**2 - 12._dl*Ht)*gpinudot &
                  & + (-64._dl*H**3 + 48._dl*H*Ht - 4._dl*Ht2)*gpinu )
            end do
         end if
         Ht4 = H**3*( &
            & grhom/6._dl + 8._dl/3._dl*grhor &
            & + grhov*( 1._dl/6._dl + 2._dl*wDE + 9._dl*wDE**2 + 18._dl*wDE**3 + 13.5_dl*wDE**4 - 1.5_dl*dwDE - 12._dl*wDE*dwDE - 27._dl*wDE**2*dwDE + 4.5_dl*dwDE**2 + 0.5_dl*d2wDE + 6._dl*wDE*d2wDE - 0.5_dl*d3wDE ) &
            & + grhon/6._dl - 2.5_dl*gpn &
            & ) &
            & + H*Ht*( &
            & -0.5_dl*grhom - 4._dl*grhor &
            & + grhov*( -0.5_dl - 4.5_dl*wDE - 13.5_dl*wDE**2 - 13.5_dl*wDE**3 + 3._dl*dwDE +13.5_dl*wDE*dwDE - 1.5_dl*d2wDE ) &
            & -0.5_dl*grhon - 4.5_dl*gpn &
            & ) &
            & + Ht2*( &
            & grhom/6._dl + 2._dl/3._dl*grhor &
            & + grhov*( 1._dl/6._dl + wDE + 1.5_dl*wDE**2 - 0.5_dl*dwDE ) &
            & + grhon/6._dl - 0.5_dl*gpn &
            & ) &
            & + gdpn*( -2._dl*Ht - 4.5_dl*H**2 ) - 2.5_dl*H*gddpn - 0.5_dl*gdddpn

         ! short vars for dnH = d^n H / da^n / H
         dH = Ht/H/a/H
         d2H = (-(Ht/H) - Ht**2/H**3 + Ht2/H**2)/a**2/H
         d3H = ((2._dl*Ht)/H + (3._dl*Ht**3)/H**5 - (3._dl*Ht2)/H**2 - (4._dl*Ht*Ht2)/H**4 + (3._dl*Ht**2 + Ht3)/H**3)/a**3/H
         d4H = ((-6._dl*Ht)/H - (15._dl*Ht**4)/H**7 + (11._dl*Ht2)/H**2 + (25._dl*Ht**2*Ht2)/H**6 + (-11._dl*Ht**2 - 6._dl*Ht3)/H**3 + (-18._dl*Ht**3 - 4._dl*Ht2**2 - 7._dl*Ht*Ht3)/H**5 + (24._dl*Ht*Ht2 + Ht4)/H**4)/a**4/H

         Ode = grhov / (grhom + grhor + grhov + grhon)
         dOde   = Ode*(-2._dl*dH + (2._dl - 3._dl*(1._dl + wDE))/a)
         d2Ode  = Ode*(-2._dl*d2H + 6._dl*dH**2 + (-8._dl*dH + 12._dl*dH*(1._dl + wDE))/a + (2._dl - 3._dl*dwDE - 9._dl*(1._dl + wDE) + 9._dl*(1._dl + wDE)**2)/a**2)
         d3Ode  = Ode*(-2._dl*d3H + 18._dl*d2H*dH - 24._dl*dH**3 + (-12._dl*d2H + 36._dl*dH**2 + 18._dl*d2H*(1._dl + wDE) - 54._dl*dH**2*(1._dl + wDE))/a + (-12._dl*dH + 18._dl*dH*dwDE + 54._dl*dH*(1._dl + wDE) - 54._dl*dH*(1._dl + wDE)**2)/a**2 + (-3._dl*d2wDE - 12._dl*dwDE - 6._dl*(1._dl + wDE) + 27._dl*dwDE*(1._dl + wDE) + 27._dl*(1._dl + wDE)**2 - 27._dl*(1._dl + wDE)**3)/a**3)
         d4Ode  = Ode*(18._dl*d2H**2 - 2._dl*d4H + 24._dl*d3H*dH - 144._dl*d2H*dH**2 + 120._dl*dH**4 + (-16._dl*d3H + 144._dl*d2H*dH - 192._dl*dH**3 + 24._dl*d3H*(1._dl + wDE) - 216._dl*d2H*dH*(1._dl + wDE) + 288._dl*dH**3*(1._dl + wDE))/a + (-24._dl*d2H + 72._dl*dH**2 + 36._dl*d2H*dwDE - 108._dl*dH**2*dwDE + 108._dl*d2H*(1._dl + wDE) - 324._dl*dH**2*(1._dl + wDE) - 108._dl*d2H*(1._dl + wDE)**2 + 324._dl*dH**2*(1._dl + wDE)**2)/a**2 + (24._dl*d2wDE*dH + 96._dl*dH*dwDE + 48._dl*dH*(1._dl + wDE) - 216._dl*dH*dwDE*(1._dl + wDE) - 216._dl*dH*(1._dl + wDE)**2 + 216._dl*dH*(1._dl + wDE)**3)/a**3 + (-15._dl*d2wDE - 3._dl*d3wDE - 6._dl*dwDE + 27._dl*dwDE**2 + 6._dl*(1._dl + wDE) + 36._dl*d2wDE*(1._dl + wDE) + 90._dl*dwDE*(1._dl + wDE) - 9._dl*(1._dl + wDE)**2 - 162._dl*dwDE*(1._dl + wDE)**2 - 54._dl*(1._dl + wDE)**3 + 81._dl*(1._dl + wDE)**4)/a**4)

         aM                = self%RPHalphaDEc_M * Ode
         daM               = self%RPHalphaDEc_M * dOde
         d2aM              = self%RPHalphaDEc_M * d2Ode
         d3aM              = self%RPHalphaDEc_M * d3Ode
         M                 = expm1(y(1))
         Mp1               = 1._dl + M
         dM                = (aM*Mp1)/a
         d2M               = ((-aM + aM**2)/a**2 + daM/a)*Mp1
         d3M               = ((2._dl*aM - 3._dl*aM**2 + aM**3)/a**3 + d2aM/a + (-2._dl*daM + 3._dl*aM*daM)/a**2)*Mp1
         d4M               = ((-6._dl*aM + 11._dl*aM**2 - 6._dl*aM**3 + aM**4)/a**4 + d3aM/a + (6._dl*daM - 14._dl*aM*daM + 6._dl*aM**2*daM)/a**3 + (-3._dl*d2aM + 4._dl*aM*d2aM + 3._dl*daM**2)/a**2)*Mp1
         aT                = self%RPHalphaDEc_T * Ode
         daT               = self%RPHalphaDEc_T * dOde
         d2aT              = self%RPHalphaDEc_T * d2Ode
         d3aT              = self%RPHalphaDEc_T * d3Ode
         d4aT              = self%RPHalphaDEc_T * d4Ode
         aK                = self%RPHalphaDEc_K * Ode
         daK               = self%RPHalphaDEc_K * dOde
         aB                = self%RPHalphaDEc_B * Ode
         daB               = self%RPHalphaDEc_B * dOde
         d2aB              = self%RPHalphaDEc_B * d2Ode
         d3aB              = self%RPHalphaDEc_B * d3Ode
         
         Omg       = M +aT*Mp1
         dOmg      = (1._dl + aT)*dM + daT*Mp1
         d2Omg     = (1._dl + aT)*d2M + 2._dl*daT*dM + d2aT*Mp1
         d3Omg     = (1._dl + aT)*d3M + 3._dl*d2M*daT + 3._dl*d2aT*dM + d3aT*Mp1
         d4Omg     = 6._dl*d2aT*d2M + (1._dl + aT)*d4M + 4._dl*d3M*daT + 4._dl*d3aT*dM + d4aT*Mp1
         
         ! store values
         self%RPHalphaDE_Omg%y(ind)       = Omg
         self%RPHalphaDE_dOmg%y(ind)      = dOmg
         self%RPHalphaDE_d2Omg%y(ind)     = d2Omg
         self%RPHalphaDE_d3Omg%y(ind)     = d3Omg
         self%RPHalphaDE_d4Omg%y(ind)     = d4Omg

         self%RPHalphaDE_c%y(ind)         = ( H**2 - Ht )*( Omg + 0.5_dl*a*dOmg ) &
            & -0.5_dl*( a*H )**2*d2Omg&
            & +0.5_dl*grhov*( 1._dl+ wDE )
         self%RPHalphaDE_dc%y(ind)      = +0.5_dl*H*grhov*( -3._dl*(1._dl +wDE)**2 + dwDE ) &
            & -Omg*( Ht2 -4._dl*H*Ht +2._dl*H**3 ) &
            & +0.5_dl*a*dOmg*( -Ht2 +H*Ht +H**3) &
            & +0.5_dl*a**2*H*d2Omg*( H**2 -3._dl*Ht ) &
            & -0.5_dl*(a*H)**3*d3Omg
         self%RPHalphaDE_Lmd%y(ind)       = -Omg*( 2._dl*Ht +H**2 ) &
            & -a*dOmg*( 2._dl*H**2 + Ht ) &
            & -( a*H )**2*d2Omg &
            & +wDE*grhov
         self%RPHalphaDE_dLmd%y(ind)      = -2._dl*Omg*( Ht2 -H*Ht -H**3 ) &
            & -a*dOmg*( +Ht2 +5._dl*H*Ht -H**3  ) &
            & -a**2*d2Omg*H*( +2._dl*H**2 +3._dl*Ht )&
            & -(a*H)**3*d3Omg &
            & +grhov*H*( dwDE -3._dl*wDE*(1._dl +wDE ))
         self%RPHalphaDE_d2Lmd%y(ind)     = -(a**4*d4Omg*H**4) + a**3*d3Omg*(-3._dl*H**4 - 6._dl*H**2*Ht) + a**2*d2Omg*(H**4 - 11._dl*H**2*Ht - 3._dl*Ht**2 - 4._dl*H*Ht2) + a*dOmg*(H**4 + 10._dl*H**2*Ht - 5._dl*Ht**2 - 6._dl*H*Ht2 - Ht3) + (-4._dl*H**4 + 2._dl*H**2*Ht + 2._dl*Ht**2 + 6._dl*H*Ht2 - 2._dl*Ht3)*Omg + grhov*(d2wDE*H**2 - 5._dl*dwDE*H**2 + dwDE*Ht + 9._dl*H**2*wDE - 9._dl*dwDE*H**2*wDE - 3._dl*Ht*wDE + 18._dl*H**2*wDE**2 - 3._dl*Ht*wDE**2 + 9._dl*H**2*wDE**3)

         ! self%RPHalphaDE_gm1%y(ind)       = 0.25_dl*( aK*Mp1*H**2 &
         !    & -2._dl*c )/(params_cache%h0_Mpc**2*a**2)
         ! self%RPHalphaDE_dgm1%y(ind)      = - 0.5_dl*( aK*Mp1*H**2 &
         !    & -2._dl*c )/(params_cache%h0_Mpc**2*a**3) &
         !    & +0.25_dl*( daK*Mp1*H**2 &
         !    & +aK*dM*H**2 &
         !    & +2._dl*aK*Mp1*Ht/a &
         !    & -4._dl*c/a -2._dl*dc/a/H )/(params_cache%h0_Mpc**2*a**2)
         self%RPHalphaDE_gm2%y(ind)       = ( +2._dl*aB*Mp1 &
            & -a*dOmg )*H/(params_cache%h0_Mpc*a)
         self%RPHalphaDE_dgm2%y(ind)      = -( +2._dl*aB*Mp1 &
            & -a*dOmg )*H/(params_cache%h0_Mpc*a**2) &
            & -( -2._dl*Mp1*( daB*H**2 &
            & + aB*Ht/a) &
            & - 2._dl*aB*H**2*dM &
            & + dOmg*( H**2 +Ht ) &
            & + a*H**2*d2Omg )/(params_cache%h0_Mpc*a*H)
         self%RPHalphaDE_d2gm2%y(ind)      = (H*(-d3Omg - 2._dl*d2Omg*dH - (4._dl*aB*dM)/a**2 + (2._dl*aB*d2M + 4._dl*daB*dM + 4._dl*aB*dH*dM)/a - d2H*dOmg + ((4._dl*aB)/a**3 + (-4._dl*daB - 4._dl*aB*dH)/a**2 + (2._dl*d2aB + 2._dl*aB*d2H + 4._dl*daB*dH)/a)*Mp1))/params_cache%h0_Mpc
         self%RPHalphaDE_d3gm2%y(ind)      = (H*(-3._dl*d2H*d2Omg - d4Omg - 3._dl*d3Omg*dH + (12._dl*aB*dM)/a**3 + (-6._dl*aB*d2M - 12._dl*daB*dM - 12._dl*aB*dH*dM)/a**2 + (2._dl*aB*d3M + 6._dl*d2M*daB + 6._dl*aB*d2M*dH + 6._dl*d2aB*dM + 6._dl*aB*d2H*dM + 12._dl*daB*dH*dM)/a - d3H*dOmg + ((-12._dl*aB)/a**4 + (12._dl*daB + 12._dl*aB*dH)/a**3 + (2._dl*d3aB + 2._dl*aB*d3H + 6._dl*d2H*daB + 6._dl*d2aB*dH)/a + (-6._dl*d2aB - 6._dl*aB*d2H - 12._dl*daB*dH)/a**2)*Mp1))/params_cache%h0_Mpc
         self%RPHalphaDE_gm3%y(ind)       = -aT*Mp1
         self%RPHalphaDE_dgm3%y(ind)      = -dM*aT -Mp1*daT
         self%RPHalphaDE_d2gm3%y(ind)     = -Mp1*d2aT - d2M*aT - 2._dl*dM*daT
         self%RPHalphaDE_d3gm3%y(ind)     = -aT*d3M - 3._dl*d2M*daT - 3._dl*d2aT*dM - d3aT*Mp1
         self%RPHalphaDE_d4gm3%y(ind)     = -6._dl*d2aT*d2M - aT*d4M - 4._dl*d3M*daT - 4._dl*d3aT*dM - d4aT*Mp1

         self%RPHalphaDE_Ode%y(ind)    = Ode
         self%RPHalphaDE_Odedot%y(ind) = dOde * a*H
         ! self%RPHalphaDE_dOde%y(ind)   = dOde
         ! self%RPHalphaDE_d2Ode%y(ind)  = d2Ode
         ! self%RPHalphaDE_d3Ode%y(ind)  = d3Ode
         ! self%RPHalphaDE_d4Ode%y(ind)  = d4Ode
         self%RPHalphaDE_M%y(ind)      = M
         self%RPHalphaDE_dM%y(ind)     = dM
         ! self%RPHalphaDE_d2M%y(ind)    = d2M
         ! self%RPHalphaDE_d3M%y(ind)    = d3M
         ! self%RPHalphaDE_d4M%y(ind)    = d4M
         self%RPHalphaDE_aK%y(ind)     = aK
         self%RPHalphaDE_daK%y(ind)    = daK
         
         

      end subroutine output

   end subroutine EFTCAMBRPHalphaDESolveDesignerEquations

   function expm1( x )

      implicit none

      real(dl), intent(in)  :: x
      real(dl) :: s, c_, e
      real(dl) :: expm1

      if (abs(x) .le. 0.6931471805599453094172321214581765680755d0) then
         !
         !   degree 13 minimax polynomial approximation
         !
         expm1=x*(0.99999999999999999301d0 &
               +x*(0.49999999999999957237d0 &
               +x*(0.16666666666666769423d0 &
               +x*(0.041666666666691966524d0 &
               +x*(0.0083333333333092372395d0 &
               +x*(0.0013888888884629134308d0 &
               +x*(0.00019841269861419530922d0 &
               +x*(0.000024801590367084381006d0 &
               +x*(2.7557312010044843920d-6 &
               +x*(2.7556249037054640535d-7 &
               +x*(2.5053112146542602098d-8 &
               +x*(2.1055737595306817199d-9 &
               +x*(1.6058956553927382096d-10)))))))))))))
      else
         !
         !   standard approach using "sinh"
         !
         s=sinh(0.5d0*x)
         c_=sqrt(s*s+1.d0)
         e=s+c_
         expm1=2.d0*e*s
      end if

   end function expm1
   
   ! ---------------------------------------------------------------------------------------------
   !> Subroutine that computes the value of the background EFT functions at a given time.
   subroutine EFTCAMBRPHalphaDEBackgroundEFTFunctions( self, a, eft_par_cache, eft_cache )

      implicit none

      class(EFTCAMB_RPH_alphaDE)                           :: self          !< the base class
      real(dl), intent(in)                         :: a             !< the input scale factor
      type(TEFTCAMB_parameter_cache), intent(inout) :: eft_par_cache !< the EFTCAMB parameter cache that contains all the physical parameters.
      type(TEFTCAMB_timestep_cache ), intent(inout) :: eft_cache     !< the EFTCAMB timestep cache that contains all the physical values.

      ! real(dl) :: M, Mp1, dM, d2M, d3M, d4M, aT, daT, d2aT, d3aT, d4aT, aK, daK, aB, daB, d2aB, d3aB, aM, daM, d2aM, d3aM
      real(dl) :: Omg, dOmg, d2Omg, d3Omg, d4Omg
      real(dl) :: H, Ht, Ht2, Ht3
      ! real(dl) :: Ode, dOde, d2Ode, d3Ode, d4Ode
      real(dl) :: wDE, dwDE, d2wDE, d3wDE
      real(dl) :: x

      ! ! precompute some functions:
      H = eft_cache%adotoa
      Ht = eft_cache%Hdot
      Ht2 = eft_cache%Hdotdot
      Ht3 = eft_cache%Hdotdotdot

      ! ! dnH = d^n H / da^n / H
      ! ! dH = Ht/H/a/H
      ! ! d2H = (-(Ht/H) - Ht**2/H**3 + Ht2/H**2)/a**2/H
      ! ! d3H = ((2._dl*Ht)/H + (3._dl*Ht**3)/H**5 - (3._dl*Ht2)/H**2 - (4._dl*Ht*Ht2)/H**4 + (3._dl*Ht**2 + Ht3)/H**3)/a**3/H
      ! ! d4H = ((-6._dl*Ht)/H - (15._dl*Ht**4)/H**7 + (11._dl*Ht2)/H**2 + (25._dl*Ht**2*Ht2)/H**6 + (-11._dl*Ht**2 - 6._dl*Ht3)/H**3 + (-18._dl*Ht**3 - 4._dl*Ht2**2 - 7._dl*Ht*Ht3)/H**5 + (24._dl*Ht*Ht2 + Ht4)/H**4)/a**4/H

      wDE   = self%RPHalphaDE_wDE%value(a)
      dwDE  = a * self%RPHalphaDE_wDE%first_derivative(a)
      d2wDE = a**2 * self%RPHalphaDE_wDE%second_derivative(a)
      d3wDE = a**3 * self%RPHalphaDE_wDE%third_derivative(a)

      ! x = log(a)
      ! if (x .le. self%x_initial) then
      !    M     = 0._dl
      !    dM  = 0._dl
      !    d2M = 0._dl
      !    d3M = 0._dl
      !    d4M = 0._dl
      !    Ode   = 0._dl
      !    dOde  = 0._dl
      !    d2Ode = 0._dl
      !    d3Ode = 0._dl
      !    d4Ode = 0._dl
      ! else
      !    M     = expm1(self%RPHalphaDE_log1pM%value(x))
      !    dM    = self%RPHalphaDE_dM%value(x)
      !    d2M   = self%RPHalphaDE_d2M%value(x)
      !    d3M   = self%RPHalphaDE_d3M%value(x)
      !    d4M   = self%RPHalphaDE_d4M%value(x)
      !    Ode   = self%RPHalphaDE_Ode%value(x)
      !    dOde  = self%RPHalphaDE_dOde%value(x)
      !    d2Ode = self%RPHalphaDE_d2Ode%value(x)
      !    d3Ode = self%RPHalphaDE_d3Ode%value(x)
      !    d4Ode = self%RPHalphaDE_d4Ode%value(x)
      !    ! Ode = eft_cache%grhov_t / H**2 / 3._dl
      !    ! dOde = Ode*(-2._dl*dH + (2._dl - 3._dl*(1._dl + wDE))/a)
      !    ! d2Ode = Ode*(-2._dl*d2H + 6._dl*dH**2 + (-8._dl*dH + 12._dl*dH*(1._dl + wDE))/a + (2._dl - 3._dl*dwDE - 9._dl*(1._dl + wDE) + 9._dl*(1._dl + wDE)**2)/a**2)
      !    ! d3Ode = Ode*(-2._dl*d3H + 18._dl*d2H*dH - 24._dl*dH**3 + (-12._dl*d2H + 36._dl*dH**2 + 18._dl*d2H*(1._dl + wDE) - 54._dl*dH**2*(1._dl + wDE))/a + (-12._dl*dH + 18._dl*dH*dwDE + 54._dl*dH*(1._dl + wDE) - 54._dl*dH*(1._dl + wDE)**2)/a**2 + (-3._dl*d2wDE - 12._dl*dwDE - 6._dl*(1._dl + wDE) + 27._dl*dwDE*(1._dl + wDE) + 27._dl*(1._dl + wDE)**2 - 27._dl*(1._dl + wDE)**3)/a**3)
      !    ! d4Ode = Ode*(18._dl*d2H**2 - 2._dl*d4H + 24._dl*d3H*dH - 144._dl*d2H*dH**2 + 120._dl*dH**4 + (-16._dl*d3H + 144._dl*d2H*dH - 192._dl*dH**3 + 24._dl*d3H*(1._dl + wDE) - 216._dl*d2H*dH*(1._dl + wDE) + 288._dl*dH**3*(1._dl + wDE))/a + (-24._dl*d2H + 72._dl*dH**2 + 36._dl*d2H*dwDE - 108._dl*dH**2*dwDE + 108._dl*d2H*(1._dl + wDE) - 324._dl*dH**2*(1._dl + wDE) - 108._dl*d2H*(1._dl + wDE)**2 + 324._dl*dH**2*(1._dl + wDE)**2)/a**2 + (24._dl*d2wDE*dH + 96._dl*dH*dwDE + 48._dl*dH*(1._dl + wDE) - 216._dl*dH*dwDE*(1._dl + wDE) - 216._dl*dH*(1._dl + wDE)**2 + 216._dl*dH*(1._dl + wDE)**3)/a**3 + (-15._dl*d2wDE - 3._dl*d3wDE - 6._dl*dwDE + 27._dl*dwDE**2 + 6._dl*(1._dl + wDE) + 36._dl*d2wDE*(1._dl + wDE) + 90._dl*dwDE*(1._dl + wDE) - 9._dl*(1._dl + wDE)**2 - 162._dl*dwDE*(1._dl + wDE)**2 - 54._dl*(1._dl + wDE)**3 + 81._dl*(1._dl + wDE)**4)/a**4)
      ! end if      

      ! ! aM                = self%RPHalphaDEc_M * Ode
      ! ! daM               = self%RPHalphaDEc_M * dOde
      ! ! d2aM              = self%RPHalphaDEc_M * d2Ode
      ! ! d3aM              = self%RPHalphaDEc_M * d3Ode
      ! Mp1               = 1._dl + M
      ! ! dM                = (aM*Mp1)/a
      ! ! d2M               = ((-aM + aM**2)/a**2 + daM/a)*Mp1
      ! ! d3M               = ((2._dl*aM - 3._dl*aM**2 + aM**3)/a**3 + d2aM/a + (-2._dl*daM + 3._dl*aM*daM)/a**2)*Mp1
      ! ! d4M               = ((-6._dl*aM + 11._dl*aM**2 - 6._dl*aM**3 + aM**4)/a**4 + d3aM/a + (6._dl*daM - 14._dl*aM*daM + 6._dl*aM**2*daM)/a**3 + (-3._dl*d2aM + 4._dl*aM*d2aM + 3._dl*daM**2)/a**2)*Mp1
      ! aT                = self%RPHalphaDEc_T * Ode
      ! daT               = self%RPHalphaDEc_T * dOde
      ! d2aT              = self%RPHalphaDEc_T * d2Ode
      ! d3aT              = self%RPHalphaDEc_T * d3Ode
      ! d4aT              = self%RPHalphaDEc_T * d4Ode

      ! ! compute the EFT functions:
      ! Omg       = M +aT*Mp1
      ! dOmg      = (1._dl + aT)*dM + daT*Mp1
      ! d2Omg     = (1._dl + aT)*d2M + 2._dl*daT*dM + d2aT*Mp1
      ! d3Omg     = (1._dl + aT)*d3M + 3._dl*d2M*daT + 3._dl*d2aT*dM + d3aT*Mp1
      ! d4Omg     = 6._dl*d2aT*d2M + (1._dl + aT)*d4M + 4._dl*d3M*daT + 4._dl*d3aT*dM + d4aT*Mp1

      ! eft_cache%EFTOmegaV    = Omg
      ! eft_cache%EFTOmegaP    = dOmg
      ! eft_cache%EFTOmegaPP   = d2Omg
      ! eft_cache%EFTOmegaPPP  = d3Omg
      ! eft_cache%EFTOmegaPPPP = d4Omg

      ! eft_cache%EFTc         = ( eft_cache%adotoa**2 - eft_cache%Hdot )*( eft_cache%EFTOmegaV + 0.5_dl*a*eft_cache%EFTOmegaP ) &
      ! & -0.5_dl*( a*eft_cache%adotoa )**2*eft_cache%EFTOmegaPP&
      ! & +0.5_dl*eft_cache%grhov_t*( 1._dl+ self%RPHalphaDE_wDE%value(a) )
      ! eft_cache%EFTLambda    = -eft_cache%EFTOmegaV*( 2._dl*eft_cache%Hdot +eft_cache%adotoa**2 ) &
      ! & -a*eft_cache%EFTOmegaP*( 2._dl*eft_cache%adotoa**2 + eft_cache%Hdot ) &
      ! & -( a*eft_cache%adotoa )**2*eft_cache%EFTOmegaPP &
      ! & +self%RPHalphaDE_wDE%value(a)*eft_cache%grhov_t
      ! eft_cache%EFTcdot      = +0.5_dl*eft_cache%adotoa*eft_cache%grhov_t*( -3._dl*(1._dl +self%RPHalphaDE_wDE%value(a))**2 + a*self%RPHalphaDE_wDE%first_derivative(a) ) &
      ! & -eft_cache%EFTOmegaV*( eft_cache%Hdotdot -4._dl*eft_cache%adotoa*eft_cache%Hdot +2._dl*eft_cache%adotoa**3 ) &
      ! & +0.5_dl*a*eft_cache%EFTOmegaP*( -eft_cache%Hdotdot +eft_cache%adotoa*eft_cache%Hdot +eft_cache%adotoa**3) &
      ! & +0.5_dl*a**2*eft_cache%adotoa*eft_cache%EFTOmegaPP*( eft_cache%adotoa**2 -3._dl*eft_cache%Hdot ) &
      ! & -0.5_dl*(a*eft_cache%adotoa)**3*eft_cache%EFTOmegaPPP
      ! eft_cache%EFTLambdadot = -2._dl*eft_cache%EFTOmegaV*( eft_cache%Hdotdot -eft_cache%adotoa*eft_cache%Hdot -eft_cache%adotoa**3 ) &
      ! & -a*eft_cache%EFTOmegaP*( +eft_cache%Hdotdot +5._dl*eft_cache%adotoa*eft_cache%Hdot -eft_cache%adotoa**3  ) &
      ! & -a**2*eft_cache%EFTOmegaPP*eft_cache%adotoa*( +2._dl*eft_cache%adotoa**2 +3._dl*eft_cache%Hdot )&
      ! & -(a*eft_cache%adotoa)**3*eft_cache%EFTOmegaPPP &
      ! & +eft_cache%grhov_t*eft_cache%adotoa*( a*self%RPHalphaDE_wDE%first_derivative(a) -3._dl*self%RPHalphaDE_wDE%value(a)*(1._dl +self%RPHalphaDE_wDE%value(a) ))
      ! ! required for positivity
      ! eft_cache%EFTLambdadotdot = -(a**4*d4Omg*H**4) + a**3*d3Omg*(-3._dl*H**4 - 6._dl*H**2*Ht) + a**2*d2Omg*(H**4 - 11._dl*H**2*Ht - 3._dl*Ht**2 - 4._dl*H*Ht2) + a*dOmg*(H**4 + 10._dl*H**2*Ht - 5._dl*Ht**2 - 6._dl*H*Ht2 - Ht3) + (-4._dl*H**4 + 2._dl*H**2*Ht + 2._dl*Ht**2 + 6._dl*H*Ht2 - 2._dl*Ht3)*Omg + eft_cache%grhov_t*(d2wDE*H**2 - 5._dl*dwDE*H**2 + dwDE*Ht + 9._dl*H**2*wDE - 9._dl*dwDE*H**2*wDE - 3._dl*Ht*wDE + 18._dl*H**2*wDE**2 - 3._dl*Ht*wDE**2 + 9._dl*H**2*wDE**3)
      x = log(a)
      if (x .le. self%x_initial) then
         Omg                       = 0._dl
         dOmg                      = 0._dl
         d2Omg                     = 0._dl
         d3Omg                     = 0._dl
         d4Omg                     = 0._dl
         eft_cache%EFTLambda       = wDE*eft_cache%grhov_t
         eft_cache%EFTLambdadot    = eft_cache%grhov_t*H*( dwDE -3._dl*wDE*(1._dl +wDE ))
         eft_cache%EFTLambdadotdot = eft_cache%grhov_t*(d2wDE*H**2 - 5._dl*dwDE*H**2 + dwDE*Ht + 9._dl*H**2*wDE - 9._dl*dwDE*H**2*wDE - 3._dl*Ht*wDE + 18._dl*H**2*wDE**2 - 3._dl*Ht*wDE**2 + 9._dl*H**2*wDE**3)
         eft_cache%EFTc            = +0.5_dl*eft_cache%grhov_t*( 1._dl+ wDE )
         eft_cache%EFTcdot         = +0.5_dl*H*eft_cache%grhov_t*( -3._dl*(1._dl +wDE)**2 + dwDE )
      else
         Omg       = self%RPHalphaDE_Omg%value(x)
         dOmg       = self%RPHalphaDE_dOmg%value(x)
         d2Omg      = self%RPHalphaDE_d2Omg%value(x)
         d3Omg     = self%RPHalphaDE_d3Omg%value(x)
         d4Omg    = self%RPHalphaDE_d4Omg%value(x)
         eft_cache%EFTLambda       = self%RPHalphaDE_Lmd%value(x)
         eft_cache%EFTLambdadot    = self%RPHalphaDE_dLmd%value(x)
         eft_cache%EFTLambdadotdot = self%RPHalphaDE_d2Lmd%value(x)
         eft_cache%EFTc            = self%RPHalphaDE_c%value(x)
         eft_cache%EFTcdot         = self%RPHalphaDE_dc%value(x)
      end if
      
      eft_cache%EFTOmegaV    = Omg
      eft_cache%EFTOmegaP    = dOmg
      eft_cache%EFTOmegaPP   = d2Omg
      eft_cache%EFTOmegaPPP  = d3Omg
      eft_cache%EFTOmegaPPPP = d4Omg
      
      ! eft_cache%EFTc         = ( H**2 - Ht )*( Omg + 0.5_dl*a*dOmg ) &
      !    & -0.5_dl*( a*H )**2*d2Omg&
      !    & +0.5_dl*eft_cache%grhov_t*( 1._dl+ wDE )
      ! eft_cache%EFTcdot      = +0.5_dl*H*eft_cache%grhov_t*( -3._dl*(1._dl +wDE)**2 + dwDE ) &
      !    & -Omg*( Ht2 -4._dl*H*Ht +2._dl*H**3 ) &
      !    & +0.5_dl*a*dOmg*( -Ht2 +H*Ht +H**3) &
      !    & +0.5_dl*a**2*H*d2Omg*( H**2 -3._dl*Ht ) &
      !    & -0.5_dl*(a*H)**3*d3Omg
      ! eft_cache%EFTLambda    = -eft_cache%EFTOmegaV*( 2._dl*eft_cache%Hdot +eft_cache%adotoa**2 ) &
      ! & -a*eft_cache%EFTOmegaP*( 2._dl*eft_cache%adotoa**2 + eft_cache%Hdot ) &
      ! & -( a*eft_cache%adotoa )**2*eft_cache%EFTOmegaPP &
      ! & +self%RPHalphaDE_wDE%value(a)*eft_cache%grhov_t
      ! eft_cache%EFTLambdadot = -2._dl*eft_cache%EFTOmegaV*( eft_cache%Hdotdot -eft_cache%adotoa*eft_cache%Hdot -eft_cache%adotoa**3 ) &
      ! & -a*eft_cache%EFTOmegaP*( +eft_cache%Hdotdot +5._dl*eft_cache%adotoa*eft_cache%Hdot -eft_cache%adotoa**3  ) &
      ! & -a**2*eft_cache%EFTOmegaPP*eft_cache%adotoa*( +2._dl*eft_cache%adotoa**2 +3._dl*eft_cache%Hdot )&
      ! & -(a*eft_cache%adotoa)**3*eft_cache%EFTOmegaPPP &
      ! & +eft_cache%grhov_t*eft_cache%adotoa*( a*self%RPHalphaDE_wDE%first_derivative(a) -3._dl*self%RPHalphaDE_wDE%value(a)*(1._dl +self%RPHalphaDE_wDE%value(a) ))
      ! eft_cache%EFTLambdadotdot = -(a**4*d4Omg*H**4) + a**3*d3Omg*(-3._dl*H**4 - 6._dl*H**2*Ht) + a**2*d2Omg*(H**4 - 11._dl*H**2*Ht - 3._dl*Ht**2 - 4._dl*H*Ht2) + a*dOmg*(H**4 + 10._dl*H**2*Ht - 5._dl*Ht**2 - 6._dl*H*Ht2 - Ht3) + (-4._dl*H**4 + 2._dl*H**2*Ht + 2._dl*Ht**2 + 6._dl*H*Ht2 - 2._dl*Ht3)*Omg + eft_cache%grhov_t*(d2wDE*H**2 - 5._dl*dwDE*H**2 + dwDE*Ht + 9._dl*H**2*wDE - 9._dl*dwDE*H**2*wDE - 3._dl*Ht*wDE + 18._dl*H**2*wDE**2 - 3._dl*Ht*wDE**2 + 9._dl*H**2*wDE**3)


   end subroutine EFTCAMBRPHalphaDEBackgroundEFTFunctions

   ! ---------------------------------------------------------------------------------------------
   !> Subroutine that computes the value of the background EFT functions at a given time.
   subroutine EFTCAMBRPHalphaDESecondOrderEFTFunctions( self, a, eft_par_cache, eft_cache )

      implicit none

      class(EFTCAMB_RPH_alphaDE)                           :: self          !< the base class
      real(dl), intent(in)                         :: a             !< the input scale factor
      type(TEFTCAMB_parameter_cache), intent(inout) :: eft_par_cache !< the EFTCAMB parameter cache that contains all the physical parameters.
      type(TEFTCAMB_timestep_cache ), intent(inout) :: eft_cache     !< the EFTCAMB timestep cache that contains all the physical values.

      ! real(dl) :: M, Mp1, dM, d2M, d3M, d4M, aT, daT, d2aT, d3aT, d4aT, aK, daK, aB, daB, d2aB, d3aB, aM, daM, d2aM, d3aM
      ! real(dl) :: Omg, dOmg, d2Omg, d3Omg, d4Omg
      real(dl) :: H, dH, d2H, d3H, d4H, ca, Ht, Ht2, Ht3, Ht4
      ! real(dl) :: Ode, dOde, d2Ode, d3Ode, d4Ode
      ! real(dl) :: wDE, dwDE, d2wDE, d3wDE, d4wDE
      real(dl) :: x, z
      real(dl) :: Mp1, dM, aK, daK
      real(dl) :: regz

      ! ! precompute some functions:
      H = eft_cache%adotoa
      Ht = eft_cache%Hdot
      ! Ht2 = eft_cache%Hdotdot
      ! Ht3 = eft_cache%Hdotdotdot
      ! ! Ht4 = eft_cache%Hdotdotdotdot

      ! ! dnH = d^n H / da^n / H
      ! dH = Ht/H/a/H
      ! d2H = (-(Ht/H) - Ht**2/H**3 + Ht2/H**2)/a**2/H
      ! d3H = ((2._dl*Ht)/H + (3._dl*Ht**3)/H**5 - (3._dl*Ht2)/H**2 - (4._dl*Ht*Ht2)/H**4 + (3._dl*Ht**2 + Ht3)/H**3)/a**3/H
      ! ! d4H = ((-6._dl*Ht)/H - (15._dl*Ht**4)/H**7 + (11._dl*Ht2)/H**2 + (25._dl*Ht**2*Ht2)/H**6 + (-11._dl*Ht**2 - 6._dl*Ht3)/H**3 + (-18._dl*Ht**3 - 4._dl*Ht2**2 - 7._dl*Ht*Ht3)/H**5 + (24._dl*Ht*Ht2 + Ht4)/H**4)/a**4/H

      ! wDE   = self%RPHalphaDE_wDE%value(a)
      ! dwDE  = a * self%RPHalphaDE_wDE%first_derivative(a)
      ! d2wDE = a**2 * self%RPHalphaDE_wDE%second_derivative(a)
      ! d3wDE = a**3 * self%RPHalphaDE_wDE%third_derivative(a)

      ! x = log(a)
      ! if (x .le. self%x_initial) then
      !    M     = 0._dl
      !    dM  = 0._dl
      !    d2M = 0._dl
      !    d3M = 0._dl
      !    d4M = 0._dl
      !    Ode   = 0._dl
      !    dOde  = 0._dl
      !    d2Ode = 0._dl
      !    d3Ode = 0._dl
      !    d4Ode = 0._dl
      ! else
      !    M     = expm1(self%RPHalphaDE_log1pM%value(x))
      !    dM    = self%RPHalphaDE_dM%value(x)
      !    d2M   = self%RPHalphaDE_d2M%value(x)
      !    d3M   = self%RPHalphaDE_d3M%value(x)
      !    d4M   = self%RPHalphaDE_d4M%value(x)
      !    Ode   = self%RPHalphaDE_Ode%value(x)
      !    dOde  = self%RPHalphaDE_dOde%value(x)
      !    d2Ode = self%RPHalphaDE_d2Ode%value(x)
      !    d3Ode = self%RPHalphaDE_d3Ode%value(x)
      !    d4Ode = self%RPHalphaDE_d4Ode%value(x)
      ! !   Ode = eft_cache%grhov_t / H**2 / 3._dl
      ! !    dOde = Ode*(-2._dl*dH + (2._dl - 3._dl*(1._dl + wDE))/a)
      ! !    d2Ode = Ode*(-2._dl*d2H + 6._dl*dH**2 + (-8._dl*dH + 12._dl*dH*(1._dl + wDE))/a + (2._dl - 3._dl*dwDE - 9._dl*(1._dl + wDE) + 9._dl*(1._dl + wDE)**2)/a**2)
      ! !    d3Ode = Ode*(-2._dl*d3H + 18._dl*d2H*dH - 24._dl*dH**3 + (-12._dl*d2H + 36._dl*dH**2 + 18._dl*d2H*(1._dl + wDE) - 54._dl*dH**2*(1._dl + wDE))/a + (-12._dl*dH + 18._dl*dH*dwDE + 54._dl*dH*(1._dl + wDE) - 54._dl*dH*(1._dl + wDE)**2)/a**2 + (-3._dl*d2wDE - 12._dl*dwDE - 6._dl*(1._dl + wDE) + 27._dl*dwDE*(1._dl + wDE) + 27._dl*(1._dl + wDE)**2 - 27._dl*(1._dl + wDE)**3)/a**3)
      ! !    d4Ode = Ode*(18._dl*d2H**2 - 2._dl*d4H + 24._dl*d3H*dH - 144._dl*d2H*dH**2 + 120._dl*dH**4 + (-16._dl*d3H + 144._dl*d2H*dH - 192._dl*dH**3 + 24._dl*d3H*(1._dl + wDE) - 216._dl*d2H*dH*(1._dl + wDE) + 288._dl*dH**3*(1._dl + wDE))/a + (-24._dl*d2H + 72._dl*dH**2 + 36._dl*d2H*dwDE - 108._dl*dH**2*dwDE + 108._dl*d2H*(1._dl + wDE) - 324._dl*dH**2*(1._dl + wDE) - 108._dl*d2H*(1._dl + wDE)**2 + 324._dl*dH**2*(1._dl + wDE)**2)/a**2 + (24._dl*d2wDE*dH + 96._dl*dH*dwDE + 48._dl*dH*(1._dl + wDE) - 216._dl*dH*dwDE*(1._dl + wDE) - 216._dl*dH*(1._dl + wDE)**2 + 216._dl*dH*(1._dl + wDE)**3)/a**3 + (-15._dl*d2wDE - 3._dl*d3wDE - 6._dl*dwDE + 27._dl*dwDE**2 + 6._dl*(1._dl + wDE) + 36._dl*d2wDE*(1._dl + wDE) + 90._dl*dwDE*(1._dl + wDE) - 9._dl*(1._dl + wDE)**2 - 162._dl*dwDE*(1._dl + wDE)**2 - 54._dl*(1._dl + wDE)**3 + 81._dl*(1._dl + wDE)**4)/a**4)
      ! end if

      ! ! aM                = self%RPHalphaDEc_M * Ode
      ! ! daM               = self%RPHalphaDEc_M * dOde
      ! ! d2aM              = self%RPHalphaDEc_M * d2Ode
      ! ! d3aM              = self%RPHalphaDEc_M * d3Ode
      ! Mp1               = 1._dl + M
      ! ! dM                = (aM*Mp1)/a
      ! ! d2M               = ((-aM + aM**2)/a**2 + daM/a)*Mp1
      ! ! d3M               = ((2._dl*aM - 3._dl*aM**2 + aM**3)/a**3 + d2aM/a + (-2._dl*daM + 3._dl*aM*daM)/a**2)*Mp1
      ! ! d4M               = ((-6._dl*aM + 11._dl*aM**2 - 6._dl*aM**3 + aM**4)/a**4 + d3aM/a + (6._dl*daM - 14._dl*aM*daM + 6._dl*aM**2*daM)/a**3 + (-3._dl*d2aM + 4._dl*aM*d2aM + 3._dl*daM**2)/a**2)*Mp1
      ! aT                = self%RPHalphaDEc_T * Ode
      ! daT               = self%RPHalphaDEc_T * dOde
      ! d2aT              = self%RPHalphaDEc_T * d2Ode
      ! d3aT              = self%RPHalphaDEc_T * d3Ode
      ! d4aT              = self%RPHalphaDEc_T * d4Ode
      ! aK                = self%RPHalphaDEc_K * Ode
      ! daK               = self%RPHalphaDEc_K * dOde
      ! aB                = self%RPHalphaDEc_B * Ode
      ! daB               = self%RPHalphaDEc_B * dOde
      ! d2aB              = self%RPHalphaDEc_B * d2Ode
      ! d3aB              = self%RPHalphaDEc_B * d3Ode

      ! ! compute the EFT functions:
      ! Omg       = eft_cache%EFTOmegaV
      ! dOmg      = eft_cache%EFTOmegaP
      ! d2Omg     = eft_cache%EFTOmegaPP
      ! d3Omg     = eft_cache%EFTOmegaPPP
      ! d4Omg     = eft_cache%EFTOmegaPPPP

      ! ! compute the EFT functions:
      ! eft_cache%EFTGamma1V    = 0.25_dl*( aK*Mp1*eft_cache%adotoa**2 &
      ! & -2._dl*eft_cache%EFTc )/(eft_par_cache%h0_Mpc**2*a**2)
      ! eft_cache%EFTGamma1P    = - 0.5_dl*( aK*Mp1*eft_cache%adotoa**2 &
      ! & -2._dl*eft_cache%EFTc )/(eft_par_cache%h0_Mpc**2*a**3) &
      ! & +0.25_dl*( daK*Mp1*eft_cache%adotoa**2 &
      ! & +aK*dM*eft_cache%adotoa**2 &
      ! & +2._dl*aK*Mp1*eft_cache%Hdot/a &
      ! & -4._dl*eft_cache%EFTc/a -2._dl*eft_cache%EFTcdot/a/eft_cache%adotoa )/(eft_par_cache%h0_Mpc**2*a**2)
      ! eft_cache%EFTGamma2V    = ( +2._dl*aB*Mp1 &
      ! & -a*eft_cache%EFTOmegaP )*eft_cache%adotoa/(eft_par_cache%h0_Mpc*a)
      ! eft_cache%EFTGamma2P    = -( +2._dl*aB*Mp1 &
      ! & -a*eft_cache%EFTOmegaP )*eft_cache%adotoa/(eft_par_cache%h0_Mpc*a**2) &
      ! & -( -2._dl*Mp1*( daB*eft_cache%adotoa**2 &
      ! & + aB*eft_cache%Hdot/a) &
      ! & - 2._dl*aB*eft_cache%adotoa**2*dM &
      ! & + eft_cache%EFTOmegaP*( eft_cache%adotoa**2 +eft_cache%Hdot ) &
      ! & + a*eft_cache%adotoa**2*eft_cache%EFTOmegaPP )/(eft_par_cache%h0_Mpc*a*eft_cache%adotoa)
      ! eft_cache%EFTGamma2PP   = (H*(-d3Omg - 2._dl*d2Omg*dH - (4._dl*aB*dM)/a**2 + (2._dl*aB*d2M + 4._dl*daB*dM + 4._dl*aB*dH*dM)/a - d2H*dOmg + ((4._dl*aB)/a**3 + (-4._dl*daB - 4._dl*aB*dH)/a**2 + (2._dl*d2aB + 2._dl*aB*d2H + 4._dl*daB*dH)/a)*Mp1))/eft_par_cache%h0_Mpc
      ! eft_cache%EFTGamma2PPP  = (H*(-3._dl*d2H*d2Omg - d4Omg - 3._dl*d3Omg*dH + (12._dl*aB*dM)/a**3 + (-6._dl*aB*d2M - 12._dl*daB*dM - 12._dl*aB*dH*dM)/a**2 + (2._dl*aB*d3M + 6._dl*d2M*daB + 6._dl*aB*d2M*dH + 6._dl*d2aB*dM + 6._dl*aB*d2H*dM + 12._dl*daB*dH*dM)/a - d3H*dOmg + ((-12._dl*aB)/a**4 + (12._dl*daB + 12._dl*aB*dH)/a**3 + (2._dl*d3aB + 2._dl*aB*d3H + 6._dl*d2H*daB + 6._dl*d2aB*dH)/a + (-6._dl*d2aB - 6._dl*aB*d2H - 12._dl*daB*dH)/a**2)*Mp1))/eft_par_cache%h0_Mpc
      ! eft_cache%EFTGamma3V    = -aT*Mp1
      ! eft_cache%EFTGamma3P    = -dM*aT -Mp1*daT
      ! eft_cache%EFTGamma3PP   = -Mp1*d2aT - d2M*aT - 2._dl*dM*daT
      ! eft_cache%EFTGamma3PPP  = -aT*d3M - 3._dl*d2M*daT - 3._dl*d2aT*dM - d3aT*Mp1
      ! eft_cache%EFTGamma3PPPP = -6._dl*d2aT*d2M - aT*d4M - 4._dl*d3M*daT - 4._dl*d3aT*dM - d4aT*Mp1

      x = log(a)
      ! z = 1._dl/a + 1._dl
      ! regz = 0.5_dl + 0.5_dl*tanh(0.25_dl*(z - 10._dl))
      if (x .le. self%x_initial) then
         aK                      = 0._dl
         daK                     = 0._dl
         Mp1                     = 1._dl
         dM                      = 0._dl
         eft_cache%EFTGamma2V    = 0._dl
         eft_cache%EFTGamma2P    = 0._dl
         eft_cache%EFTGamma2PP   = 0._dl
         eft_cache%EFTGamma2PPP  = 0._dl
         eft_cache%EFTGamma3V    = 0._dl
         eft_cache%EFTGamma3P    = 0._dl
         eft_cache%EFTGamma3PP   = 0._dl
         eft_cache%EFTGamma3PPP  = 0._dl
         eft_cache%EFTGamma3PPPP = 0._dl
         eft_cache%alphaK        = 0._dl
         eft_cache%alphaKdot     = 0._dl
         eft_cache%alphaB        = 0._dl
         eft_cache%alphaBdot     = 0._dl
         eft_cache%alphaM        = 0._dl
         eft_cache%alphaMdot     = 0._dl
         eft_cache%alphaT        = 0._dl
         eft_cache%alphaTdot     = 0._dl
      else
         aK                      = self%RPHalphaDE_aK%value(x)
         daK                     = self%RPHalphaDE_daK%value(x)
         Mp1                     = self%RPHalphaDE_M%value(x) + 1._dl
         dM                      = self%RPHalphaDE_dM%value(x)
         eft_cache%EFTGamma2V    = self%RPHalphaDE_gm2%value(x)
         eft_cache%EFTGamma2P    = self%RPHalphaDE_dgm2%value(x)
         eft_cache%EFTGamma2PP   = self%RPHalphaDE_d2gm2%value(x)
         eft_cache%EFTGamma2PPP  = self%RPHalphaDE_d3gm2%value(x)
         eft_cache%EFTGamma3V    = self%RPHalphaDE_gm3%value(x)
         eft_cache%EFTGamma3P    = self%RPHalphaDE_dgm3%value(x)
         eft_cache%EFTGamma3PP   = self%RPHalphaDE_d2gm3%value(x)
         eft_cache%EFTGamma3PPP  = self%RPHalphaDE_d3gm3%value(x)
         eft_cache%EFTGamma3PPPP = self%RPHalphaDE_d4gm3%value(x)
         eft_cache%alphaK        = self%RPHalphaDEc_K * self%RPHalphaDE_Ode%value(x)
         eft_cache%alphaKdot     = self%RPHalphaDEc_K * self%RPHalphaDE_Odedot%value(x)
         eft_cache%alphaB        = self%RPHalphaDEc_B * self%RPHalphaDE_Ode%value(x)
         eft_cache%alphaBdot     = self%RPHalphaDEc_B * self%RPHalphaDE_Odedot%value(x)
         eft_cache%alphaT        = self%RPHalphaDEc_T * self%RPHalphaDE_Ode%value(x)
         eft_cache%alphaTdot     = self%RPHalphaDEc_T * self%RPHalphaDE_Odedot%value(x)
         eft_cache%alphaM        = self%RPHalphaDEc_M * self%RPHalphaDE_Ode%value(x)
         eft_cache%alphaMdot     = self%RPHalphaDEc_M * self%RPHalphaDE_Odedot%value(x)
      end if
      eft_cache%EFTGamma1V    = 0.25_dl*( aK*Mp1*H**2 &
         & -2._dl*eft_cache%EFTc )/(eft_par_cache%h0_Mpc**2*a**2)
      eft_cache%EFTGamma1P    = - 0.5_dl*( aK*Mp1*H**2 )/(eft_par_cache%h0_Mpc**2*a**3) &
         & +0.25_dl*( daK*Mp1*H**2 &
         & +aK*dM*H**2 &
         & +2._dl*aK*Mp1*Ht/a &
         & -2._dl*eft_cache%EFTcdot/a/H )/(eft_par_cache%h0_Mpc**2*a**2)
      eft_cache%EFTGamma4V    = -eft_cache%EFTGamma3V
      eft_cache%EFTGamma4P    = -eft_cache%EFTGamma3P
      eft_cache%EFTGamma4PP   = -eft_cache%EFTGamma3PP
      eft_cache%EFTGamma5V    = +0.5_dl*eft_cache%EFTGamma3V
      eft_cache%EFTGamma5P    = +0.5_dl*eft_cache%EFTGamma3P
      eft_cache%EFTGamma6V    = 0._dl
      eft_cache%EFTGamma6P    = 0._dl

    !   if (a>1e-3) then
    !      write (*,*) a, Omg, eft_cache%EFTGamma1V/(H/a)**2*eft_par_cache%h0_Mpc**2, eft_cache%EFTGamma2V, eft_cache%EFTGamma3V
    !   end if

   end subroutine EFTCAMBRPHalphaDESecondOrderEFTFunctions

   ! ---------------------------------------------------------------------------------------------
   !> Function that computes dtauda = 1/sqrt(a^2H^2).For pure EFT std this has to be overridden
   !! for performance reasons.
   function EFTCAMBRPHalphaDEComputeDtauda( self, a, eft_par_cache, eft_cache )

      implicit none

      class(EFTCAMB_RPH_alphaDE)                           :: self          !< the base class
      real(dl), intent(in)                         :: a             !< the input scale factor
      type(TEFTCAMB_parameter_cache), intent(inout) :: eft_par_cache !< the EFTCAMB parameter cache that contains all the physical parameters.
      type(TEFTCAMB_timestep_cache ), intent(inout) :: eft_cache     !< the EFTCAMB timestep cache that contains all the physical values.

      real(dl) :: EFTCAMBRPHalphaDEComputeDtauda                           !< the output dtauda

      real(dl) :: temp

      temp = eft_cache%grhoa2 +eft_par_cache%grhov*a*a*self%RPHalphaDE_wDE%integral(a)
      EFTCAMBRPHalphaDEComputeDtauda = sqrt(3._dl/temp)

   end function EFTCAMBRPHalphaDEComputeDtauda

   ! ---------------------------------------------------------------------------------------------
   !> Subroutine that computes adotoa = H and its two derivatives wrt conformal time.
   subroutine EFTCAMBRPHalphaDEComputeAdotoa( self, a, eft_par_cache, eft_cache )

      implicit none

      class(EFTCAMB_RPH_alphaDE)                           :: self          !< the base class
      real(dl), intent(in)                         :: a             !< the input scale factor
      type(TEFTCAMB_parameter_cache), intent(inout) :: eft_par_cache !< the EFTCAMB parameter cache that contains all the physical parameters.
      type(TEFTCAMB_timestep_cache ), intent(inout) :: eft_cache     !< the EFTCAMB timestep cache that contains all the physical values.

      eft_cache%grhov_t = eft_par_cache%grhov*self%RPHalphaDE_wDE%integral(a)
      eft_cache%adotoa  = sqrt( ( eft_cache%grhom_t +eft_cache%grhov_t )/3._dl )

   end subroutine EFTCAMBRPHalphaDEComputeAdotoa

   ! ---------------------------------------------------------------------------------------------
   !> Subroutine that computes the two derivatives wrt conformal time of H.
   subroutine EFTCAMBRPHalphaDEComputeHubbleDer( self, a, eft_par_cache, eft_cache )

      implicit none

      class(EFTCAMB_RPH_alphaDE)                           :: self          !< the base class
      real(dl), intent(in)                         :: a             !< the input scale factor
      type(TEFTCAMB_parameter_cache), intent(inout) :: eft_par_cache !< the EFTCAMB parameter cache that contains all the physical parameters.
      type(TEFTCAMB_timestep_cache ), intent(inout) :: eft_cache     !< the EFTCAMB timestep cache that contains all the physical values.

      real(dl) :: wDE, dwDE, d2wDE, d3wDE, H, Ht, Ht2, grhor, grhom, grhov, grhon, gpn, gdpn, gddpn, gdddpn
      real(dl) :: grhonu, gpinu, gpinudot, gpinudotdot, gpinudotdotdot, grhormass_t
      integer :: nu_i

      wDE = self%RPHalphaDE_wDE%value(a)
      dwDE = a*self%RPHalphaDE_wDE%first_derivative(a)
      d2wDE = a*a*self%RPHalphaDE_wDE%second_derivative(a)
      H = eft_cache%adotoa
      grhom = eft_cache%grhob_t +eft_cache%grhoc_t
      grhor = eft_cache%grhor_t +eft_cache%grhog_t
      grhov = eft_cache%grhov_t

      grhon = 0._dl
      gpn = 0._dl
      gdpn = 0._dl
      gddpn = 0._dl
      ! gdddpn = 0._dl
      if ( eft_par_cache%Num_Nu_Massive /= 0 ) then
         grhon = eft_cache%grhonu_tot
         gpn = eft_cache%gpinu_tot
         gdpn = eft_cache%gpinudot_tot
         gddpn = eft_cache%gpinudotdot_tot
      end if

      eft_cache%gpiv_t  = wDE*grhov
      eft_cache%Hdot    = -0.5_dl*( H**2 +eft_cache%gpresm_t +eft_cache%gpiv_t )

      Ht = eft_cache%Hdot

      eft_cache%Hdotdot = H*( grhom/6._dl + 2._dl*grhor/3._dl ) &
      & +H*grhov*( 1._dl/6._dl +wDE +1.5_dl*wDE**2 -0.5_dl*dwDE ) &
      & +H*grhon/6._dl -0.5_dl*H*gpn -0.5_dl*gdpn

      ! needed by positivity in general:
      eft_cache%Hdotdotdot = H**2*( -grhom/6._dl - 4._dl*grhor/3._dl &
      & + grhov*(-1._dl/6._dl - 1.5_dl*wDE - 4.5_dl*wDE**2 - 4.5_dl*wDE**3 + dwDE + 4.5_dl*wDE*dwDE - 0.5_dl*d2wDE) &
      & -1._dl/6._dl*grhon - 1.5_dl*gpn ) &
      & + Ht*( grhom/6._dl + 2._dl/3._dl*grhor &
      & + grhov*(1._dl/6._dl + wDE +1.5_dl*wDE**2 - 0.5_dl*dwDE) &
      & + grhon/6._dl - 0.5_dl*gpn ) &
      & -1.5_dl*H*gdpn - 0.5_dl*gddpn

      ! only needed by this particular parameterization:
      ! Ht2 = eft_cache%Hdotdot
      ! if ( eft_par_cache%Num_Nu_Massive /= 0 ) then
      !    do nu_i = 1, eft_par_cache%Nu_mass_eigenstates
      !       grhonu      = 0._dl
      !       gpinu       = 0._dl
      !       gpinudot    = 0._dl
      !       gpinudotdot = 0._dl
      !       gpinudotdotdot = 0._dl
      !       grhormass_t = eft_par_cache%grhormass(nu_i)/a**2
      !       call ThermalNuBack%rho_P( a*eft_par_cache%nu_masses(nu_i), grhonu, gpinu )
      !       gpinudot = ThermalNuBack%pidot(a*eft_par_cache%nu_masses(nu_i),H,gpinu)
      !       gpinudotdot = ThermalNuBack%pidotdot(a*eft_par_cache%nu_masses(nu_i),H,Ht,gpinu,gpinudot)
      !       gpinudotdotdot = ThermalNuBack%pidotdotdot(a*eft_par_cache%nu_masses(nu_i),H,Ht,Ht2,gpinu,gpinudot,gpinudotdot)
      !       gdddpn = gdddpn + grhormass_t*(gpinudotdotdot &
      !       & - 12._dl*H*gpinudotdot &
      !       & + (48._dl*H**2 - 12._dl*Ht)*gpinudot &
      !       & + (-64._dl*H**3 + 48._dl*H*Ht - 4._dl*Ht2)*gpinu )
      !    end do
      ! end if
      ! d3wDE = a*a*a*self%RPHalphaDE_wDE%third_derivative(a)
      ! eft_cache%Hdotdotdotdot = H**3*( &
      ! & grhom/6._dl + 8._dl/3._dl*grhor &
      ! & + grhov*( 1._dl/6._dl + 2._dl*wDE + 9._dl*wDE**2 + 18._dl*wDE**3 + 13.5_dl*wDE**4 - 1.5_dl*dwDE - 12._dl*wDE*dwDE - 27._dl*wDE**2*dwDE + 4.5_dl*dwDE**2 + 0.5_dl*d2wDE + 6._dl*wDE*d2wDE - 0.5_dl*d3wDE ) &
      ! & + grhon/6._dl - 2.5_dl*gpn &
      ! & ) &
      ! & + H*Ht*( &
      ! & -0.5_dl*grhom - 4._dl*grhor &
      ! & + grhov*( -0.5_dl - 4.5_dl*wDE - 13.5_dl*wDE**2 - 13.5_dl*wDE**3 + 3._dl*dwDE +13.5_dl*wDE*dwDE - 1.5_dl*d2wDE ) &
      ! & -0.5_dl*grhon - 4.5_dl*gpn &
      ! & ) &
      ! & + Ht2*( &
      ! & grhom/6._dl + 2._dl/3._dl*grhor &
      ! & + grhov*( 1._dl/6._dl + wDE + 1.5_dl*wDE**2 - 0.5_dl*dwDE ) &
      ! & + grhon/6._dl - 0.5_dl*gpn &
      ! & ) &
      ! & + gdpn*( -2._dl*Ht - 4.5_dl*H**2 ) - 2.5_dl*H*gddpn - 0.5_dl*gdddpn

   end subroutine EFTCAMBRPHalphaDEComputeHubbleDer

   ! ---------------------------------------------------------------------------------------------
   !> Function that computes model specific stability requirements.
   function EFTCAMBRPHalphaDEAdditionalModelStability( self, a, eft_par_cache, eft_cache )

      implicit none

      class(EFTCAMB_RPH_alphaDE)                           :: self          !< the base class
      real(dl), intent(in)                         :: a             !< the input scale factor.
      type(TEFTCAMB_parameter_cache), intent(inout) :: eft_par_cache !< the EFTCAMB parameter cache that contains all the physical parameters.
      type(TEFTCAMB_timestep_cache ), intent(inout) :: eft_cache     !< the EFTCAMB timestep cache that contains all the physical values.

      logical :: EFTCAMBRPHalphaDEAdditionalModelStability          !< the return value of the stability computation. True if the model specific stability criteria are met, false otherwise.

      EFTCAMBRPHalphaDEAdditionalModelStability = .True.
      if ( self%RPHalphaDE_wDE%value(a) > -1._dl/3._dl ) EFTCAMBRPHalphaDEAdditionalModelStability = .False.

   end function EFTCAMBRPHalphaDEAdditionalModelStability

   ! ---------------------------------------------------------------------------------------------

end module EFTCAMB_Reparametrized_Horndeski_alphaDE
