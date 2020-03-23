!                                               30/09/2019
! SerraNA: ANALYSIS
! Calculates overall elastic constants
! and different definitions of persistence lengths
!
!    -----------------------------------------------------------------------
!    Copyright (C) 2019 Victor Velasco
!
!    SerraNA is free software: you can redistribute it and/or modify
!    it under the terms of the GNU Lesser General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    SerraNA is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with this program.  If not, see <https://www.gnu.org/licenses/>.
!    -----------------------------------------------------------------------

  program Analysis

  use io_mod
  use parms
  use functions_mod

  implicit none
  !nbp : number of base-pairs
  !frames: number of snapshots analysed
  !strands : type of structure (double or single stranded)
  !str : type of structure (closed or opened)
  !ranges(:,1) : contains the range in bp in which
  ! the overall ctes will be obtained.
  !ranges(:,2) : range of length in which ctes will be
  ! calculated
  integer :: nbp, frames, strands, str, n_bsp, ierror, &
           & ranges_A(2,2), ranges_T(2,2), &
           & ranges_S(2,2)

  !seq: sequence of strand I or II
  character(1), allocatable :: seq_I(:), seq_II(:)

  character(360) :: elasparms, strucparms

! Input parameters
!_____________________________________________________________________________________
  !elasp: elastic parameters
  !strucp: structural parameters
  !OV: Overall values
  real(dp), allocatable :: elasp(:,:), OV_elasp(:,:,:), &
                         & strucp(:,:,:), OV_strucp(:,:,:), &
                         & avstrp(:,:), OV_avstrp(:,:,:)

! Rearranged parameters
!_____________________________________________________________________________________
  !av:  averages
  !avstr : from average structure
  !b2 : Directional decay (or corrrelation)
  !E_E : Partial variance of end to end distance
  !av_Ad_c : Dynamic persitence length from dynamic tilt and roll
  real(dp), allocatable :: av_tilt(:), av_roll(:), av_b2(:), b2_avstr(:), &
                         & b2_d(:), av_twist(:), av_E_E(:), av_Ad_c(:)

! Elastic constants
!_____________________________________________________________________________________
 
  !Stretch(i) : =1 stretch, =2 intercept, =3 slope, =4 residual 
  !A : Persistence length. As : Static component. Ad : Dynamic component
  !_a(i) : same as strech
  !_bb : method b
  !_c : method c (dynamic through elastic ctes of tilt and roll)
  !_d : Persistence length using As_a and Ad_c
  real(dp) :: Tilt(2), Roll(2), Twist(2), Stretch(4), &
            & A_a(4), As_a(4), Ad_a(4), A_bb(2), Ad_c(2), &
            & A_d(2), &
            & r_aux !auxiliar variable

! Other...
!_____________________________________________________________________________________
  !These integers are auxiliars...
  integer :: A_n, T_n, S_n, i, l, aux_r, a, b, s, l1, l2

  real(dp), allocatable :: temp(:) !Auxiliar variable for linear fittings

  !Auxiliars helping to perform option for processing defaults
  logical :: whole_Af, whole_Al, whole_Tf, whole_Tl, whole_Sf, whole_Sl 

  !---------------------ANALYSIS------------------------------------------------------
  !-----------------------------------------------------------------------------------

  !INPUT SECTION
  !-----------------------------------------------------------------------------------

  !Get ranges and files
  call Analysis_inputs(ranges_A, ranges_T, ranges_S, &
                     & elasparms, strucparms)

  !READING SECTION
  !-----------------------------------------------------------------------------------

  call read_elastic_parms(elasparms, nbp, frames, strands, str, n_bsp, seq_I, seq_II, &
                        & elasp, OV_elasp)

  !Deallocate sequence to read again
  if (strands ==2) then
    deallocate(seq_I, seq_II, stat=ierror)
  else 
    deallocate(seq_I, stat=ierror)
  end if
  if (ierror/=0) stop "Error in deallocating sequence"
 
  call read_structural_parms(strucparms, nbp, frames, strands, str, n_bsp, seq_I, seq_II, &
                        & strucp, avstrp, OV_strucp, OV_avstrp)

  if (size(elasp,2) /= size(strucp,3)) &
  & stop "Are you sure you are giving me parameters for same system?"

  !CHECK INPUT DATA AND ALLOCATE AVERAGES SECTION
  !-----------------------------------------------------------------------------------

  !Check if elastic parameters have same length of structural parameters. If not
  !different structures could being analysed.
  if (size(strucp,3) /= size(elasp,2) ) then
    print*, "Are you sure elastic and structural parameters are from same structure?"
    print*, "Because I'm getting different sizes in each file."
    stop
  end if

  !Warnings (Check ranges)
  !-----------------------------------------------------------------------------------------------  

  !Persistence length's dependencies
  !---------------------------------
  a  = ranges_A(1,1)
  b  = ranges_A(2,1)
  l1 = ranges_A(1,2)
  l2 = ranges_A(2,2)

  ![a,b]
  if ( a == 0 .and. b == 0) then
    whole_Af = .true.
  else
    whole_Af = .false.
  end if

  if ( a == b .and. b /= 0 .and. a /= 0 ) stop "Invalid subfragment selection for Persistence length"
  if ( a < 1  .and. b /= 0 ) stop "Invalid subfragment selection for Persistence length"
  if ( b < 1  .and. a /= 0 ) stop "Invalid subfragment selection for Persistence length"
  if ( a > nbp ) stop "Invalid subfragment selection for Persistence length"
  if ( b > nbp ) stop "Invalid subfragment selection for Persistence length"

  if ( str == 1 ) then !This warning is only invalid for linear structures
    if ( b < a ) stop "Invalid subfragment selection for Persistence length"
  end if

  !(l1,l2)
  if ( l1 == 0 .and. l2 == 0) then
    whole_Al = .true.
  else
    whole_Al = .false.
  end if

  if ( l1 == l2 .and. l1 /= 0 .and. l2 /= 0 ) stop "Invalid sublength selection for Persistence length"
  if ( l1 < 1 .and. l2 /= 0 ) stop "Invalid sublength selection for Persistence length"
  if ( l2 < 1 .and. l1 /= 0 ) stop "Invalid sublength selection for Persistence length"
  if ( l2 < l1) stop "Invalid sublength selection for Persistence length"
  
  if ( a /= 0 .and. b /= 0) then

    if (b < a) then !s is an auxiliar that will help us to define the number of sublengths
      s = b + nbp   !This can only be true for closed structures
    else
      s = b
    end if

    if (l1 /= 0 .and. l2 /=0 ) then
      if ( l1 < 1 .or. l1 > abs( s - a ) ) stop "Invalid sublength selection for Persistence length"
      if ( l2 < 1 .or. l2 > abs( s - a ) ) stop "Invalid sublength selection for Persistence length"  
    end if

  end if

  !Twist
  !---------------------------------
  a  = ranges_T(1,1)
  b  = ranges_T(2,1)
  l1 = ranges_T(1,2)
  l2 = ranges_T(2,2)

  ![a,b]
  if ( a == 0 .and. b == 0) then
    whole_Tf = .true.
  else
    whole_Tf = .false.
  end if

  if ( a == b .and. b /= 0 .and. a /= 0 ) stop "Invalid subfragment selection for Twist"
  if ( a < 1  .and. b /= 0 ) stop "Invalid subfragment selection for Twist"
  if ( b < 1  .and. a /= 0 ) stop "Invalid subfragment selection for Twist"
  if ( a > nbp ) stop "Invalid subfragment selection for Twist"
  if ( b > nbp ) stop "Invalid subfragment selection for Twist"

  if ( str == 1 ) then !This warning is only invalid for linear structures
    if ( b < a ) stop "Invalid subfragment selection for Twist"
  end if

  !(l1,l2)
  if ( l1 == 0 .and. l2 == 0) then
    whole_Tl = .true.
  else
    whole_Tl = .false.
  end if

  if ( l1 == l2 .and. l1 /= 0 .and. l2 /= 0 ) stop "Invalid sublength selection for Twist"
  if ( l1 < 1 .and. l2 /= 0 ) stop "Invalid sublength selection for Twist"
  if ( l2 < 1 .and. l1 /= 0 ) stop "Invalid sublength selection for Twist"
  if ( l2 < l1) stop "Invalid sublength selection for Twist"
  
  if ( a /= 0 .and. b /= 0) then

    if (b < a) then !s is an auxiliar that will help us to define the number of sublengths
      s = b + nbp   !This can only be true for closed structures
    else
      s = b
    end if

    if (l1 /= 0 .and. l2 /=0 ) then
      if ( l1 < 1 .or. l1 > abs( s - a ) ) stop "Invalid sublength selection for Twist"
      if ( l2 < 1 .or. l2 > abs( s - a ) ) stop "Invalid sublength selection for Twist"  
    end if

  end if

  !Stretch
  !---------------------------------
  a  = ranges_S(1,1)
  b  = ranges_S(2,1)
  l1 = ranges_S(1,2)
  l2 = ranges_S(2,2)

  ![a,b]
  if ( a == 0 .and. b == 0) then

    !Check if it is circular and fix a,b
    if ( str == 2 ) then
 
      whole_Sf = .true. !If circular, then the whole fragment is considered
 
    else  !If linear

      if (nbp >= 18) then !Only if it is larger than a 17mer

        s = nbp/2 !Mid point, who cares if it is slightly on the left (by one bp)

        ranges_S(1,1) = s - 8 !a
        ranges_S(2,1) = s + 9 !b
        a  = ranges_S(1,1)
        b  = ranges_S(2,1)

        whole_Sf = .false.

      else

        whole_Sf = .true. !If short DNA

      end if !nbp >=18

    end if !str==2

  else

    whole_Sf = .false.

  end if  ![a,b]

  if ( a == b .and. b /= 0 .and. a /= 0 ) stop "Invalid subfragment selection for Stretch"
  if ( a < 1  .and. b /= 0 ) stop "Invalid subfragment selection for Stretch"
  if ( b < 1  .and. a /= 0 ) stop "Invalid subfragment selection for Stretch"
  if ( a > nbp ) stop "Invalid subfragment selection for Stretch"
  if ( b > nbp ) stop "Invalid subfragment selection for Stretch"

  if ( str == 1 ) then !This warning is only invalid for linear structures
    if ( b < a ) stop "Invalid subfragment selection for Stretch"
  end if

  !(l1,l2)
  if ( l1 == 0 .and. l2 == 0) then
    whole_Sl = .true.
  else
    whole_Sl = .false.
  end if

  if ( l1 == l2 .and. l1 /= 0 .and. l2 /= 0 ) stop "Invalid sublength selection for Stretch"
  if ( l1 < 1 .and. l2 /= 0 ) stop "Invalid sublength selection for Stretch"
  if ( l2 < 1 .and. l1 /= 0 ) stop "Invalid sublength selection for Stretch"
  if ( l2 < l1) stop "Invalid sublength selection for Stretch"
  
  if ( a /= 0 .and. b /= 0) then

    if (b < a) then !s is an auxiliar that will help us to define the number of sublengths
      s = b + nbp   !This can only be true for closed structures
    else
      s = b
    end if

    if (l1 /= 0 .and. l2 /=0 ) then
      if ( l1 < 1 .or. l1 > abs( s - a ) ) stop "Invalid sublength selection for Stretch"
      if ( l2 < 1 .or. l2 > abs( s - a ) ) stop "Invalid sublength selection for Stretch"  
    end if

  end if

  !Allocate Data
  !--------------------------------------  
  ! The second if's correspond to closed structures

  !Persistence length's dependencies
  a  = ranges_A(1,1)
  b  = ranges_A(2,1)

  if (a < b ) then 
    A_n = b - a
  else if (b > a) then
    A_n = nbp + b - a 
  else            !This is only possible if a == b == 0 which means that the whole fragmetn will be
    A_n = nbp - 1 !considered
  end if

  allocate(av_tilt(A_n), av_roll(A_n), av_Ad_c(A_n), av_b2(A_n), b2_avstr(A_n), b2_d(A_n), stat=ierror)
  if ( ierror /= 0 ) stop "Error in allocating persistences lengths parameters"

  !For twist
  a  = ranges_T(1,1)
  b  = ranges_T(2,1)

  if (a < b ) then 
    T_n = b - a
  else if (b > a) then
    T_n = nbp + b - a 
  else            !This is only possible if a == b == 0 which means that the whole fragmetn will be
    T_n = nbp - 1 !considered
  end if

  allocate(av_twist(T_n), stat=ierror)
  if ( ierror /= 0 ) stop "Error in allocating twist parameters"
 
  !For stretch
  a  = ranges_S(1,1)
  b  = ranges_S(2,1)

  if (a < b ) then 
    S_n = b - a
  else if (b > a) then
    S_n = nbp + b - a 
  else            !This is only possible if a == b == 0 
    S_n = nbp - 1 !If circular, then the whole fragment is considered
  end if

  allocate(av_E_E(S_n), stat=ierror)
  if ( ierror /= 0 ) stop "Error in allocating stretch parameters"

  !REARRENGING DATA SECTION
  !-----------------------------------------------------------------------------------
  ! centra_fragment subroutine obtaines the parameters in the desired range. This
  ! subroutine only outputs a vector (dimension one)
  ! 

  !TILT (4th element)
  if ( whole_Af ) then
    av_tilt = OV_elasp(1,4,:) !No need to calculate anything 
  else
    call central_fragment( av_tilt(:), elasp(4,:), ranges_A(1,1), & 
                                 & ranges_A(2,1), nbp, str )
  end if


  !ROLL (3rd element)
  if ( whole_Af ) then
    av_roll = OV_elasp(1,3,:) !No need to calculate anything 
  else
    call central_fragment( av_roll(:), elasp(3,:), ranges_A(1,1), & 
                               & ranges_A(2,1), nbp, str )
  end if


  !BEND (11th element) [It is actually the dynamic persistence length from tilt and roll]
  if ( whole_Af ) then
    av_Ad_c = OV_elasp(1,11,:) !No need to calculate anything 
  else
    call central_fragment( av_Ad_c(:), elasp(11,:), ranges_A(1,1), & 
                               & ranges_A(2,1), nbp, str )
  end if


  !TWIST (2nd element)
  if ( whole_Tf ) then
    av_twist = OV_elasp(1,2,:) !No need to calculate anything 
  else
    call central_fragment( av_twist(:), elasp(2,:), ranges_T(1,1), & 
                               & ranges_T(2,1), nbp, str )
  end if


  !PARTIAL VARIANCE OF END TO END DISTANCE (12th element)
  if ( whole_Sf ) then
    av_E_E = OV_elasp(1,13,:) !No need to calculate anything 
  else
    call central_fragment( av_E_E(:), elasp(13,:), ranges_S(1,1), & 
                               & ranges_S(2,1), nbp, str )
  end if


  !SQUARED BENDING ANGLE (10th element)
  !But we want 1-<a**2>/2
  !Squared bendings are in degree**2. We need to change to radian**2
  if ( whole_Af ) then
    av_b2 = 1.0_dp-0.5_dp*OV_strucp(1,10,:)*deg_to_rad*deg_to_rad !No need to calculate anything 
  else
    call central_fragment( av_b2(:), 1.0_dp-0.5_dp*strucp(1,10,:)*deg_to_rad*deg_to_rad, &
                       & ranges_A(1,1), ranges_A(2,1), nbp, str )
  end if


  !AVERAGE STRUCTURE BEND**2 (2nd element)
  if ( whole_Af ) then
    b2_avstr = 1.0_dp-0.5_dp*OV_avstrp(1,2,:)*deg_to_rad*deg_to_rad !No need to calculate anything 
  else
    call central_fragment( b2_avstr(:), 1.0_dp-0.5_dp*avstrp(2,:)*deg_to_rad*deg_to_rad, & 
                       & ranges_A(1,1), ranges_A(2,1), nbp, str )
  end if


  !==================================================================================
  !==================================================================================
  !==================================================================================
  !CALCULATE ELASTIC CONSTANTS SECTION
  !==================================================================================
  !==================================================================================
  !==================================================================================

  !-----------------------------------------------------------------------------------
  !LET US DEFINE THE RANGES FIRTS

  !PERSISTENCE LENGTH CRITERIA
  !-------------------------------
 
  if (whole_Al) then

    if (A_n >= 21) then  !We'll consider bulk behaviour and discard 10 longest sub-lengths
      ranges_A(1,2) = 11
      ranges_A(2,2) = A_n-10
    else
      ranges_A(1,2) = 1  !Or we'll consider everything
      ranges_A(2,2) = A_n
    end if

  end if

  !TWIST PERSISTENCE LENGTH CRITERIA
  !-------------------------------

  if (whole_Tl) then

    if (T_n >= 21) then  !We'll consider bulk behaviour and discard 10 longest sub-lengths
      ranges_T(1,2) = 11
      ranges_T(2,2) = A_n-10
    else
      ranges_T(1,2) = 1  !Or we'll consider everything
      ranges_T(2,2) = A_n
    end if

  end if

  !STRETCH MODULUS CRITERIA
  !-------------------------------

  if (whole_Sl) then

    if (S_n >= 17) then                     !Fit from 8 to 17
      ranges_S(1,2) = 8
      ranges_S(2,2) = 17
    else if (S_n >= 10 .and. S_n < 17) then !Fit from 8 to S_n
      ranges_S(1,2) = 8  
      ranges_S(2,2) = S_n
    else                                    !Fit of whatever (will yield something weird)
      ranges_S(1,2) = 1
      ranges_S(2,2) = S_n
    end if

  end if

  !-----------------------------------------------------------------------------------
  !NOW DO THE ACTUAL CALCULATIONS
  !-----------------------------------------------------------------------------------
  !-----------------------------------------------------------------------------------
  !TILT
  !-----------------------------------------------------------------------------------
  !let's recycle l. + 1 because we are counting the first element

  l = ranges_A(2,2) - ranges_A(1,2) + 1

  Tilt = average_std( av_tilt(ranges_A(1,2):ranges_A(2,2)), l )

  Tilt(2) = Tilt(2)/sqrt(real(l,dp))   !Change standard deviation to standard error
  !-----------------------------------------------------------------------------------
  !ROLL
  !-----------------------------------------------------------------------------------

  Roll = average_std( av_roll(ranges_A(1,2):ranges_A(2,2)), l)

  Roll(2) = Roll(2)/sqrt(real(l,dp))

  !-----------------------------------------------------------------------------------
  !BEND (Dynamic persistence length
  !-----------------------------------------------------------------------------------

  Ad_c = average_std( av_Ad_c(ranges_A(1,2):ranges_A(2,2)), l)

  Ad_c(2) = Ad_c(2)/sqrt(real(l,dp))
 
  !-----------------------------------------------------------------------------------
  !TWIST
  !-----------------------------------------------------------------------------------

  l = ranges_T(2,2) - ranges_T(1,2) + 1

  Twist = average_std( av_twist(ranges_T(1,2):ranges_T(2,2)), l)

  Twist(2) = Twist(2)/sqrt(real(l,dp))
  !-----------------------------------------------------------------------------------
  !STRETCH
  !-----------------------------------------------------------------------------------

  !Allocate temporal variable
  allocate(temp(ranges_S(2,2)-ranges_S(1,2)+1), stat=ierror)
  if (ierror /= 0 ) stop "Error in allocating auxiliar array in Stretch"

  !Set values of temporal variable
  l = 0
  do i = ranges_S(1,2),ranges_S(2,2)
    l = l + 1
    temp(l) = i
  end do

  !Do fit
  call simple_linear_regression( av_E_E(ranges_S(1,2):ranges_S(2,2)), temp(:), &
                               & Stretch(2), Stretch(3), Stretch(4) )  

  Stretch(1) = bnm*1.0E23_dp*bKTpN/Stretch(3)
  r_aux = bnm*1.0E23_dp*bKTpN*Stretch(4)
  Stretch(4) = abs( r_aux/Stretch(3)**2 ) !Confidence interval transformation:
                                          !f(x+$x) = f(x) + df(x)/dx *$x + ...
                                          !$f = f(x+$x) - f(x) ~ df(x)/dx *$x
                                          !Remember that our real stretch is a
                                          !non-linear function of the slope
  !clear temp
  deallocate(temp, stat=ierror)
  if (ierror /= 0 ) stop "Error in deallocating auxiliar array in Stretch"
  !-----------------------------------------------------------------------------------
  !PERSISTENCE LENGTHS
  !-----------------------------------------------------------------------------------
  !The persistence length and its static and dynamic contributions are calculated, by
  !four different methods:
  !(a) -> Linear Fit
  !(b) -> check
  !(c) -> Dynamic through roll and tilt
  !(d) -> Using static (a) and dynamic (c)

  !Allocate temporal variable
  aux_r = ranges_A(1,2) !So we can recover it
  ranges_A(1,2) = 1 !This is to perform the linear fit from 1 to r_A
  allocate(temp(ranges_A(2,2)-ranges_A(1,2)+1), stat=ierror)
  if (ierror /= 0 ) stop "Error in allocating auxiliar array in Persistence Lengths"

  !Set values of temporal variable
  l = 0
  do i = ranges_A(1,2),ranges_A(2,2)
    l = l + 1
    temp(l) = i
  end do
 
  !Persistence Length (a)
  !-----------------------------------------------------------------------------------
  
  !Linear fit with intercept 1 (we force it)
  call simple_linear_regression_a_1( av_b2(ranges_A(1,2):ranges_A(2,2)), temp(:), & 
                                   & A_a(2), A_a(3), A_a(4) )

  A_a(1) = -bnm/A_a(3)
  r_aux = bnm*A_a(4)
  A_a(4) = abs( r_aux/A_a(3)**2 ) !Confidence interval transformation
  
  !Static Persistence Length (a)
  !-----------------------------------------------------------------------------------

  call simple_linear_regression_a_1( b2_avstr(ranges_A(1,2):ranges_A(2,2)), temp(:), & 
                                   & As_a(2), As_a(3), As_a(4) )

  As_a(1) = -bnm/As_a(3)
  r_aux = bnm*As_a(4)
  As_a(4) = abs( r_aux/As_a(3)**2 )   !transformation
 
  !Dynamic Persistence Length (a)
  !-----------------------------------------------------------------------------------
  
  !Get difference of logarithms (Dynamic component)
  do i = 1, A_n
    b2_d(i) = 1.0_dp + av_b2(i) - b2_avstr(i)
  end do

  call simple_linear_regression_a_1( b2_d(ranges_A(1,2):ranges_A(2,2)), temp(:), & 
                                   & Ad_a(2), Ad_a(3), Ad_a(4) )

  Ad_a(1) = -bnm/Ad_a(3)
  r_aux = bnm*Ad_a(4)
  Ad_a(4) =abs( r_aux/Ad_a(3)**2 ) !Confidence interval transformation
 
  !Persistence Length (b)
  !-----------------------------------------------------------------------------------

  A_bb(1) = Ad_a(1)*As_a(1)/( Ad_a(1) + As_a(1) )
  
  !Persistence Length (d)
  !-----------------------------------------------------------------------------------

  A_d(1) = Ad_c(1)*As_a(1)/( Ad_c(1) + As_a(1) )


  !PRINT ELASTIC CONSTANTS SECTION
  !-----------------------------------------------------------------------------------
  
  ranges_A(1,2) = aux_r !recovered complete!

  call print_elastic_constants(Tilt, Roll, Twist, Stretch, A_a, As_a, Ad_a, A_bb, &
                             & Ad_c, A_d, ranges_A, ranges_S, ranges_T, &
                             & str, frames, nbp ) 

  end program Analysis
  !-----------------------------------------------------------------------------------
  !-----------------------------------------------------------------------------------

