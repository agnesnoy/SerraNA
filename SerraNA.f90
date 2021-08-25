!						25/08/2021
!
!
! SerraNA: MAIN
! Calculates BPP, BSP, structural parameters and elastic parameters
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
!
  program SerraNA

  use io_mod
  use parms
  use functions_mod

  implicit none

! R__1 = 1/n_R, Y__1 = 1/n_R, R__1_1 = 1/(n_R-1), Y__1_1 = 1/(n_Y-1)
! This will help to perform faster subroutine Get_RotationR_OriginO
  real(dp), parameter :: R__1 = 1.0_dp/real(n_R,dp), Y__1 = 1.0_dp/real(n_Y,dp), &
                       & R__1_1 = 0.125_dp, Y__1_1 = 0.2_dp
 
! nbp -> number of bp
! nb -> number of bases. If single stranded nb = nbp, if double nb = 2nbp
! box -> 1 if there's a box in trajectory file
! strands -> 1 if is a single stranded structure, 2 if double
! str -> 1 if is a linear structure, 2 if circular
! frames -> number of frames in trajectory file
! n_atoms -> total number of atoms
! ierror -> used in allocating and deallocating
! i_case -> for organizing by cases (see BPP and BSP section)
  integer :: nbp, nb, box, strands, str, frames, n_atoms, n_bsp, &
           & i, j, k, l, w, &
           & ierror, n, i_case !auxilar variables

! ring_index(j) -> array with locations of ring atom "j", in traj file
  integer, allocatable :: ring_index(:)

! seq(j) -> sequence of nucleotide j, the array has dimension str*nbp
  character(1), allocatable :: seq(:)

! top -> topology file
! traj -> trajectory file
  character(360) :: top, traj

! n__1, n__1_1 are auxiliar variables
! t_twist -> total twist: sum of twists from bp i to bp j
  real(dp) :: n__1, n__1_1, t_twist

! coord(:,j,k) -> j ring atom coordinates in time k
! R(:,i,j) -> axis i of base j
! O(:,j) -> position vector of base j
! Rmbt(:,i,j) -> axis i of basepair j (mid-base-triad)
! Ombt(:,j) -> position of basepair j
  real(dp), allocatable :: coords(:,:,:), R(:,:,:), O(:,:), & 
                           & Rmbt(:,:,:), Ombt(:,:)

! Structural parameters
!_____________________________________________________________________________________
! BPP(s,j,k) -> Base-pair parameters of basepair j in frame k
! s=1, Shear (Sx)           
! s=2, Streatch (Sy)       
! s=3, Stagger (Sz)       
! s=4, Buckle (kappa)      
! s=5, Propeller (omega)  
! s=6, Opening (sigma) 

! BSP(s,l,k) -> Base-step parameters between bp i and bp i+j, being j length in bp
! s=1, Shift (Dx)
! s=2, Slide (Dy)
! s=3, Rise (Dz)
! s=4, Tilt (tao)
! s=5, Roll (rho)
! s=6, Twist (Omega [Upper case])
! s=7, Bending (theta)
! s=8, Directional correlation  (dot product between z1 and z2 "cosb", b=bending)
! s=9, Bending squared (theta**2)

! E_E_dist(l,k) ->  End to end distance between bp i and bp i+j
! C_length(l,k) ->  Contour length between bp i and bp i+j
! added_bsp(s,l,k) -> Component s of contour length
  real(dp), allocatable :: BPP(:,:,:), BSP(:,:,:), E_E_dist(:,:), C_length(:,:), &
                         & added_bsp(:,:,:)

! Elastic parameters
!_____________________________________________________________________________________
 
! V(:,:,l) ->    Covariance matrix between bp i and bp i+j
! V__1(:,:,l) -> Invearse of covariance matrix
! F(s,l) ->      Coefficient s of Elastic matrix
! Ad2(l) -> Dynamic persistence length through Tilt and Roll from F
  real(dp), allocatable :: V(:,:,:), V__1(:,:,:), F(:,:), Ad2(:)

! Average and standard deviation parameters
!_____________________________________________________________________________________
! avstd_ parameters are averages and standard deviations of BPP, BSP, end to end
! distance, contour length and added parameters. 1 of first dimension is the average
! and 2 is the standard deviation.
! BSP_avstr are BSP parameters of the average structure
  real(dp), allocatable :: avstd_BPP(:,:,:), avstd_BSP(:,:,:), avstd_E_E(:,:), &
                         & avstd_C_l(:,:), avstd_added(:,:,:), BSP_avstr(:,:)

! Overall parameters
!_____________________________________________________________________________________
! These parameters are averages of all parameters calculated by length.
! First dimension corresponds to average and standard deviation
! For arrays of three dimensions, the second one represents a parameter.
! The last one is the bp distance (1bp, 2bp, ..., nbp-1 bp).
! V_E_E is the variance of the end to end distance and pV_E_E is the partial variance.
! OV_BSP_avstr are the overall BSP of the average structure
! OV_Ad2 overall of dynamic persistence length
  real(dp), allocatable :: OV_BPP(:,:), OV_BSP(:,:,:), OV_E_E(:,:), OV_C_l(:,:), &
           & OV_added(:,:,:), OV_V_E_E(:,:), OV_pV_E_E(:,:), OV_F(:,:,:), OV_Ad2(:,:), & 
           & OV_BSP_avstr(:,:,:)
 

  !-----------------------------MAIN--------------------------------------------------
  !-----------------------------------------------------------------------------------

  !INPUT SECTION
  !-----------------------------------------------------------------------------------
  !Get the directories of trajectory and topology files. Also, read what type of 
  !structure is being analysed, single or double stranded, linear or circular.

  call SerraNA_inputs(traj,top,strands,str)

  !READING  SECTION
  !----------------------------------------------------------------------------------- 
  write(6,*) "Reading topology file"
  call topology_amber(nbp,n_atoms,box,seq,ring_index,top,strands)
  
  write(6,*) "Reading trajectory file"
  call coordinates_amber_crd(coords, frames, ring_index, n_atoms, seq, nbp, traj, & 
                           & box,strands, str )

  deallocate(ring_index, stat=ierror)
  if (ierror/=0) stop "Error in deallocating ring_index"

  !STRUCTURAL PARAMETERS  SECTION
  !----------------------------------------------------------------------------------- 
  write(6,*) "Calculating structural parameters"

  !Number of bases
  nb = strands*nbp
  !Number of basestep parameters  
  if (str == 2) then
    n_bsp = nbp*(nbp-1)         !If it's circular, then we have more parameters
  else                          !nbp-1 for each bp
    n_bsp = nbp*(nbp-1)/2
  end if

  !Let's organise the cases
  i_case = 0
  if ( (str == 1) .and. (strands == 2) ) then
    i_case = 1                                     ! Linear double stranded
  else if ( (str == 1) .and. (strands == 1) ) then
    i_case = 2                                     ! Linear single stranded
  else if ( (str == 2) .and. (strands == 2) ) then
    i_case = 3                                     ! Circular double stranded
  else if ( (str == 2) .and. (strands == 1) ) then
    i_case = 4                                     ! Circular single stranded
  else 
    stop "Error in determining case"
  end if

  !Allocate orientations matrices and origin vectors
  allocate(R(3,3,nb), O(3,nb), stat=ierror)
  if (ierror/=0) stop "Error in allocating orientations and origins"

  if (strands == 2) then  !MBT is not needed if structure is single stranded
    allocate(Rmbt(3,3,nbp), Ombt(3,nbp), stat=ierror)
    if (ierror/=0) stop "Error in allocating MBT"
  end if

  !Allocate basepair parameters and basestep parameters
  allocate(BSP(9,n_bsp,frames), stat=ierror)
  if (ierror/=0) stop "Error in allocating basestep parameters"

  if (strands ==2) then      !Not needed if single stranded structure
    allocate(BPP(6,nbp,frames), stat=ierror)
    if (ierror/=0) stop "Error in allocating basepair parameters"
  end if

  !Allocate end to end distance and contour length, and added parameters (
  allocate(E_E_dist(n_bsp,frames), C_length(n_bsp,frames), added_bsp(3,n_bsp,frames), stat=ierror)
  if (ierror/=0) stop "Error in allocating end to end distance, contour length and added parameters"


  !Time to select case
  !The triad method for calculating bp positions and orientations are 
  !the same in all cases. We have a bigger code but better 
  !organised (hopefully)
  select case( i_case )

  ! CASE 1! LINEAR DOUBLE STRANDED STRUCTURE
  !----------------------------------------------------------------------
  case(1) 

  !k loop will cover all the calculations 
  do k=1,frames
    l=1
    do i=1,nb
      !BASE origins O and orientitions R
      !Check how many atoms will be taken per base i.
      ! n_R for purines and n_Y for pyrimidines, these are stored in parms.f90
      if ( (seq(i) .eq. "G") .or. (seq(i) .eq. "A") ) then
        n = n_R !purines
        n__1 = R__1
        n__1_1 = R__1_1
      else
        n = n_Y !pyrimidines
        n__1 = Y__1
        n__1_1 = Y__1_1
      end if  
      !calculate orientation matrices "R" and origins "O" for each base
      ! G_b, A_b, C_b, T_b are parameter bases stored in parms.f90
      if (seq(i) .eq. 'G') then
        call Get_Rotation_R_Origin_O(R(:,:,i),O(:,i),coords(:,l:l+n-1,k),G_b,n,n__1,n__1_1)
      else if(seq(i) .eq. 'A') then
        call Get_Rotation_R_Origin_O(R(:,:,i),O(:,i),coords(:,l:l+n-1,k),A_b,n,n__1,n__1_1)
      else if(seq(i) .eq. 'C') then
        call Get_Rotation_R_Origin_O(R(:,:,i),O(:,i),coords(:,l:l+n-1,k),C_b,n,n__1,n__1_1)
      else if(seq(i) .eq. 'T') then
        call Get_Rotation_R_Origin_O(R(:,:,i),O(:,i),coords(:,l:l+n-1,k),T_b,n,n__1,n__1_1)
      else if(seq(i) .eq. 'U') then !Uracil
        call Get_Rotation_R_Origin_O(R(:,:,i),O(:,i),coords(:,l:l+n-1,k),U_b,n,n__1,n__1_1)
      end if
      l=l+n
    end do !close i

    !BASEPAIR AND BASESTEP PARAMETERS
    !----------------------------------------------------------------------
    !Basepair parameters only for double stranded structures 
    !Reverse the y- and z- axis from second strand.
    do i=nbp+1,nb
      do j=2,3
        R(:,j,i) = -R(:,j,i)
      end do
    end do
    !Get BPP
    do i=1,nbp
      call basepair_parameters(O(:,i), O(:,nb+1-i), R(:,:,i), R(:,:,nb+1-i), &
                             & BPP(:,i,k), Rmbt(:,:,i), Ombt(:,i) )
    end do 

    !Base-step Parameters, end to end distances, contour lengths and added parameters
    l = 0                                        !index
    
    !******************************************Double Stranded Structure
    !First BSP, lengths of 1 bp
    j=1 
    t_twist = 0.0_dp                          
    do i=1,nbp-j                               !first bp
      l = l + 1
      call basestep_parameters(Ombt(:,i), Ombt(:,i+j), Rmbt(:,:,i), Rmbt(:,:,i+j), &
                             & BSP(:,l,k), t_twist)
      E_E_dist(l,k) = absv( Ombt(:,i+j) - Ombt(:,i) )
      C_length(l,k) = E_E_dist(l,k)            !at lengths of 1bp, these are the same
      added_bsp(:,l,k) = BSP(1:3,i,k)          !i=l
    end do      !close i
    !Lengths > 1 bp
    do j=2,nbp-1                               !from length of two bp
      do i=1,nbp-j                             !first bp
        l = l + 1
        t_twist = BSP(6,l-nbp+j-1,k)        !Twist from BSP i to bp i+j-1
        call basestep_parameters(Ombt(:,i), Ombt(:,i+j), Rmbt(:,:,i), Rmbt(:,:,i+j), &
                               & BSP(:,l,k), t_twist)
        E_E_dist(l,k) = absv( Ombt(:,i+j) - Ombt(:,i) )
        C_length(l,k) = C_length(l-nbp+j-1,k) + C_length(i+j-1,k)
        added_bsp(:,l,k) = added_bsp(:,l-nbp+j-1,k) + BSP(1:3,i+j-1,k)
      end do
    end do  

  end do !close k


  ! CASE 2! LINEAR SINGLE STRANDED STRUCTURE
  !----------------------------------------------------------------------
  case(2)      

  !k loop will cover all the calculations 
  do k=1,frames
    l=1
    do i=1,nb
      !BASE origins O and orientitions R
      !Check how many atoms will be taken per base i.
      ! n_R for purines and n_Y for pyrimidines, these are stored in parms.f90
      if ( (seq(i) .eq. "G") .or. (seq(i) .eq. "A") ) then
        n = n_R !purines
        n__1 = R__1
        n__1_1 = R__1_1
      else
        n = n_Y !pyrimidines
        n__1 = Y__1
        n__1_1 = Y__1_1
      end if  
      !calculate orientation matrices "R" and origins "O" for each base
      ! G_b, A_b, C_b, T_b are parameter bases stored in parms.f90
      if (seq(i) .eq. 'G') then
        call Get_Rotation_R_Origin_O(R(:,:,i),O(:,i),coords(:,l:l+n-1,k),G_b,n,n__1,n__1_1)
      else if(seq(i) .eq. 'A') then
        call Get_Rotation_R_Origin_O(R(:,:,i),O(:,i),coords(:,l:l+n-1,k),A_b,n,n__1,n__1_1)
      else if(seq(i) .eq. 'C') then
        call Get_Rotation_R_Origin_O(R(:,:,i),O(:,i),coords(:,l:l+n-1,k),C_b,n,n__1,n__1_1)
      else if(seq(i) .eq. 'T') then
        call Get_Rotation_R_Origin_O(R(:,:,i),O(:,i),coords(:,l:l+n-1,k),T_b,n,n__1,n__1_1)
      else if(seq(i) .eq. 'U') then !Uracil
        call Get_Rotation_R_Origin_O(R(:,:,i),O(:,i),coords(:,l:l+n-1,k),U_b,n,n__1,n__1_1)
      end if
      l=l+n
    end do !close i

    ! BASESTEP PARAMETERS (NO BASEPAIR PARAMETERS FOR SINGLE STRANDED STRUCTURES)
    !----------------------------------------------------------------------

    !Base-step Parameters, end to end distances, contour lengths and added parameters
    l = 0                                        !index
    !Same process than double stranded, but using O and R instead of 
    !Ombt and Rmbt
    j=1 
    t_twist = 0.0_dp                          
    do i=1,nbp-j                               !first bp
      l = l + 1
      call basestep_parameters(O(:,i), O(:,i+j), R(:,:,i), R(:,:,i+j), &
                             & BSP(:,l,k), t_twist)
      E_E_dist(l,k) = absv( O(:,i+j) - O(:,i) )
      C_length(l,k) = E_E_dist(l,k)            !at lengths of 1bp, these are the same
      added_bsp(:,l,k) = BSP(1:3,i,k)          !i=l
    end do      !close i
    !Lengths > 1 bp
    do j=2,nbp-1                               !from length of two bp
      do i=1,nbp-j                             !first bp
        l = l + 1
        t_twist = BSP(6,l-nbp+j-1,k)        !Twist from BSP i to bp i+j-1
        call basestep_parameters(O(:,i), O(:,i+j), R(:,:,i), R(:,:,i+j), &
                                 BSP(:,l,k), t_twist)
        E_E_dist(l,k) = absv( O(:,i+j) - O(:,i) )
        C_length(l,k) = C_length(l-nbp+j-1,k) + C_length(i+j-1,k)
        added_bsp(:,l,k) = added_bsp(:,l-nbp+j-1,k) + BSP(1:3,i+j-1,k)
      end do
    end do  

  end do !close k


  ! CASE 3! CIRCULAR DOUBLE STRANDED STRUCTURE
  !----------------------------------------------------------------------
  case(3)   

  !k loop will cover all the calculations 
  do k=1,frames
    l=1
    do i=1,nb
      !BASE origins O and orientitions R
      !Check how many atoms will be taken per base i.
      ! n_R for purines and n_Y for pyrimidines, these are stored in parms.f90
      if ( (seq(i) .eq. "G") .or. (seq(i) .eq. "A") ) then
        n = n_R !purines
        n__1 = R__1
        n__1_1 = R__1_1
      else
        n = n_Y !pyrimidines
        n__1 = Y__1
        n__1_1 = Y__1_1
      end if  
      !calculate orientation matrices "R" and origins "O" for each base
      ! G_b, A_b, C_b, T_b are parameter bases stored in parms.f90
      if (seq(i) .eq. 'G') then
        call Get_Rotation_R_Origin_O(R(:,:,i),O(:,i),coords(:,l:l+n-1,k),G_b,n,n__1,n__1_1)
      else if(seq(i) .eq. 'A') then
        call Get_Rotation_R_Origin_O(R(:,:,i),O(:,i),coords(:,l:l+n-1,k),A_b,n,n__1,n__1_1)
      else if(seq(i) .eq. 'C') then
        call Get_Rotation_R_Origin_O(R(:,:,i),O(:,i),coords(:,l:l+n-1,k),C_b,n,n__1,n__1_1)
      else if(seq(i) .eq. 'T') then
        call Get_Rotation_R_Origin_O(R(:,:,i),O(:,i),coords(:,l:l+n-1,k),T_b,n,n__1,n__1_1)
      else if(seq(i) .eq. 'U') then !Uracil
        call Get_Rotation_R_Origin_O(R(:,:,i),O(:,i),coords(:,l:l+n-1,k),U_b,n,n__1,n__1_1)
      end if
      l=l+n
    end do !close i

    !BASEPAIR AND BASESTEP PARAMETERS
    !----------------------------------------------------------------------
    !Basepair parameters only for double stranded structures 
    !Reverse the y- and z- axis from second strand.
    do i=nbp+1,nb
      do j=2,3
        R(:,j,i) = -R(:,j,i)
      end do
    end do
    !Get BPP
    do i=1,nbp
      call basepair_parameters(O(:,i), O(:,nb+1-i), R(:,:,i), R(:,:,nb+1-i), &
                             & BPP(:,i,k), Rmbt(:,:,i), Ombt(:,i) )
    end do 

    !Base-step Parameters, end to end distances, contour lengths and added parameters
    l = 0                                        !index
    
    !******************************************Double Stranded Structure
    !First BSP, lengths of 1 bp
    j=1 
    t_twist = 0.0_dp                          
    do i=1,nbp-j                               !first bp
      l = l + 1
      w = i+j
      call basestep_parameters(Ombt(:,i), Ombt(:,w), Rmbt(:,:,i), Rmbt(:,:,w), &
                             & BSP(:,l,k), t_twist)
      E_E_dist(l,k) = absv( Ombt(:,w) - Ombt(:,i) )
      C_length(l,k) = E_E_dist(l,k)            !at lengths of 1bp, these are the same
      added_bsp(:,l,k) = BSP(1:3,i,k)          !i=l
    end do      !close i
    !Other side
    do i=nbp-j+1,nbp
      w = i+j-nbp
      l = l + 1
      call basestep_parameters(Ombt(:,i), Ombt(:,w), Rmbt(:,:,i), Rmbt(:,:,w), &
                             & BSP(:,l,k), t_twist)
      E_E_dist(l,k) = absv( Ombt(:,w) - Ombt(:,i) )
      C_length(l,k) = E_E_dist(l,k)            !at lengths of 1bp, these are the same
      added_bsp(:,l,k) = BSP(1:3,i,k)          !i=l
    end do
    !Lengths > 1 bp
    do j=2,nbp-1                               !from length of two bp
      do i=1,nbp-j                             !first bp
        w = i+j
        l = l + 1
        t_twist = BSP(6,l-nbp,k)        !Twist from BSP i to bp i+j-1
        call basestep_parameters(Ombt(:,i), Ombt(:,w), Rmbt(:,:,i), Rmbt(:,:,w), &
                               & BSP(:,l,k), t_twist)
        E_E_dist(l,k) = absv( Ombt(:,w) - Ombt(:,i) )
        C_length(l,k) = C_length(l-nbp,k) + C_length(w-1,k)
        added_bsp(:,l,k) = added_bsp(:,l-nbp,k) + BSP(1:3,w-1,k)
      end do
     !Other side
      i = nbp-j+1
      w = 1
      l = l + 1
        t_twist = BSP(6,l-nbp,k)        !Twist from BSP i to bp i+j-1
      call basestep_parameters(Ombt(:,i), Ombt(:,w), Rmbt(:,:,i), Rmbt(:,:,w), &
                             & BSP(:,l,k), t_twist)
      E_E_dist(l,k) = absv( Ombt(:,w) - Ombt(:,i) )
      C_length(l,k) = C_length(l-nbp,k) + C_length(nbp,k)
      added_bsp(:,l,k) = added_bsp(:,l-nbp,k) + BSP(1:3,nbp,k)
      do i=nbp-j+2,nbp     
        w = i+j-nbp
        l = l + 1
        t_twist = BSP(6,l-nbp,k)        !Twist from BSP i to bp i+j-1
        call basestep_parameters(Ombt(:,i), Ombt(:,w), Rmbt(:,:,i), Rmbt(:,:,w), &
                               & BSP(:,l,k), t_twist)
        E_E_dist(l,k) = absv( Ombt(:,w) - Ombt(:,i) )
        C_length(l,k) = C_length(l-nbp,k) + C_length(w-1,k)
        added_bsp(:,l,k) = added_bsp(:,l-nbp,k) + BSP(1:3,w-1,k)
      end do
    end do  

  end do !close k


  ! CASE 4! CIRCULAR SINGLE STRANDED STRUCTURE
  !----------------------------------------------------------------------
  case(4) 

  !k loop will cover all the calculations 
  do k=1,frames
    l=1
    do i=1,nb
      !BASE origins O and orientitions R
      !Check how many atoms will be taken per base i.
      ! n_R for purines and n_Y for pyrimidines, these are stored in parms.f90
      if ( (seq(i) .eq. "G") .or. (seq(i) .eq. "A") ) then
        n = n_R !purines
        n__1 = R__1
        n__1_1 = R__1_1
      else
        n = n_Y !pyrimidines
        n__1 = Y__1
        n__1_1 = Y__1_1
      end if  
      !calculate orientation matrices "R" and origins "O" for each base
      ! G_b, A_b, C_b, T_b are parameter bases stored in parms.f90
      if (seq(i) .eq. 'G') then
        call Get_Rotation_R_Origin_O(R(:,:,i),O(:,i),coords(:,l:l+n-1,k),G_b,n,n__1,n__1_1)
      else if(seq(i) .eq. 'A') then
        call Get_Rotation_R_Origin_O(R(:,:,i),O(:,i),coords(:,l:l+n-1,k),A_b,n,n__1,n__1_1)
      else if(seq(i) .eq. 'C') then
        call Get_Rotation_R_Origin_O(R(:,:,i),O(:,i),coords(:,l:l+n-1,k),C_b,n,n__1,n__1_1)
      else if(seq(i) .eq. 'T') then
        call Get_Rotation_R_Origin_O(R(:,:,i),O(:,i),coords(:,l:l+n-1,k),T_b,n,n__1,n__1_1)
      else if(seq(i) .eq. 'U') then !Uracil
        call Get_Rotation_R_Origin_O(R(:,:,i),O(:,i),coords(:,l:l+n-1,k),U_b,n,n__1,n__1_1)
      end if
      l=l+n
    end do !close i

    ! BASESTEP PARAMETERS (NO BASEPAIR PARAMETERS FOR SINGLE STRANDED STRUCTURES)
    !----------------------------------------------------------------------

    !Base-step Parameters, end to end distances, contour lengths and added parameters
    l = 0                                        !index
    !Same process than double stranded, but using O and R instead of 
    !Ombt and Rmbt
    j=1 
    t_twist = 0.0_dp                          
    do i=1,nbp-j                               !first bp
      l = l + 1
      w = i+j
      call basestep_parameters(O(:,i), O(:,w), R(:,:,i), R(:,:,w), &
                             & BSP(:,l,k), t_twist)
      E_E_dist(l,k) = absv( O(:,w) - O(:,i) )
      C_length(l,k) = E_E_dist(l,k)            !at lengths of 1bp, these are the same
      added_bsp(:,l,k) = BSP(1:3,i,k)          !i=l
    end do      !close i
    !Other side
    do i=nbp-j+1,nbp
     w = i+j-nbp
      l = l + 1
      call basestep_parameters(O(:,i), O(:,w), R(:,:,i), R(:,:,w), &
                             & BSP(:,l,k), t_twist)
      E_E_dist(l,k) = absv( O(:,w) - O(:,i) )
      C_length(l,k) = E_E_dist(l,k)            !at lengths of 1bp, these are the same
      added_bsp(:,l,k) = BSP(1:3,i,k)          !i=l
    end do
    
    !Lengths > 1 bp
    do j=2,nbp-1                               !from length of two bp
      do i=1,nbp-j                             !first bp
        l = l + 1
        w = i+j
        t_twist = BSP(6,l-nbp,k)        !Twist from BSP i to bp i+j-1
        call basestep_parameters(O(:,i), O(:,w), R(:,:,i), R(:,:,w), &
                                 BSP(:,l,k), t_twist)
        E_E_dist(l,k) = absv( O(:,w) - O(:,i) )
        C_length(l,k) = C_length(l-nbp,k) + C_length(w-1,k)
        added_bsp(:,l,k) = added_bsp(:,l-nbp,k) + BSP(1:3,w-1,k)
      end do
      !Other side
      i = nbp-j+1
      w = 1
      l = l + 1
      t_twist = BSP(6,l-nbp,k)        !Twist from BSP i to bp i+j-1
      call basestep_parameters(O(:,i), O(:,w), R(:,:,i), R(:,:,w), &
                             & BSP(:,l,k), t_twist)
      E_E_dist(l,k) = absv( O(:,w) - O(:,i) )
      C_length(l,k) = C_length(l-nbp,k) + C_length(nbp,k)
      added_bsp(:,l,k) = added_bsp(:,l-nbp,k) + BSP(1:3,nbp,k)
 
      do i=nbp-j+2,nbp     
        w = i+j-nbp
        l = l + 1
        t_twist = BSP(6,l-nbp,k)        !Twist from BSP i to bp i+j-1
        call basestep_parameters(O(:,i), O(:,w), R(:,:,i), R(:,:,w), &
                               & BSP(:,l,k), t_twist)
        E_E_dist(l,k) = absv( O(:,w) - O(:,i) )
        C_length(l,k) = C_length(l-nbp,k) + C_length(w-1,k)
        added_bsp(:,l,k) = added_bsp(:,l-nbp,k) + BSP(1:3,w-1,k)
      end do
    end do  

  end do !close k


  case default   
    stop "Error in select case, structural parameters section"
  end select
  !END SELECT CASES...
  !----------------------------------------------------------------------------------- 

  !Structural parameters calculated
  !Get rid of what no longer is needed
  deallocate(coords, R, O, stat=ierror)
  if (ierror/= 0) stop "Error in deallocating coordinates, R and O"

  if (strands ==2) then
    deallocate(Ombt,Rmbt, stat=ierror)
    if (ierror/= 0) stop "Error in deallocating MBT arrays"
  end if

  !ELASTIC PARAMETERS  SECTION
  !----------------------------------------------------------------------------------- 
  if (frames > 1) then
    write(6,*) "Calculating elastic parameters"

    allocate(V(4,4,n_bsp), V__1(4,4,n_bsp), F(10,n_bsp), Ad2(n_bsp), stat=ierror)
    if (ierror /= 0) stop "Error in allocating elastic parameters"

    !calculate covariance matrix, inverse of covariance matrix and elastic coeficients
    l=0 
    do j=1,nbp-1 

      if (str == 2) then  !If it is circular

        if ( j < nbp - j ) then
          w = j
        else
          w = nbp -j
        end if

      else     !Linear
        w = j      
      end if 

      do i=1,nbp+j*(str-2)  !if str ==2, then i =1, nbp. if str ==1 then i =1, nbp-j
        l=l+1
        call deformation_covariance( V(:,:,l), BSP(4,l,:), BSP(5,l,:), &
                                   & BSP(6,l,:), E_E_dist(l,:), frames )
        call inverse_matrix_analytic4x4( V(:,:,l), V__1(:,:,l) )

        F(:,l) = elastic_F( V__1(:,:,l), w)

        !And the dynamic persistence length through tilt and roll | 4 -> tilt, 3 -> roll
        Ad2(l) = dynamic_persistence_length2( F(4,l), F(3,l) )
      end do
    end do
  else
    write(6,*) "Only one frame read"
    write(6,*) "Cannot calculate elastic parameters"
  end if


  !AVERAGES AND STANDARD DEVIATION SECTION
  !----------------------------------------------------------------------------------- 
  write(6,*) "Calculating averages and standard deviations"

  !Allocate arrays
  if (strands==2) then
    allocate(avstd_BPP(2,6,nbp), avstd_BSP(2,9,n_bsp), avstd_E_E(2,n_bsp), &
           & avstd_C_l(2,n_bsp), avstd_added(2,3,n_bsp), stat=ierror)
  else
    allocate(avstd_BSP(2,9,n_bsp), avstd_E_E(2,n_bsp), &
           & avstd_C_l(2,n_bsp), avstd_added(2,3,n_bsp), stat=ierror)
  end if
  if (ierror/=0) stop "Error in allocating average parameters"

  !BPP if there are BPP
  if (strands==2) then
    do i=1,nbp
      do j=1,6
        avstd_BPP(:,j,i) = average_std( BPP(j,i,:), frames )
      end do
    end do
    deallocate(BPP, stat=ierror)
    if (ierror/=0) stop "Error in deallocating base-pair parameters"
  end if
  
  !BSP
  do i=1,n_bsp
    do j=1,9
      avstd_BSP(:,j,i) = average_std( BSP(j,i,:), frames )
    end do
  end do
  
  !End to end distance
  do i=1,n_bsp
    avstd_E_E(:,i) = average_std( E_E_dist(i,:), frames )
  end do

  !Contour length
  do i=1,n_bsp
    avstd_C_l(:,i) = average_std( C_length(i,:), frames )
  end do

  !Added parameters
  do i=1,n_bsp
    do j=1,3
      avstd_added(:,j,i) = average_std( added_bsp(j,i,:), frames )
    end do
  end do

  deallocate(BSP, E_E_dist, C_length, added_bsp, stat=ierror)
  if (ierror/=0) stop "Error in deallocating structural parameters"
  
  !AVERAGE STRUCTURE SECTION
  !----------------------------------------------------------------------------------- 
  write(6,*) "Calculating average structure parameters"

  !First, we need to obtain the triads for the average structure,
  !these will be obtained through the average values of BSP [avstd_BSP(1,:)]

  !let's re-use R and O and get the base triads of the average structure
  allocate(Rmbt(3,3,nbp), Ombt(3,nbp), BSP_avstr(9,n_bsp), stat=ierror)
  if (ierror/=0) stop "Error in allocating average structure parameters"

  call reverse_algorithm_base_triads(Rmbt,Ombt,avstd_BSP(1,:,1:nbp),nbp) !get R and O

  if (str /= 2) then
  !LINEAR STRUCTURE
 
    !First BSP, lengths of 1 bp
    l = 0
    j = 1 
    t_twist = 0.0_dp                           !First twists are 0
    do i=1,nbp-j                               !first bp
      l = l + 1
      call basestep_parameters(Ombt(:,i), Ombt(:,i+j), Rmbt(:,:,i), Rmbt(:,:,i+j), &
                            &  BSP_avstr(:,l), t_twist)
    end do      !close i
    !Lengths > 1 bp
    do j=2,nbp-1                               !from length of two bp
      do i=1,nbp-j                             !first bp
        l = l + 1
        t_twist = BSP_avstr(6,l-nbp+j-1)        !Twist from BSP i to bp i+j-1
        call basestep_parameters(Ombt(:,i), Ombt(:,i+j), Rmbt(:,:,i), Rmbt(:,:,i+j), &
                              &  BSP_avstr(:,l), t_twist)
      end do
    end do  

  else
  !CLOSED STRUCTURE
    !First BSP, lengths of 1 bp
    j=1
    l = 0 
    t_twist = 0.0_dp                          
    do i=1,nbp-j                               !first bp
      l = l + 1
      w = i+j
      call basestep_parameters(Ombt(:,i), Ombt(:,w), Rmbt(:,:,i), Rmbt(:,:,w), &
                             & BSP_avstr(:,l), t_twist)
    end do      !close i
    !Other side
    do i=nbp-j+1,nbp
      w = i+j-nbp
      l = l + 1
      call basestep_parameters(Ombt(:,i), Ombt(:,w), Rmbt(:,:,i), Rmbt(:,:,w), &
                             & BSP_avstr(:,l), t_twist)
    end do
    !Lengths > 1 bp
    do j=2,nbp-1                               !from length of two bp
      do i=1,nbp-j                             !first bp
        w = i+j
        l = l + 1
        t_twist = BSP_avstr(6,l-nbp)        !Twist from BSP i to bp i+j-1
        call basestep_parameters(Ombt(:,i), Ombt(:,w), Rmbt(:,:,i), Rmbt(:,:,w), &
                               & BSP_avstr(:,l), t_twist)
      end do
      !Other side
      i = nbp-j+1
      w = 1
      l = l + 1
      t_twist = BSP_avstr(6,l-nbp)        !Twist from BSP i to bp i+j-1
      call basestep_parameters(Ombt(:,i), Ombt(:,w), Rmbt(:,:,i), Rmbt(:,:,w), &
                             & BSP_avstr(:,l), t_twist)
 
      do i=nbp-j+2,nbp     
        w = i+j-nbp
        l = l + 1
        t_twist = BSP_avstr(6,l-nbp)        !Twist from BSP i to bp i+j-1
        call basestep_parameters(Ombt(:,i), Ombt(:,w), Rmbt(:,:,i), Rmbt(:,:,w), &
                               & BSP_avstr(:,l), t_twist)
      end do
    end do  

  end if
 
  deallocate(Rmbt, Ombt, stat=ierror) !No longer needed

  !OVERAL PARAMETERS SECTION
  !----------------------------------------------------------------------------------- 
  write(6,*) "Calculating overall parameters"
  
  !Allocate arrays
  if (strands==2) then

    !If only one frame, then we don't have elastic parameters
    if (frames > 1) then
      allocate(OV_BPP(2,6), OV_BSP(2,9,nbp-1), OV_E_E(2,nbp-1), OV_C_l(2,nbp-1), &
             & OV_added(2,3,nbp-1), OV_V_E_E(2,nbp-1), OV_pV_E_E(2,nbp-1), & 
             & OV_F(2,10,nbp-1), OV_Ad2(2,nbp-1), OV_BSP_avstr(2,9,nbp-1), stat=ierror)
    else
      allocate(OV_BPP(2,6), OV_BSP(2,9,nbp-1), OV_E_E(2,nbp-1), OV_C_l(2,nbp-1), &
             & OV_added(2,3,nbp-1), OV_BSP_avstr(2,9,nbp-1), stat=ierror)
    end if

  else
    !As in double stranded structures, if only one frame, then we don't have elastic parameters
    !The only difference from double stranded structure, is that there are not BPP
    if (frames > 1) then
      allocate(OV_BSP(2,9,nbp-1), OV_E_E(2,nbp-1), OV_C_l(2,nbp-1), &
             & OV_added(2,3,nbp-1), OV_V_E_E(2,nbp-1), OV_pV_E_E(2,nbp-1), & 
             & OV_F(2,10,nbp-1), OV_Ad2(2,nbp-1), OV_BSP_avstr(2,9,nbp-1), stat=ierror)
    else
      allocate(OV_BSP(2,9,nbp-1), OV_E_E(2,nbp-1), OV_C_l(2,nbp-1), &
             & OV_added(2,3,nbp-1),  OV_BSP_avstr(2,9,nbp-1), stat=ierror)
    end if
  end if

  if (ierror/=0) stop "Error in allocating average parameters"

  !Calculate overall length averages and standard deviations

  !BPP
  if (strands ==2) then
    do i=1,6
      OV_BPP(:,i) = average_std(avstd_BPP(1,i,:), nbp)
    end do
  end if  

  !BSP
  l=1
  if (str /=2) then
  !Linear
    do j=1,nbp-1
      do i=1,9
        OV_BSP(:,i,j) = average_std(avstd_BSP(1,i, l:l-1+nbp-j ), nbp-j)
      end do
      l = l+nbp-j
    end do
  else
  !Closed
    do j=1,nbp-1
      do i=1,9
        OV_BSP(:,i,j) = average_std(avstd_BSP(1,i, l:l-1+nbp ), nbp)
      end do
      l = l+nbp
    end do
  end if
 
  !End to end distance
  l=1
  if (str /=2) then
  !Linear
    do j=1,nbp-1
      OV_E_E(:,j) = average_std(avstd_E_E(1, l:l-1+nbp-j ), nbp-j)
      l = l+nbp-j
    end do
  else
  !Closed
    do j=1,nbp-1
      OV_E_E(:,j) = average_std(avstd_E_E(1, l:l-1+nbp ), nbp)
      l = l+nbp
    end do
  end if

  !Contour length
  l=1
  if (str /=2) then
  !Linear
    do j=1,nbp-1
      OV_C_l(:,j) = average_std(avstd_C_l(1, l:l-1+nbp-j ), nbp-j)
      l = l+nbp-j
    end do 
  else
  !Closed
    do j=1,nbp-1
      OV_C_l(:,j) = average_std(avstd_C_l(1, l:l-1+nbp ), nbp)
      l = l+nbp
    end do 
  end if

  !Added parameters
  l=1
  if (str /=2) then
  !Linear
    do j=1,nbp-1
      do i=1,3
        OV_added(:,i,j) = average_std(avstd_added(1,i, l:l-1+nbp-j ), nbp-j)
      end do
      l = l+nbp-j
    end do
  else
  !Closed
    do j=1,nbp-1
      do i=1,3
        OV_added(:,i,j) = average_std(avstd_added(1,i, l:l-1+nbp ), nbp)
      end do
      l = l+nbp
    end do
  end if
 
  !If only one frame, we don't have elastic parameters 
  if (frames > 1) then

    !Variance and partial Variance of End to end distance
    l=1
    if (str /=2) then
    !Linear
      do j=1,nbp-1
        OV_V_E_E(:,j) = average_std(V(4,4, l:l-1+nbp-j ), nbp-j)
        OV_pV_E_E(:,j) = average_std(1.0_dp/V__1(4,4, l:l-1+nbp-j ), nbp-j)
        l = l+nbp-j
      end do 
    else
    !Closed
      do j=1,nbp-1
        OV_V_E_E(:,j) = average_std(V(4,4, l:l-1+nbp ), nbp)
        OV_pV_E_E(:,j) = average_std(1.0_dp/V__1(4,4, l:l-1+nbp ), nbp)
        l = l+nbp
      end do 
    end if
 
    !Elastic coeficients F
    l=1
    if (str /=2) then
    !Linear
      do j=1,nbp-1
        do i=1,10
          OV_F(:,i,j) = average_std( F(i, l:l-1+nbp-j ), nbp-j)
        end do
        l = l+nbp-j
      end do 
    else
    !Closed
      do j=1,nbp-1
        do i=1,10
          OV_F(:,i,j) = average_std( F(i, l:l-1+nbp ), nbp)
        end do
        l = l+nbp
      end do 
    end if

    !Dynamic persistence length (Ad2)
    l=1
    if (str /=2) then
    !Linear
      do j=1,nbp-1
        OV_Ad2(:,j) = average_std( Ad2(l:l-1+nbp-j ), nbp-j)
        l = l+nbp-j
      end do 
    else
    !Closed
      do j=1,nbp-1
        OV_Ad2(:,j) = average_std( Ad2(l:l-1+nbp ), nbp)
        l = l+nbp
      end do 
    end if

  end if !close frames > 1

  !BSP of average structure
  l=1
  if (str /=2) then
  !Linear
    do j=1,nbp-1
      do i=1,9
        OV_BSP_avstr(:,i,j) = average_std( BSP_avstr(i, l:l-1+nbp-j ), nbp-j)
      end do
      l = l+nbp-j
    end do
  else
  !Closed
    do j=1,nbp-1
      do i=1,9
        OV_BSP_avstr(:,i,j) = average_std( BSP_avstr(i, l:l-1+nbp ), nbp)
      end do
      l = l+nbp
    end do
  end if

  !WRITING SECTION
  !----------------------------------------------------------------------------------- 
  write(6,*) "Writing output files"

  if (strands == 2) then
    call write_BPP(avstd_BPP, OV_BPP, seq, nbp, frames, strands, str)
  end if

  call write_BSP(avstd_BSP, OV_BSP, seq, nbp, n_bsp, frames, strands, str)

  call write_structural(avstd_BSP, avstd_added, avstd_E_E, avstd_C_l, BSP_avstr, &
                      & OV_BSP, OV_added, OV_E_E, OV_C_l, OV_BSP_avstr, &
                      & seq, nbp, n_bsp, frames, strands, str)
 
  !As always, one frame then no elastic parameters
  if (frames > 1) then
    call write_elastic_parms(F, Ad2, V, V__1, OV_F, OV_Ad2, OV_V_E_E, OV_pV_E_E, &
                           & seq, nbp, n_bsp, frames, strands, str)
  end if

  !CLEANING SECTION
  !----------------------------------------------------------------------------------- 
  
  if (strands == 2) then
    deallocate(avstd_BPP, OV_BPP, stat=ierror)
    if (ierror /= 0) print*, "Couldn't deallocate overall BPP"
  end if
  
  deallocate( avstd_BSP, OV_BSP, avstd_added, avstd_E_E, avstd_C_L, &
            & BSP_avstr, OV_added, OV_E_E, OV_C_l, OV_BSP_avstr, stat=ierror )
  if (ierror /= 0) print*, "Couldn't deallocate overall structural parameters"
  
  if (frames > 1) then
    deallocate(F, Ad2, V, V__1, OV_F, OV_Ad2, OV_V_E_E, OV_pV_E_E, stat=ierror)
    if (ierror /= 0) print*, "Couldn't deallocate overall BPP"
  end if

  contains

  !----------------------------------------------------------------------------------- 
  !Calculates elastic matrix F (just to save lines in main section)
  !Values f elastic matrix are returned in nanometers (nm)
  !IV = Inverse covariance matrix.
  !N*bnm = rise
  function elastic_F(IV,N)
  implicit none
  real(dp) :: elastic_F(10), IV(4,4)
  integer :: N                      !bp length
  
  elastic_F(1)  = 1.0E23_dp*bKTpN*IV(4,4)*real(N,dp)*bnm       !Stretch
  elastic_F(2)  = IV(3,3)*real(N,dp)*bnm*rad_to_deg*rad_to_deg !Twist
  elastic_F(3)  = IV(2,2)*real(N,dp)*bnm*rad_to_deg*rad_to_deg !Roll
  elastic_F(4)  = IV(1,1)*real(N,dp)*bnm*rad_to_deg*rad_to_deg !Tilt
  elastic_F(5)  = IV(4,3)*real(N,dp)*bnm*rad_to_deg/w_0        !Stretch-Twist
  elastic_F(6)  = IV(4,2)*real(N,dp)*bnm*rad_to_deg/w_0        !Stretch-Roll
  elastic_F(7)  = IV(4,1)*real(N,dp)*bnm*rad_to_deg/w_0        !Stretch-Tilt
  elastic_F(8)  = IV(3,2)*real(N,dp)*bnm*rad_to_deg*rad_to_deg !Twist-Roll
  elastic_F(9)  = IV(3,1)*real(N,dp)*bnm*rad_to_deg*rad_to_deg !Twist-Tilt
  elastic_F(10) = IV(2,1)*real(N,dp)*bnm*rad_to_deg*rad_to_deg !Tilt-Roll

  end function elastic_F

  end program SerraNA
