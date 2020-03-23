!						31/07/2019
! SerraNA: Extract
! Extracts lengths from BSP, structural or elastic parameters.
! It can also calculate overalls (by length) in a region [a,b]
! Along the fragment
! It can also extract BPP, but the same values are printed (ready to plot)
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


  program Extract

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
  !sublength : length to extract data
  ![a,b] : ranges to obtain overalls
  !type_ext : type of extraction = 0, 1
  !type_parm : type of parameter to extract (BSP, elastic, structural) = 1,2,3
  integer :: nbp, frames, strands, str, n_bsp, ierror, &
           & sublength, a, b, &
           & type_ext, type_parm, &
           & i, j, l, N

  !file_in : input file (can be BSP, structural or elastic parameters)
  !cl, ca, cb -> auxiliar variables used in output file file_out
  !cF -> auxiliar format
  character(360) :: file_in, file_out, cl, ca, cb, cF

  !Even though the sequence would not matter, the subroutines read them, and 
  !anyway, we might need it in the future...
  !seq: sequence of strand I or II
  character(1), allocatable :: seq_I(:), seq_II(:)

  !subov_ ->  sub-overall, which is the overall values calculated from a to b
  !            all possible lengths are printed.
  !subl_ -> sub-length of interest
  !OV_ -> overalls read from input files
  !BPP -> Base-pair parameters
  !BSP -> Base-step parameters
  !elasp -> elastic parameters
  !strucp -> structural parameters
  !avstrp -> paremeters from the average structure
  !aux_parm -> Auxiliar parameters
  !mid -> mid point in fragment (between bp i and j)
  real(dp), allocatable :: BSP(:,:,:), OV_BSP(:,:,:), & !subl_BSP(:,:,:), subov_BSP(:,:,:), &
                         & elasp(:,:), OV_elasp(:,:,:), subl_elasp(:,:), subov_elasp(:,:,:), &
                         & strucp(:,:,:), avstrp(:,:), OV_strucp(:,:,:), OV_avstrp(:,:,:), &
                         & subl_strucp(:,:,:), subl_avstrp(:,:), &
                         & subov_strucp(:,:,:), subov_avstrp(:,:,:), &
                         & BPP(:,:,:), OV_BPP(:,:), aux_parm(:,:), aux_parm1(:), mid(:)

  real(dp) :: aux_r !auxiliar

  !---------------------EXTRACTION----------------------------------------------------
  !-----------------------------------------------------------------------------------

  !INPUT SECTION
  !-----------------------------------------------------------------------------------

  !Get input file, type of extraction and sublength or subfragment section
  call Extract_inputs( file_in, type_ext, sublength, a, b )

  !READING SECTION
  !-----------------------------------------------------------------------------------

  !Now read the input file (the whole thing).
  !The next subroutine figures out what type of data is going to be read
  call read_extract_input_file(file_in, type_parm, nbp, frames, strands, str, n_bsp, & 
                        & seq_I, seq_II, BSP, OV_BSP, elasp, OV_elasp, &
                        & strucp, avstrp, OV_strucp, OV_avstrp, BPP, OV_BPP)

  if (allocated(BPP) ) type_ext = 0 !We will output the parameters in this condition
  if (allocated(BSP) ) type_ext = 0 !Same for BSP

  !CHECK INPUT DATA, ALLOCATE AND WARNINGS
  !-----------------------------------------------------------------------------------

  !First, let's check for the case extracting a particular length
  !--------------------------------------------------------------
  if (type_ext == 0) then

    !Warnings!
    if (type_parm /= 4 .and. type_parm /=1 ) then !For not BPP.out
      if (sublength < 2) stop "Error, sublength must be greater than 1"
      if (sublength > nbp) stop "Error, sublength must not be greater than number of bps"

      write(6,*) "Extracting sublength", sublength
    end if
 
    !Now, decide what are we going to allocate
    select case ( type_parm )

    !BSP
    !--------------
    case(1)

      a = 1 !Nothing to do here
      sublength = 1 !Actually, do something. This is dinucleotide step
!I won't erase it just in case...
      !Linear
      !-----------------
!      if (str == 1) then

!        allocate(subl_BSP(2,7,nbp-sublength), stat=ierror)
!        if (ierror /= 0) stop "Error in allocating sub-lengths" 

      !Circular
      !-----------------
!      else 

!        allocate(subl_BSP(2,7,nbp), stat=ierror)
!        if (ierror /= 0) stop "Error in allocating sub-lengths" 

!      end if


    !Elastic
    !--------------
    case(2)


      !Linear
      !-----------------
      if (str == 1) then

        allocate(subl_elasp(13,nbp-sublength), stat=ierror)
        if (ierror /= 0) stop "Error in allocating sub-lengths" 

      !Circular
      !-----------------
      else 
 
        allocate(subl_elasp(13,nbp), stat=ierror)
        if (ierror /= 0) stop "Error in allocating sub-lengths" 

      end if


    !Structural
    !--------------
    case(3)


      !Linear
      !-----------------
      if (str == 1) then

        allocate(subl_strucp(2,11,nbp-sublength), stat=ierror)
        if (ierror /= 0) stop "Error in allocating sub-lengths" 

        allocate(subl_avstrp(3,nbp-sublength), stat=ierror)
        if (ierror /= 0) stop "Error in allocating sub-lengths" 

      !Circular
      !-----------------
      else 
 
        allocate(subl_strucp(2,11,nbp), stat=ierror)
        if (ierror /= 0) stop "Error in allocating sub-lengths" 

        allocate(subl_avstrp(3,nbp), stat=ierror)
        if (ierror /= 0) stop "Error in allocating sub-lengths" 

      end if

    !BPP
    !--------------
    case(4)

      a = 1 !Nothing to do here

    case default !Impossible

      stop "Different data from BSP, elas, struc? How could this be?"
    end select 


    !There is no need to allocate anything if we want to extract a subfragment

  !--------------------------------------------------------------
  else !Note that we already checked that type_ext is eather 0 or 1
  !--------------------------------------------------------------

    !Warnings!
    
    if (str == 1) then  !If linear...

      if (a == b .and. b /= 0 .and. a /= 0 ) stop "For opened structures: 0 < a, 0 < b or b == a == 0"
      if (a < 0 ) stop "For opened structures: 0 < a, 0 < b or b == a == 0" 
      if (b < 0 ) stop "For opened structures: 0 < a, 0 < b or b == a == 0"
      if (b < a ) stop "For opened structures: 0 < a, 0 < b or b == a == 0"
      if (a > nbp ) stop "For opened structures: 0 < a, 0 < b or b == a == 0"
      if (b > nbp ) stop "For opened structures: 0 < a, 0 < b or b == a == 0"
 
    else if (str == 2) then  !If closed...

      if ( a == b .and. b /= 0 .and. a /= 0 ) stop "For closed structures: if a = b, then a = b =0"
      if ( a < 0 ) stop  "For closed structures: 0 < a and 0 < b"
      if ( b < 0 ) stop  "For closed structures: 0 < a and 0 < b"
      if ( a > nbp ) stop  "For closed structures: a <= nbp and b <= nbp"
      if ( b > nbp ) stop  "For closed structures: a <= nbp and b <= nbp"
      if ( a == 0 .and. b > 0 ) stop "For closed structures: if a = b, then a = b = 0"
      if ( b == 0 .and. b > 0 ) stop "For closed structures: if a = b, then a = b = 0"

    else !Impossible case...

      stop "Not linear or closed? How could that be?"

    end if

    if ( a==0 .and. b==0) then
      write(6,*) "Writing overalls considering the whole fragment"
    else
      write(6,*) "Calculating overalls from subfragment [",a,b,"]"
    end if

    !Now, decide what are we going to allocate
    select case ( type_parm )

    !BSP
    !--------------
    case(1)
   
      a = 1 !Nothing to do here

      sublength = 1 !Actually, do something
      !Linear
      !-----------------
!      if (str == 1) then

!        allocate(subov_BSP(2,7,b-a), stat=ierror)
!        if (ierror /= 0) stop "Error in allocating sub-overalls" 

      !Circular
      !-----------------
!      else 
        
!        if (b ==0 .and. a == 0 ) then    !In case we want all the overalls

!          allocate(subov_BSP(2,7,nbp-1), stat=ierror)
!          if (ierror /= 0) stop "Error in allocating sub-overalls"

!        else                             !This condition covers both a>b and b>a

!          allocate(subov_BSP(2, 7, abs(b-a) ), stat=ierror)
!          if (ierror /= 0) stop "Error in allocating sub-overalls"

!        end if

!      end if


    !Elastic
    !--------------
    case(2)


      !Linear
      !-----------------
      if (str == 1) then

        allocate(subov_elasp(2,13,b-a), stat=ierror)
        if (ierror /= 0) stop "Error in allocating sub-overalls" 

      !Circular
      !-----------------
      else 
        
        if (b ==0 .and. a == 0 ) then    !In case we want all the overalls

          allocate(subov_elasp(2,13,nbp-1), stat=ierror)
          if (ierror /= 0) stop "Error in allocating sub-overalls"

        else                             !This condition covers both a>b and b>a

          allocate(subov_elasp(2, 13, abs(b-a) ), stat=ierror)
          if (ierror /= 0) stop "Error in allocating sub-overalls"

        end if

      end if


    !Structural
    !--------------
    case(3)


      !Linear
      !-----------------
      if (str == 1) then

        allocate(subov_strucp(2,11,b-a), stat=ierror)
        if (ierror /= 0) stop "Error in allocating sub-overalls" 

        allocate(subov_avstrp(2,3,b-a), stat=ierror)
        if (ierror /= 0) stop "Error in allocating sub-overalls" 

      !Circular
      !-----------------
      else 
        
        if (b ==0 .and. a == 0 ) then    !In case we want all the overalls

          allocate(subov_strucp(2,11,nbp-1), stat=ierror)
          if (ierror /= 0) stop "Error in allocating sub-overalls"

          allocate(subov_avstrp(2,3,nbp-1), stat=ierror)
          if (ierror /= 0) stop "Error in allocating sub-overalls"

       else                             !This condition covers both a>b and b>a

          allocate(subov_strucp(2, 11, abs(b-a) ), stat=ierror)
          if (ierror /= 0) stop "Error in allocating sub-overalls"

          allocate(subov_avstrp(2, 3, abs(b-a) ), stat=ierror)
          if (ierror /= 0) stop "Error in allocating sub-overalls"

        end if

      end if

    !BPP
    !--------------
    case(4)

      a = 1 !Nothing to do here

    case default !Impossible

      stop "Different data from BSP, elas, struc? How could this be?"
    end select 


  end if  
  ! End of type_ext
  !--------------------------------------------------------------
 
  !EXTRACTING SECTION  
  !-----------------------------------------------------------------------------------

  !If we want to extract a particular sublength
  if (type_ext == 0) then

    !What type of parameter?
    select case ( type_parm )

!Nothing happens if case(1)
    !BSP
    !--------------
    case(1)

      a = 1 !Do nothing, but just to make me feel that I'm doing something
      !This subroutine will figure out if linear or circular   
!      call extract_sublength_3d( BSP, subl_BSP, nbp, sublength, str )

    !Elastic
    !--------------
    case(2)

      call extract_sublength_2d( elasp, subl_elasp, nbp, sublength, str )

    !Structural
    !--------------
    case(3)

      !structural
      call extract_sublength_3d( strucp, subl_strucp, nbp, sublength, str )

      !from average structure
      call extract_sublength_2d( avstrp, subl_avstrp, nbp, sublength, str )

    end select 

  !--------------------------------------------------------------
  else !Note that we already checked that type_ext is eather 0 or 1
  !--------------------------------------------------------------
  !Here, central_fragment_mean_std does all the job. It figures out
  !the ranges and everything for us.

    !Now, type of parameter
    select case ( type_parm )

    !BSP
    !--------------
    case(1)
   
      a = 1 !Do nothing, but just to make me feel that I'm doing something

      !In the first two cases, there is no need to do anything, since
      !the we would only return values that we already have...
!      if ( a == 0 .and. b == 0) then
!        subov_BSP = OV_BSP
!      else if ( str==1 .and. a==1 .and. b==nbp ) then
!        subov_BSP = OV_BSP
!      else
!        do i = 1,7
!          call central_fragment_mean_std( subov_BSP(1,i,:), subov_BSP(2,i,:), &
!                                        & BSP(1,i,:), a, b, nbp, str )
!        end do
!      end if

    !Elastic
    !--------------
    case(2)
   
      !In the first two cases, there is no need to do anything, since
      !the we would only return values that we already have...
      if (a == 0 .and. b == 0) then
        subov_elasp = OV_elasp
      else if ( str==1 .and. a==1 .and. b==nbp ) then
        subov_elasp = OV_elasp
      else
        do i = 1,13
          call central_fragment_mean_std( subov_elasp(1,i,:), subov_elasp(2,i,:), &
                                        & elasp(i,:), a, b, nbp, str )
        end do
      end if

    !Structural
    !--------------
    case(3)


      !In the first two cases, there is no need to do anything, since
      !the we would only return values that we already have...
      if ( a == 0 .and. b == 0) then
        subov_strucp = OV_strucp
        subov_avstrp = OV_avstrp
      else if ( str==1 .and. a==1 .and. b==nbp ) then
        subov_strucp = OV_strucp
        subov_avstrp = OV_avstrp
      else
        do i = 1,11  !structural parms
          call central_fragment_mean_std( subov_strucp(1,i,:), subov_strucp(2,i,:), &
                                        & strucp(1,i,:), a, b, nbp, str )
        end do
        do i = 1,3  !average structure parms
          call central_fragment_mean_std( subov_avstrp(1,i,:), subov_avstrp(2,i,:), &
                                        & avstrp(i,:), a, b, nbp, str ) 
        end do
      end if


    end select 

  end if  
  ! End of type_ext
 
  !SORTING SECTION  (only performed for sublengths [if type_ext == 0])
  !-----------------------------------------------------------------------------------
  if (type_ext == 0) then

    !If linear
    if (str == 1) then
      N = nbp - sublength
    !If circular
    else
      N = nbp
    end if

    !What type of parameter?
    select case ( type_parm )

    !BSP
    !--------------
    case(1)
 
      allocate(aux_parm( 2, size(BSP,2) ), mid(N), stat=ierror )
      if (ierror /= 0) stop "Error in allocating mid point" 

      !Let's rearrange (this works for both circular and linear
      do i=1,N

        !Calculate mid point
        mid(i) = real( 2*i + sublength, dp) / 2.0_dp
        if (  mid(i) - real(nbp,dp)  > eps ) mid(i) = mid(i) - real(nbp,dp)

      end do

      !NOTE: Previous version was subl_BSP instead or BSP
      !Sort data
      do j = N-1, 1, -1
        do i = 1, j
          if (mid(i) > mid(i+1)) then
            aux_r = mid(i)
            aux_parm = BSP(:,:,i)
            mid(i) = mid(i+1)
            BSP(:,:,i) = BSP(:,:,i+1)
            mid(i+1) = aux_r
            BSP(:,:,i+1) = aux_parm
          end if
        end do
      end do

    !Elastic
    !--------------
    case(2)

      allocate(aux_parm1( size(subl_elasp,1) ), mid(N), stat=ierror )
      if (ierror /= 0) stop "Error in allocating mid point" 

      !Let's rearrange (this works for both circular and linear
      do i=1,N

        !Calculate mid point
        mid(i) = real( 2*i + sublength, dp) / 2.0_dp
        if (  mid(i) - real(nbp,dp)  > eps ) mid(i) = mid(i) - real(nbp,dp)

      end do

      !Sort data
      do j = N-1, 1, -1
        do i = 1, j
          if (mid(i) > mid(i+1)) then
            aux_r = mid(i)
            aux_parm1 = subl_elasp(:,i)
            mid(i) = mid(i+1)
            subl_elasp(:,i) = subl_elasp(:,i+1)
            mid(i+1) = aux_r
            subl_elasp(:,i+1) = aux_parm1
          end if
        end do
      end do


    !Structural
    !--------------
    case(3)

      allocate(aux_parm( 2, size(subl_strucp,2) ), aux_parm1(size(subl_avstrp,1)), mid(N), stat=ierror )
      if (ierror /= 0) stop "Error in allocating mid point" 

      !Let's rearrange (this works for both circular and linear
      do i=1,N

        !Calculate mid point
        mid(i) = real( 2*i + sublength, dp) / 2.0_dp
        if (  mid(i) - real(nbp,dp)  > eps ) mid(i) = mid(i) - real(nbp,dp)

      end do

      !Sort data for strucp
      do j = N-1, 1, -1
        do i = 1, j
          if (mid(i) > mid(i+1)) then
            aux_r = mid(i)
            aux_parm = subl_strucp(:,:,i)
            mid(i) = mid(i+1)
            subl_strucp(:,:,i) = subl_strucp(:,:,i+1)
            mid(i+1) = aux_r
            subl_strucp(:,:,i+1) = aux_parm
          end if
        end do
      end do

      !Sort data for avstrp
      do j = N-1, 1, -1
        do i = 1, j
          if (mid(i) > mid(i+1)) then
            aux_r = mid(i)
            aux_parm1 = subl_avstrp(:,i)
            mid(i) = mid(i+1)
            subl_avstrp(:,i) = subl_avstrp(:,i+1)
            mid(i+1) = aux_r
            subl_avstrp(:,i+1) = aux_parm1
          end if
        end do
      end do

    end select 

  end if  !type_ext ==0

  !OUTPUT SECTION  
  !-----------------------------------------------------------------------------------
  !Print sublengths
  if (type_ext == 0) then

    !Let's prepare the format before selecting a case
    write(cl,*) sublength
    cl = adjustl(cl)

    !What type of parameter?
    select case ( type_parm )
 
    !BSP
    !--------------
    case(1)

      !Commented lines in case we regret and want to use this method
!      !Let's first prepare the output file
!      ! 1 digit
!      if (sublength < 10) then
!        write(file_out,"(A14,I1,A4)") "BSP_sublength_", sublength, ".out"
!      ! 2 digits
!      else if (sublength > 9 .and. sublength < 100) then
!      write(file_out,"(A14,I2,A4)") "BSP_sublength_", sublength, ".out"
!      ! 3 digits
!      else if (sublength > 99 .and. sublength < 1000) then
!      write(file_out,"(A14,I3,A4)") "BSP_sublength_", sublength, ".out"
!      ! 4 difits
!      else if (sublength > 999 .and. sublength < 10000) then
!      write(file_out,"(A14,I4,A4)") "BSP_sublength_", sublength, ".out"
!      ! too big
!      else
!        stop "sublength too big, modify Extract.f90"
!      end if

      !Let's first prepare the output file and format
      file_out = "BSP_plot.out"
      do i=1,50
        if (F_BSP_2(i:i) .eq. "F" .or. F_BSP_2(i:i) .eq. "E") then
          l = i - 2
          cF =  "(1F10.1,"//F_BSP_2(l:)
          exit
        end if
      end do
     
      open(unit=10, file=trim(file_out), status="replace",action="write",iostat=ierror)
      if ( ierror/=0 ) stop "Error in opening output file"

        !Write!
        !------------------
                                            !Again, previous version was subl_BSP
        !Linear------------------

        if (str == 1 ) then
          do i = 1,nbp-sublength
            write(10,trim(cF)) mid(i), BSP(:,:,i)
          end do
        !Circular----------------
        else
          do i = 1,nbp
            write(10,trim(cF)) mid(i), BSP(:,:,i)
          end do
        end if

      close(unit=10,iostat=ierror)
      if ( ierror/=0 ) stop "Error in closing output file"
 
    !Elastic
    !--------------
    case(2)

      !Let's first prepare the output file
      file_out = "elastic_"//trim(cl)//"mer.out"
      do i=1,50
        if (F_ELAP_2(i:i) .eq. "F" .or. F_ELAP_2(i:i) .eq. "E") then
          l = i - 2
          cF =  "(1F10.1,"//F_ELAP_2(l:)
          exit
        end if
      end do
 
      open(unit=10, file=trim(file_out), status="replace",action="write",iostat=ierror)
      if ( ierror/=0 ) stop "Error in opening output file"

        !Write!
        !------------------
         
        !Linear------------------
        if (str == 1 ) then
          do i = 1,nbp-sublength
            write(10,trim(cF)) mid(i), subl_elasp(:,i)
          end do
        !Circular----------------
        else
          do i = 1,nbp
            write(10,trim(cF)) mid(i), subl_elasp(:,i)
          end do
        end if

      close(unit=10,iostat=ierror)
      if ( ierror/=0 ) stop "Error in closing output file"

    !Structural
    !--------------
    case(3)

      !Let's first prepare the output file
      file_out = "structural_"//trim(cl)//"mer.out"
      do i=1,50
        if (F_STRP_2(i:i) .eq. "F" .or. F_STRP_2(i:i) .eq. "E") then
          l = i - 2
          cF =  "(1F10.1,"//F_STRP_2(l:)
          exit
        end if
      end do
 
      open(unit=10, file=trim(file_out), status="replace",action="write",iostat=ierror)
      if ( ierror/=0 ) stop "Error in opening output file"

        !Write!
        !------------------
         
        !Linear------------------
        if (str == 1 ) then
          do i = 1,nbp-sublength
            write(10,trim(cF)) mid(i), subl_strucp(:,:,i), subl_avstrp(:,i)
          end do
        !Circular----------------
        else
          do i = 1,nbp
            write(10,trim(cF)) mid(i), subl_strucp(:,:,i), subl_avstrp(:,i)
          end do
        end if

      close(unit=10,iostat=ierror)
      if ( ierror/=0 ) stop "Error in closing output file"

    !BPP
    !--------------
    case(4)

      !Let's first prepare the output file and format
      file_out = "BPP_plot.out"
      do i=1,50
        if (F_BPP_2(i:i) .eq. "F" .or. F_BPP_2(i:i) .eq. "E") then
          l = i - 2
          cF =  "(I5,"//F_BPP_2(l:)
          exit
        end if
      end do
     
      open(unit=10, file=trim(file_out), status="replace",action="write",iostat=ierror)
      if ( ierror/=0 ) stop "Error in opening output file"

        !Write
        do i = 1,nbp
          write(10,trim(cF)) i, BPP(:,:,i)
        end do

      close(unit=10,iostat=ierror)
      if ( ierror/=0 ) stop "Error in closing output file"



     end select 

  !--------------------------------------------------------------
  else !If printing overalls...
  !--------------------------------------------------------------

    !Prepare a bit of format before deciding the case
    write(ca,*) a
    write(cb,*) b
    ca = adjustl(ca)
    cb = adjustl(cb)
 
    !Type of parameter
    select case ( type_parm )

    !BSP
    !--------------
!    case(1)

      !Let's first prepare the output file and format
!      if (a == 0 .and. b == 0) then
!        file_out = "ov_BSP_complete.out"
 !     else
 !       file_out = "ov_BSP_["//trim(ca)//":"//trim(cb)//"].out"
 !     end if 

 !     do i=1,50
 !       if (F_BSP_3(i:i) .eq. "F" .or. F_BSP_3(i:i) .eq. "E") then
 !         l = i - 2
 !         cF =  "(I5,"//F_BSP_3(l:)
 !         exit
 !       end if
 !     end do
     
 !     open(unit=10, file=trim(file_out), status="replace",action="write",iostat=ierror)
 !     if ( ierror/=0 ) stop "Error in opening output file"

        !Write!
        !------------------
        
 !       !Only one condition
 !       if (a == 0 .and. b == 0) then
 !         do i = 1, nbp-1
 !           write(10,trim(cF) ) i, subov_BSP(:,:,i)
 !         end do
 !       else                !This will work for any other condition
 !         do i = 1, abs(b-a)
 !           write(10,trim(cF) ) i, subov_BSP(:,:,i)
 !         end do
 !       end if
!
 !     close(unit=10,iostat=ierror)
 !     if ( ierror/=0 ) stop "Error in closing output file"
 
    !Elastic
    !--------------
    case(2)

      !Let's first prepare the output file and format
      if (a == 0 .and. b == 0) then
        file_out = "elastic_plot.out"
      else
        file_out = "elastic_["//trim(ca)//":"//trim(cb)//"].out"
      end if 

      do i=1,50
        if (F_ELAP_3(i:i) .eq. "F" .or. F_ELAP_3(i:i) .eq. "E") then
          l = i - 2
          cF =  "(I5,"//F_ELAP_3(l:)
          exit
        end if
      end do
     
      open(unit=10, file=trim(file_out), status="replace",action="write",iostat=ierror)
      if ( ierror/=0 ) stop "Error in opening output file"

        !Write!
        !------------------
        
        !Only one condition
        if (a == 0 .and. b == 0) then
          do i = 1, nbp-1
!            write(10,* ) i, subov_elasp(:,:,i)
            write(10,trim(cF) ) i+1, subov_elasp(:,:,i)
          end do
        else                !This will work for any other condition
          do i = 1, abs(b-a)
!            write(10,* ) i, subov_elasp(:,:,i)
            write(10,trim(cF) ) i+1, subov_elasp(:,:,i)
          end do
        end if

      close(unit=10,iostat=ierror)
      if ( ierror/=0 ) stop "Error in closing output file"
   
    !Structural
    !--------------
    case(3)

      !Let's first prepare the output file and format
      if (a == 0 .and. b == 0) then
        file_out = "structural_plot.out"
      else
        file_out = "structural_["//trim(ca)//":"//trim(cb)//"].out"
      end if 

      do i=1,50
        if (F_STRP_3(i:i) .eq. "F" .or. F_STRP_3(i:i) .eq. "E") then
          l = i - 2
          cF =  "(I5,"//F_STRP_3(l:)
          exit
        end if
      end do
     
      open(unit=10, file=trim(file_out), status="replace",action="write",iostat=ierror)
      if ( ierror/=0 ) stop "Error in opening output file"

        !Write!
        !------------------
        
        !Only one condition
        if (a == 0 .and. b == 0) then
          do i = 1, nbp-1
            write(10,trim(cF) ) i+1, subov_strucp(:,:,i), subov_avstrp(:,:,i)
          end do
        else                !This will work for any other condition
          do i = 1, abs(b-a)
            write(10,trim(cF) ) i+1, subov_strucp(:,:,i), subov_avstrp(:,:,i)
          end do
        end if

      close(unit=10,iostat=ierror)
      if ( ierror/=0 ) stop "Error in closing output file"
 
    end select 

  end if  
 
  !That's it! nearly 700 lines just for this...
  !Bye!


  !---------------------END-OF-EXTRACTION---------------------------------------------
  !-----------------------------------------------------------------------------------

  end program Extract
