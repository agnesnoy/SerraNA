!------------------------------------------------------------- 10/09/2019
! This module contains functions and subroutines needed for reading
! and writing.
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

  module io_mod

  use parms
  implicit none

  !Everything public
!  private arrdims
!  public topology_amber,  &
!       & coordinates_amber_crd

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

  contains

!-----------------------------------------------------------------------
!  Reads topology file from amber topology file
!    Gets number bp (nbp), checks if there's a "box" and the sequence on
!    both strands

  subroutine topology_amber(nbp,n_atoms,box,seq,ring_index,top,str)
  implicit none
  integer ::  nbp, n_atoms, box, str                !if box == 1 then there is a box
  integer, allocatable :: res_index(:), seq_index(:), ring_index(:)
  integer :: i, j, ierror, n_res, iaux(10), nlines, l, s1, s2, f_lines, &
           & pur, pyr  !purines, pyrimidines and type of structure
  character(360) :: top, caux                 !caux helps in reading
  character(1), allocatable, intent(out) :: seq(:)
  character(4), allocatable :: atom_names(:), res_names(:)

  !Get dimension of file top. Here, "l" is a dummy variable
  call arrdims(top,f_lines,l)

! READING TOPOLOGY----------------------------------------------------------------
  !Open topology file
  open(unit=10, file=trim(top), status="old",action="read",iostat=ierror)
    if (ierror/=0) stop 'Error in opening topology file'

    do j=1,f_lines
      read(10,"(A)") caux

      !Escape statement: if three quarters of file have been readed, then it
      !  it would automatically exit the loop. Hopefully, this will never happend,
      !  hence, it only is an emergency escape before collapsing.
      if (j > 3*f_lines / 4) exit

      !FLAG POINTERS
      if ( trim(caux) == "%FLAG POINTERS" ) then

        do                                       !In case there are any comments
          read(10,"(A)") caux                              
          if (caux(1:7) == "%FORMAT" ) exit                 
        end do

        read(10, 1001) iaux                                 !Has dimension 10 
        n_atoms=iaux(1)                                     !n_atoms
        read(10, 1001) iaux  
        n_res=iaux(2)                                       !n_res
        read(10, 1001) iaux 
        box=iaux(8)                                         !box

        !If there are no atoms stop.
        if (n_atoms <= 0) stop "No atoms in topology"
        !In case of no residues the program will try to
        !identify them (it could go wrong)

        !Print info
        write(6,"(1A20,1I10)") "Number of atoms = ", n_atoms
        write(6,"(1A20,1I10)") "Number of residues = ", n_res
        write(6,"(1A20,1I10)") "BOX = ", box

        !allocate atom_names
        allocate(atom_names(n_atoms), stat=ierror)
        if(ierror/=0) stop 'Error in allocating atom_names'
        !Not point in allocating residues since we can try identify them 
      end if
      !End of POINTERS section 


      !FLAG ATOM_NAME
      if ( trim(caux) == "%FLAG ATOM_NAME" ) then

         do                                       !In case there are any comments
          read(10,"(A)") caux                              
          if (caux(1:7) == "%FORMAT" ) exit                 
        end do

        !calculate lines in ATOM_NAME section
        if ( mod(n_atoms,20) == 0 ) then
          nlines = n_atoms / 20                       !20 atoms per line
        else
          nlines = n_atoms / 20 + 1        
        end if

        !READ atoms names
        l=1
        if (nlines > 1) then                          !in case (n_atoms < 20) which is unlikely
          do i=1,nlines-1
            read(10, 1002) atom_names(l:19+l)
            l=l+20
          end do
        end if
        read(10, 1002) atom_names(l:n_atoms)
      end if
      !End of ATOM_NAME SECTION


      !FLAG RESIDUE_LABEL
      if ( trim(caux) == "%FLAG RESIDUE_LABEL") then

         do                                       !In case there are any comments
          read(10,"(A)") caux                              
          if (caux(1:7) == "%FORMAT" ) exit                 
        end do

        !calculate lines in RESIDUE_LABEL section
        if ( mod(n_res,20) == 0 ) then
          nlines = n_res / 20                        !20 residues per line
        else
          nlines = n_res / 20 + 1                 
        end if
        
        !allocate res_names (is a temporary variable)
        allocate(res_names(n_res), stat=ierror)
        if (ierror /= 0 ) stop "Error in allocating res_names"

        !Read residues names
        l=1
        if (nlines > 1) then                         !in case (nlines < 20) which is unlikely
          do i=1,nlines-1
            read(10, 1002) res_names(l:19+l)
            l=l+20
          end do
        end if
        read(10, 1002) res_names(l:n_res)
      end if
      !End of RESIDUE_LABEL SECTION
      
      !FLAG RESIDUE_POINTER
       if ( trim(caux) == "%FLAG RESIDUE_POINTER") then
 
         do                                       !In case there are any comments
          read(10,"(A)") caux                              
          if (caux(1:7) == "%FORMAT" ) exit                 
        end do

        !calculate lines in RESIDUE_POINTER section
        if ( mod(n_res,10) == 0 ) then
          nlines = n_res / 10                         !10 residues per line
        else
          nlines = n_res / 10 + 1                 
        end if
        
        !allocate res_index (is a temporary variable)
        allocate(res_index(n_res), stat=ierror)
        if (ierror /= 0 ) stop "Error in allocating res_names"

        !Read residues indices (start)
        l=1
        if (nlines > 1) then                          !in case (nlines < 10) which is unlikely
          do i=1,nlines-1
            read(10, 1001) res_index(l:9+l)
            l=l+10
          end do
        end if
        read(10, 1001) res_index(l:n_res)

        !!!!!!!!!!!!!!!!!We have all we need!!!!!!!!!!!!!!!
        exit        

      end if
      !End of RESIDUE POINTER SECTION 
      
    end do
  close(10, iostat=ierror)
  if (ierror /= 0) stop "Error in closing topology file"
! END OF READING -------------------------------------------------------


  !Check if data obtained from topology file is enough to go on
  if ( .not. allocated(atom_names) ) then
    stop "Atoms names could not be identified"
  end if 

! COUNT NBP -------------------------------------------------------
  !The structure will be treated as a double stranded structure or a single strandedi one.
  !In case of dsDNA, if the number of residues in both strands is not the same, then
  !this will be treated as an error. 

  nbp=0
  do i=1,n_res
    !Identify Adenine
    if (any( A_l(:) .eq. res_names(i) ) ) then
      nbp=nbp+1
    end if
    !Identify Guanine
    if (any( G_l(:) .eq. res_names(i) ) ) then
      nbp=nbp+1
    end if
    !Identify Cytosine
    if (any( C_l(:) .eq. res_names(i) ) ) then
      nbp=nbp+1
    end if
    !Identify Thymine
    if (any( T_l(:) .eq. res_names(i) ) ) then
      nbp=nbp+1
    end if
    !Identify Uracil
    if (any( U_l(:) .eq. res_names(i) ) ) then
      nbp=nbp+1
    end if
  end do   !close i
 
       
! GET SEQUENCE -------------------------------------------------------

  !Allocate sequence and sequence index
  allocate(seq(nbp),seq_index(nbp), stat=ierror) 
  if(ierror/=0) stop 'Error in allocating sequence'

  !Get the correct number of nbp and check if structure is incomplete in case of 
  !a double stranded structure
  if ( str == 2 ) then
    if ( mod(nbp,2) /= 0 ) stop "Incomplete double stranded structure"
    nbp = nbp/2                 !We counted bases from both strands
  end if  
  !if str=1 then its a single stranded structure. If its 2 then is double stranded
   
  !Check if molecule is longer than 4bp, if not then stop
  if ( nbp/2 .le. 4 ) stop "Invalid DNA fragment, the fragment must be larger than 4 bp"

  !Similar loop as before but now identifying residues
  l=0 !will help us count bps
  pur = 0 ; pyr = 0 !purines and pyrimidines counters
  do i=1,n_res
    !Identify Adenine
    if (any( A_l(:) .eq. res_names(i) ) ) then
      l=l+1
      seq(l) = "A" 
      seq_index(l) = res_index(i)
      pur = pur +1
    end if
    !Identify Guanine
    if (any( G_l(:) .eq. res_names(i) ) ) then
      l=l+1
      seq(l) = "G" 
      seq_index(l) = res_index(i) 
      pur = pur +1
    end if
    !Identify Cytosine
    if (any( C_l(:) .eq. res_names(i) ) ) then
      l=l+1
      seq(l) = "C" 
      seq_index(l) = res_index(i) 
      pyr = pyr +1
    end if
    !Identify Thymine
    if (any( T_l(:) .eq. res_names(i) ) ) then
      l=l+1
      seq(l) = "T" 
      seq_index(l) = res_index(i) 
      pyr = pyr +1
    end if
    !Identify Uracil
    if (any( U_l(:) .eq. res_names(i) ) ) then
      l=l+1
      seq(l) = "U" 
      seq_index(l) = res_index(i) 
      pyr = pyr +1
    end if
  end do   !close i
     
  !Check we read all correctly the sequence and indices
  if (l /= nbp*str) stop 'Error in extracting sequence and res pointers'

  write(6,*) "Base pairs read = ", nbp

! INDECES OF RING ATOMS -------------------------------------------------------

  !Allocate ring index
  allocate(ring_index(pyr*6+pur*9), stat=ierror)
  if(ierror/=0) stop 'Error in allocating ring_index'

  ring_index = 0 !initializing, this will be true for SerraLINE
  l = 0          !This will be an auxiliar

  do i=1,nbp*str !l is number of nucleotides or pyr+pur

    !Set boundaries
    s1 = seq_index(i)
    if  (i < nbp*str) then
      s2 = seq_index(i+1)
    else
      s2 = n_atoms
    end if

    !For purines
    if (seq(i) .eq. "A" .or. seq(i) .eq. "G" ) then
      do j=s1,s2
        if (atom_names(j) .eq. N9) then
          ring_index(l+1) = j
        end if
        if (atom_names(j) .eq. C8) then
          ring_index(l+2) = j
        end if
        if (atom_names(j) .eq. N7) then
          ring_index(l+3) = j
        end if
        if (atom_names(j) .eq. C5) then
          ring_index(l+4) = j
        end if
        if (atom_names(j) .eq. C6) then
          ring_index(l+5) = j
        end if
        if (atom_names(j) .eq. N1) then
          ring_index(l+6) = j
        end if
        if (atom_names(j) .eq. C2) then
          ring_index(l+7) = j
        end if
        if (atom_names(j) .eq. N3) then
          ring_index(l+8) = j
        end if
        if (atom_names(j) .eq. C4) then
          ring_index(l+9) = j
        end if
     end do
     l=l+9
    else 
      !If its pyrimidine
      do j=s1,s2
        if (atom_names(j) .eq. N1) then
          ring_index(l+1) = j
        end if
        if (atom_names(j) .eq. C2) then
          ring_index(l+2) = j
        end if
        if (atom_names(j) .eq. N3) then
          ring_index(l+3) = j
        end if
        if (atom_names(j) .eq. C4) then
          ring_index(l+4) = j
        end if
        if (atom_names(j) .eq. C5) then
          ring_index(l+5) = j
        end if
        if (atom_names(j) .eq. C6) then
          ring_index(l+6) = j
        end if
      end do
      l=l+6 
    end if
  end do

 1001 format (10I8)   !Format used in FLAG pointers and residue pointers
 1002 format (20A4)   !Format used in atom and residue names

  end subroutine topology_amber
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------

 !Reads the atom coordinates from amber trajectory file
 !Also, in this subroutine the first two and last two bp will be ignored,
 !if it is not a closed structure (c_str), so nbp will be modified (nbp-4)
 !Only works with complete structures

 subroutine coordinates_amber_crd(coords, frames, ring_index, n_atoms, seq, nbp, &
                                & traj, box, s_str, c_str)

 !coords will save the coordinates of ring atoms.
 !b and a_coords are auxiliar arrays for reading
 !N is the total number of atoms
 !r_atoms -> number of ring atoms
 real(dp), allocatable :: coords(:,:,:), a_coords(:,:), b(:)
 integer :: ring_index(:), n_atoms, nbp, box, s_str, c_str, &
          & nlin, ncol, frames, r_atoms, ierror, i, k, l, r, rl, &
          & rn, ends(2*s_str) !ends
 character(360) :: traj, aux
 character(1), allocatable :: seq(:), aseq(:) !aseq -> auxiliar sequence 

  !bp that will be ignored
  l=0 !auxiliars
  r=1
  do i=1,nbp*s_str
    if (seq(i) .eq. "G" .or. seq(i) .eq. "A") then
      l=l+9  !purine
    else
      l=l+6  !pyrimidine
    end if
    if (i==2 .or. i==nbp-2 .or. i==nbp+2 .or. i==nbp*2-2) then
      ends(r) = l+1
      r=r+1
    end if
  end do

  !correct ends if linear structure
  if (c_str /= 2) then
    ends(2) = ends(2) -1 
    r_atoms = ends(2)-ends(1)+1 !correct number of atoms
    if (s_str == 2) then
      ends(4) = ends(4) -1
      r_atoms = r_atoms+ends(4)-ends(3)+1
    end if
  else
    r_atoms = l
  end if

 !Correct new sequence if linear structure
 if (c_str /= 2) then
   allocate(aseq(nbp*s_str), stat=ierror)
   if (ierror /= 0) stop "Error in allocating auxiliar sequence"
 
   aseq=seq
   deallocate(seq, stat=ierror)
   if (ierror /= 0) stop "Error in deallocating sequence"

   allocate(seq( s_str*(nbp-4) ), stat=ierror)
   if (ierror /= 0) stop "Error in allocating sequence - 4bp"

   !rearrange sequence - 4bp
   do i=3,nbp-2
     seq(i-2)= aseq(i)
   end do 

   !if double structure
   if (s_str ==2) then
     do i=nbp+3, 2*nbp-2
       seq(i-6)= aseq(i)
     end do 
   end if

   !correct nbp 
   nbp = nbp - 4

 end if

 !READ TRAJECTORY
 !get dimensions of trajectory file
  call arrdims(traj,nlin,ncol)

  ncol=3*n_atoms/10 !Let's recycle ncol
  if (3*n_atoms-10*ncol > 0) then
    ncol = ncol+1   !the extra line
  end if

  !Get number of frames in file

  if (box /= 0) then           !Just in case let's put it here
    frames = nlin/(ncol+1)
  else
    frames = nlin/ncol
  end if

  write(6,"(1A20,1I10)") "Number of frames =", frames
  !Lets recycle nlin, it will be the number of lines that can be readed completly by 3
  nlin = n_atoms/10  

  !rl is the number of remaining lines
  rl=3*mod(n_atoms,10)
  if ( mod(rl,10) == 0 ) then
    rl = rl/10
  else
    rl = 1+ rl/10
  end if

  !rn is number of atoms remaining
  rn = n_atoms - 10*nlin    !Its 10 atoms per nlin (3 lines)

  allocate(coords(3,r_atoms,frames), a_coords(3,n_atoms), b(rn*3), stat=ierror)
  if (ierror /= 0) stop "Error in allocating coordinates array"

  open(unit=10, file=trim(traj), status="old",action="read",iostat=ierror)
    if (ierror/=0) stop 'Error in opening trajectory file'

  !Begin reading
  read(10,*) aux !skips first line

  do k=1,frames

  l = 0

    do i=1,nlin
      read(10, 1003) a_coords(1:3,l+1), a_coords(1:3,l+2), &
                   & a_coords(1:3,l+3), a_coords(1,l+4)

      read(10, 1003) a_coords(2:3,l+4), a_coords(1:3,l+5), &
                   & a_coords(1:3,l+6), a_coords(1:2,l+7)

      read(10, 1003) a_coords(3,l+7),   a_coords(1:3,l+8), &
                   & a_coords(1:3,l+9), a_coords(1:3,l+10)
      l = l+10
    end do !close i

    ! If there are more lines to read
    if (rl > 0) then
      !If there's one more
      if (rl >= 1) then
        if (rl < 2) then
          read(10, 1003) b(1:rn*3)
        else
          read(10, 1003) b(1:10)
        end if
      end if 
    
      !If there's more than one
      if (rl >= 2) then
        if (rl < 3) then
          read(10, 1003) b(11:rn*3)
        else
          read(10, 1003) b(11:20)
        end if
      end if 
    
      !If there's 3 left
      if (rl == 3) then
        read(10, 1003) b(21:rn*3)
      end if 
      
      !Now, collect the remaining coordinates
      do i=1,rn
        a_coords(1,l+i) = b((i-1)*3+1)
        a_coords(2,l+i) = b((i-1)*3+2)
        a_coords(3,l+i) = b((i-1)*3+3)
      end do

    end if !close rl>0 
    
    !Just in case that there is a box, skip this line
    if (box /= 0) then
      read(10,*) aux
    end if

    !Now, let's take only the atoms of interest
    if (c_str ==2) then 
      !if closed structure
      l = 1
      do i=1,r_atoms
        coords(:,l,k) = a_coords(:,ring_index(i))
        l = l+1
      end do

    else 

      !if open...
      l = 0
      do i=ends(1),ends(2)
        l = l+1
        coords(:,l,k) = a_coords(:,ring_index(i))
      end do

      if (s_str==2) then
        do i=ends(3),ends(4)
          l = l+1
          coords(:,l,k) = a_coords(:,ring_index(i))
        end do
      end if
    end if 

    !all if's within this loop are probably bad for performance,
    !but they are called once per frame so it is not that bad.
  end do !close k

  close(10) !close trajectory

 1003 format (10F8.3) !Format used in amber trajectory file 

 end subroutine coordinates_amber_crd
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! Obtains the lines and columns of the file "filen"
 subroutine arrdims(filen,nlin,ncol)
 character(260), intent(in) :: filen
 integer, intent(out) :: nlin,ncol
 character(260) :: fwc
 integer :: ierr

 fwc='.functionsmod_wc-ou.tmp'   ! auxilliary file, where we store number
                                   ! of lines and columns separated by newline

 call system('head -n 1 '//filen//' | wc -w  > '//trim(fwc),ierr)
 if (ierr == 0) then
   call system('cat '//filen//' | wc -l  >> '//trim(fwc),ierr)
   if (ierr == 0) then
     open(unit=10,file=trim(fwc),status='old',action='read')
     read(10,*) ncol
     read(10,*) nlin        ! We obtain the number of lines and number of
     close(10)              ! columns of the data file
     call system( 'rm ' //trim(fwc),  ierr )
     if (ierr /= 0) print*, "Couldn't erase temporary file ", trim(fwc)
   else                
     write(6,*) "Could not determine number of lines in file: ", trim(filen)                
     stop
   end if
 else
   write(6,*) "Could not determine number of columns in file: ", trim(filen)                

 end if

 if (max(nlin,ncol) == 0) then
   write(6,'(a)') 'Error in arrdims. Did not get valid dimensions'
   write(6,'(a,i10,a,i10)') '  nlin=',nlin,'   ncol=',ncol
   stop
 end if

 end subroutine arrdims

!-----------------------------------------------------------------------


!----------------------------------------------------------------------
  !Subroutine used in SerraNA
  !Writes base-pair parameters.
  !BPP contains averages and standard deviations of the six parameters.
  !OV_BPP is the overall values by bp length
  subroutine write_BPP(BPP, OV_BPP, seq, nbp, frames, strands, str)
  implicit none

  integer, intent(in) :: nbp, frames, strands, str
  character(1), intent(in) :: seq(strands*nbp)
  real(dp), intent(in) :: BPP(2,6,nbp), OV_BPP(2,6)

  integer :: i, ierror
  character(1) :: seq2(nbp) !second strand sequence

  if (strands==2) then
    do i=1,nbp
      seq2(i) = seq(2*nbp - i+1) !backwards, this is the second strand
    end do
  end if
  open(unit=10, file="BPP.out", iostat=ierror)
  if (ierror/=0) stop "Error in opening BPP output file"
    write(10,*) "BASE-PAIR PARAMETERS"
    write(10,*) ""
    if (str==2) then
      write(10,*) "CLOSED STRUCTURE"
    else
      write(10,*) "LINEAR STRUCTURE"
    end if 
    if (strands==2) then
      write(10,*) "DOUBLE-STRANDED STRUCTURE ANALYSED"
    else
      write(10,*) "SINGLE-STRANDED STRUCTURE ANALYSED"
    end if
    write(10,*) "BASE-PAIRS: ", nbp
    write(10,*) "FRAMES:     ", frames
    write(10,*) "SEQUENCE:   "
    write(10,*) "STRAND 1:   ", seq(1:nbp)
    if (strands==2) then
      write(10,*) "STRAND 2:   ", seq2(1:nbp)
    else
      write(10,*) "STRAND 2:   "
    end if
    write(10,*) ""
    write(10,*) "First column averages, second column standard deviations"
    write(10,trim(F_BPP_1) ) " base-pair ", "Shear", "Stretch", "Stagger", "Buckle", "Propeller", "Openning"

    write(10,*) "-----------------------------------------------------------------", &
              & "-----------------------------------------------------------------" 
    do i=1,nbp
      if (strands == 2) then 
        write(10, trim(F_BPP_2) ) i, " ", seq(i), "-", seq2(i), BPP(:,1,i), BPP(:,2,i), & 
                     & BPP(:,3,i), BPP(:,4,i), BPP(:,5,i), BPP(:,6,i) 
      else
        write(10,trim (F_BPP_2) ) i, " ", seq(i), " ", "", BPP(:,1,i), BPP(:,2,i), &
                     & BPP(:,3,i), BPP(:,4,i), BPP(:,5,i), BPP(:,6,i) 
      end if
    end do
    write(10,*) "-----------------------------------------------------------------", &
              & "-----------------------------------------------------------------"
    write(10,trim(F_BPP_3) ) " AVG-STD= ", OV_BPP(:,1),  OV_BPP(:,2),  OV_BPP(:,3),  OV_BPP(:,4),  &
                 & OV_BPP(:,5),  OV_BPP(:,6)

  close(10, iostat=ierror)
  if (ierror/=0) stop "Error in closing BPP output file"
  
  end subroutine write_BPP
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
  !Writes base-step parameters.
  !BSP contains the averages and std for each base-step
  !OV_BSP the averages and std of the averga of BSP by length.
  subroutine write_BSP(BSP, OV_BSP, seq, nbp, n_bsp, frames, strands, str)
  implicit none

  integer, intent(in) :: nbp, n_bsp, frames, strands, str
  character(1), intent(in) :: seq(strands*nbp)
  real(dp), intent(in) :: BSP(2,9,n_bsp), OV_BSP(2,9,nbp-1)

  integer :: i, j, l, w, ierror
  character(1) :: seq2(nbp) !second strand sequence

  if (strands==2) then
    do i=1,nbp
      seq2(i) = seq(2*nbp - i+1) !backwards, this is the second strand
    end do
  end if

  open(unit=10, file="BSP.out", iostat=ierror)
  if (ierror/=0) stop "Error in opening BSP output file"
    write(10,*) "BASE-STEP PARAMETERS"
    write(10,*) ""
    if (str==2) then
      write(10,*) "CLOSED STRUCTURE"
    else
      write(10,*) "LINEAR STRUCTURE"
    end if 
    if (strands==2) then
      write(10,*) "DOUBLE-STRANDED STRUCTURE ANALYSED"
    else
      write(10,*) "SINGLE-STRANDED STRUCTURE ANALYSED"
    end if
   write(10,*) "BASE-PAIRS: ", nbp
    write(10,*) "FRAMES:     ", frames
    write(10,*) "SEQUENCE:   "
    write(10,*) "STRAND 1:   ", seq(1:nbp)
    if (strands==2) then
      write(10,*) "STRAND 2:   ", seq2(1:nbp)
    else
      write(10,*) "STRAND 2:   "
    end if
    write(10,*) ""
    write(10,*) "First column averages, second column standard deviations"

    l = 0
    do j=1,1 !nbp-1
  !    write(10,*) j, "bp"
      write(10,trim(F_BSP_1)) " base-step ", "Shift", "Slide", "Rise", "Tilt", "Roll", "Twist", "Bending"
      write(10,*) "-----------------------------------------------------------------", &
                & "-----------------------------------------------------------------", &
                & "-------------------------" 

      !Write parameters
      if (str /= 2) then  !Linear or closed structure

        !Linear structure-------------------------------------------------------------------
        if (strands==2 ) then
          do i=1,nbp-j
            l=l+1
            write(10,trim(F_BSP_2)) i, "-", i+j, " ", seq(i), seq(i+j), "/", seq(2*nbp-j-i+1), & 
                         & seq(2*nbp+1-i), BSP(:,1,l), BSP(:,2,l), BSP(:,3,l), & 
                         & BSP(:,4,l), BSP(:,5,l), BSP(:,6,l), BSP(:,7,l)
          end do !close i
        else 
          do i=1,nbp-j
            l=l+1
            write(10,trim(F_BSP_2)) i, "-", i+j, " ", seq(i), seq(i+j), "/", "#", & 
                         & "#", BSP(:,1,l), BSP(:,2,l), BSP(:,3,l), & 
                         & BSP(:,4,l), BSP(:,5,l), BSP(:,6,l), BSP(:,7,l)
          end do !close i
        end if
      else
        !closed structure-------------------------------------------------------------------
        if (strands==2 ) then
          do i=1,nbp
            l=l+1

            !Check correct bp w (remember, we are on the other side when i > nbp-j)
            if (i > nbp-j) then
              w = i+j-nbp
            else
              w = i+j
            end if 

            write(10,trim(F_BSP_2)) i, "-", w, " ", seq(i), seq(w), "/", seq2(w), & 
                         & seq2(i), BSP(:,1,l), BSP(:,2,l), BSP(:,3,l), & 
                         & BSP(:,4,l), BSP(:,5,l), BSP(:,6,l), BSP(:,7,l)

          end do !close i

        else 
          do i=1,nbp
            l=l+1

            !Check correct bp w (remember, we are on the other side when i > nbp-j)
            if (i > nbp-j) then
              w = i+j
            else
              w = i+j-nbp
            end if 

            write(10,trim(F_BSP_2)) i, "-", w, " ", seq(i), seq(w), "/", "#", & 
                         & "#", BSP(:,1,l), BSP(:,2,l), BSP(:,3,l), & 
                         & BSP(:,4,l), BSP(:,5,l), BSP(:,6,l), BSP(:,7,l)
          end do !close i

        end if !Close if strands
 
      end if !Close IF str

      !Write overalls
      write(10,*) "-----------------------------------------------------------------", &
                & "-----------------------------------------------------------------", &
                & "-------------------------" 
     write(10,trim(F_BSP_3)) " AVG-STD = ", OV_BSP(:,1,j), OV_BSP(:,2,j), OV_BSP(:,3,j), &
                  & OV_BSP(:,4,j), OV_BSP(:,5,j), OV_BSP(:,6,j), OV_BSP(:,7,j)
      write(10,*) ""
    end do   !close j

  close(10, iostat=ierror)
  if (ierror/=0) stop "Error in closing BPP output file"
  

  end subroutine write_BSP
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
  !Write structural parameters.
  !First row are averages and std of BSP, added parameters, End to end 
  !distance, Contour length, and BSP of average structure.
  !Arrays with OV are the overall values
  subroutine write_structural(BSP, added, E_E, C_l, BSP_avstr, &
                           & OV_BSP, OV_added, OV_E_E, OV_C_l, OV_BSP_avstr, &
                           & seq, nbp, n_bsp, frames, strands, str)
  implicit none

  integer, intent(in) :: nbp, n_bsp, frames, strands, str
  character(1), intent(in) :: seq(strands*nbp)
  real(dp), intent(in) :: BSP(2,9,n_bsp), added(2,3,n_bsp), &
                        & E_E(2,n_bsp), C_l(2,n_bsp), BSP_avstr(9,n_bsp), &
                        & OV_BSP(2,9,nbp-1), OV_added(2,3,nbp-1), &
                        & OV_E_E(2,nbp-1), OV_C_l(2,nbp-1), &
                        & OV_BSP_avstr(2,9,nbp-1)

  integer :: i, j, l, w, ierror
  character(1) :: seq2(nbp) !second strand sequence

  if (strands==2) then
    do i=1,nbp
      seq2(i) = seq(2*nbp - i+1) !backwards, this is the second strand
    end do
  end if

  open(unit=10, file="structural_parameters.out", iostat=ierror)
  if (ierror/=0) stop "Error in opening structural parameters output file"
    write(10,*) "STRUCTURAL PARAMETERS"
    write(10,*) ""
    if (str==2) then
      write(10,*) "CLOSED STRUCTURE"
    else
      write(10,*) "LINEAR STRUCTURE"
    end if  
    if (strands==2) then
      write(10,*) "DOUBLE-STRANDED STRUCTURE ANALYSED"
    else
      write(10,*) "SINGLE-STRANDED STRUCTURE ANALYSED"
    end if
    write(10,*) "BASE-PAIRS: ", nbp
    write(10,*) "FRAMES:     ", frames
    write(10,*) "SEQUENCE:   "
    write(10,*) "STRAND 1:   ", seq(1:nbp)
    if (strands==2) then
      write(10,*) "STRAND 2:   ", seq2(1:nbp)
    else
      write(10,*) "STRAND 2:   "
    end if
    write(10,*) ""
    write(10,*) "First column averages, second column standard deviations"

    l = 0
    do j=1,nbp-1
      write(10,*) j+1, "mer"
      write(10,trim(F_STRP_1)) " base-step ", "Added Shift", "Added Slide", "Added Rise", &
                   & "End-to-End L", "Contour L", "Twist", "Roll", "Tilt", "Bending", &
                   &  "Bending**2", "D correlation", "AVSTR B", "AVSTR B**2", "AVSTR D C"
      write(10,*) "-----------------------------------------------------------------", &
                & "-----------------------------------------------------------------", &
                & "-----------------------------------------------------------------", &
                & "-----------------------------------------------------------------", &
                & "---------------" 

      !Write parameters
      if (str /= 2) then  !Linear or closed structure

        !Linear structure-------------------------------------------------------------------
        if (strands==2 ) then
          do i=1,nbp-j
            l=l+1
            write(10,trim(F_STRP_2)) i, "-", i+j, " ", seq(i), seq(i+j), "/", seq(2*nbp-j-i+1), & 
                         & seq(2*nbp+1-i), added(:,1,l), added(:,2,l), added(:,3,l), & 
                         & E_E(:,l), C_l(:,l), BSP(:,6,l), BSP(:,5,l), BSP(:,4,l), &
                         & BSP(:,7,l), BSP(:,9,l), BSP(:,8,l), BSP_avstr(7,l), BSP_avstr(9,l), &
                         & BSP_avstr(8,l)
          end do !close i
        else 
          do i=1,nbp-j
            l=l+1
            write(10,trim(F_STRP_2)) i, "-", i+j, " ", seq(i), seq(i+j), "/", "#", & 
                         & "#", added(:,1,l), added(:,2,l), added(:,3,l), & 
                         & E_E(:,l), C_l(:,l), BSP(:,6,l), BSP(:,5,l), BSP(:,4,l), &
                         & BSP(:,7,l), BSP(:,9,l), BSP(:,8,l), BSP_avstr(7,l), BSP_avstr(9,l), &
                         & BSP_avstr(8,l)
          end do !close i
        end if !Linear structure end -------------------------------------------------------

      else
        !Closed structure-------------------------------------------------------------------
        if (strands==2 ) then
          do i=1,nbp
            l=l+1

            !Check correct bp w (remember, we are on the other side when i > nbp-j)
            if (i > nbp-j) then
              w = i+j-nbp
            else
              w = i+j
            end if 

            write(10,trim(F_STRP_2)) i, "-", w, " ", seq(i), seq(w), "/", seq2(w), & 
                         & seq2(i), added(:,1,l), added(:,2,l), added(:,3,l), & 
                         & E_E(:,l), C_l(:,l), BSP(:,6,l), BSP(:,5,l), BSP(:,4,l), &
                         & BSP(:,7,l), BSP(:,9,l), BSP(:,8,l), BSP_avstr(7,l), &
                         & BSP_avstr(9,l), BSP_avstr(8,l)
 
          end do !close i

        else
          do i=1,nbp
            l=l+1

            !Check correct bp w (remember, we are on the other side when i > nbp-j)
            if (i > nbp-j) then
              w = i+j
            else
              w = i+j-nbp
            end if

            write(10,trim(F_STRP_2)) i, "-", w, " ", seq(i), seq(w), "/", "#", & 
                         & "#", added(:,1,l), added(:,2,l), added(:,3,l), & 
                         & E_E(:,l), C_l(:,l), BSP(:,6,l), BSP(:,5,l), BSP(:,4,l), &
                         & BSP(:,7,l), BSP(:,9,l), BSP(:,8,l), BSP_avstr(7,l), & 
                         & BSP_avstr(9,l), BSP_avstr(8,l)
 
          end do !close i

        end if !Closed structure end ------------------------------------------------------

      end if !Close IF str

      !Write overalls--------
      write(10,*) "-----------------------------------------------------------------", &
                & "-----------------------------------------------------------------", &
                & "-----------------------------------------------------------------", &
                & "-----------------------------------------------------------------", &
                & "---------------" 
     write(10,trim(F_STRP_3)) " AVG-STD = ", OV_added(:,1,j), OV_added(:,2,j), OV_added(:,3,j), &
                  & OV_E_E(:,j), OV_C_l(:,j), OV_BSP(:,6,j), OV_BSP(:,5,j), & 
                  & OV_BSP(:,4,j), OV_BSP(:,7,j), OV_BSP(:,9,j), OV_BSP(:,8,j), OV_BSP_avstr(:,7,j), &
                  & OV_BSP_avstr(:,9,j), OV_BSP_avstr(:,8,j)
      write(10,*) ""
    end do   !close j

  close(10, iostat=ierror)
  if (ierror/=0) stop "Error in closing BPP output file"
  
  end subroutine write_structural

!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
  !Subroutine used in SerraNA
  !Writes elastic parameters: Coeficients of elastic matrix F, End to 
  !end variance V(4,4,:), and partial variance 1/V__1(4,4,:).
  !OV_V_E_E only carries the overall end to end variance and OV_pV_E_E
  !the overalls of the end to end partial variance
  subroutine write_elastic_parms(F, Ad2, V, V__1, OV_F, OV_Ad2, OV_V_E_E, OV_pV_E_E, & 
                               & seq, nbp, n_bsp, frames, strands, str)
  implicit none

  integer, intent(in) :: nbp, n_bsp, frames, strands, str
  character(1), intent(in) :: seq(strands*nbp)
  real(dp), intent(in) :: F(10,n_bsp), Ad2(n_bsp), V(4,4,n_bsp), &
                        & V__1(4,4,n_bsp), OV_F(2,10,nbp-1), OV_Ad2(2,n_bsp), &
                        & OV_V_E_E(2,nbp-1), OV_pV_E_E(2,nbp-1)

  integer :: i, j, l, w, ierror
  character(1) :: seq2(nbp) !second strand sequence

  if (strands==2) then
    do i=1,nbp
      seq2(i) = seq(2*nbp - i+1) !backwards, this is the second strand
    end do
  end if

  open(unit=10, file="elastic_parameters.out", iostat=ierror)
  if (ierror/=0) stop "Error in opening elastic_parameters output file"
    write(10,*) "ELASTIC PARAMETERS"
    write(10,*) ""
    if (str==2) then
      write(10,*) "CLOSED STRUCTURE"
    else
      write(10,*) "LINEAR STRUCTURE"
    end if 
    if (strands==2) then
      write(10,*) "DOUBLE-STRANDED STRUCTURE ANALYSED"
    else
      write(10,*) "SINGLE-STRANDED STRUCTURE ANALYSED"
    end if
    write(10,*) "BASE-PAIRS: ", nbp
    write(10,*) "FRAMES:     ", frames
    write(10,*) "SEQUENCE:   "
    write(10,*) "STRAND 1:   ", seq(1:nbp)
    if (strands==2) then
      write(10,*) "STRAND 2:   ", seq2(1:nbp)
    else
      write(10,*) "STRAND 2:   "
    end if
    write(10,*) ""
!    write(10,*) "First column averages, second column standard deviations"

    l = 0
    do j=1,nbp-1

      write(10,*) j+1, "mer"
      write(10,trim(F_ELAP_1)) " base-step ", "Stretch", "Twist", "Roll", "Tilt", & 
                   & "Stretch-Twist", "Stretch-Roll", "Stretch-Tilt", &
                   & "Twist-Roll",  "Twist-Tilt", "Tilt-Roll", "Dynamic PL", &
                   & "Variance End-End", "pVariance End-End"
      write(10,*) "-----------------------------------------------------------------", &
                & "-----------------------------------------------------------------", &
                & "-----------------------------------------------------------------", &
                & "-----------------------------------------------------------------", &
                & "--------------" 

      !Write parameters
      if (str /= 2) then  !Linear or closed structure

        !Linear structure-------------------------------------------------------------------
        if (strands==2 ) then
          do i=1,nbp-j
            l=l+1
            write(10,trim(F_ELAP_2)) i, "-", i+j, " ", seq(i), seq(i+j), "/", seq(2*nbp-j-i+1), & 
                         & seq(2*nbp+1-i), F(:,l), Ad2(l), V(4,4,l), 1.0_dp/V__1(4,4,l) 
          end do !close i
        else 
          do i=1,nbp-j
            l=l+1
            write(10,trim(F_ELAP_2)) i, "-", i+j, " ", seq(i), seq(i+j), "/", "#", & 
                         & "#", F(:,l), Ad2(l), V(4,4,l), 1.0_dp/V__1(4,4,l) 
          end do !close i
        end if !Linear structure end ------------------------------------------------------

      else
        !closed structure-------------------------------------------------------------------
        if (strands==2 ) then !double stranded
          do i=1,nbp
            l=l+1

            !Check correct bp w (remember, we are on the other side when i > nbp-j)
            if (i > nbp-j) then
              w = i+j-nbp
            else
              w = i+j
            end if

            write(10,trim(F_ELAP_2)) i, "-", w, " ", seq(i), seq(w), "/", seq2(w), & 
                         & seq2(i), F(:,l), Ad2(l),  V(4,4,l), 1.0_dp/V__1(4,4,l) 
 
          end do !close i
 
        else !single stranded
          do i=1,nbp
            l=l+1

            !Check correct bp w (remember, we are on the other side when i > nbp-j)
            if (i > nbp-j) then
              w = i+j-nbp
            else
              w = i+j
            end if

            write(10,trim(F_ELAP_2)) i, "-", w, " ", seq(i), seq(w), "/", "#", & 
                         & "#", F(:,l), Ad2(l), V(4,4,l), 1.0_dp/V__1(4,4,l) 
 
          end do !close i

        end if ! strands

      end if !Closed structure end ------------------------------------------------------

      !Write overalls


      write(10,*) "-----------------------------------------------------------------", &
                & "-----------------------------------------------------------------", &
                & "-----------------------------------------------------------------", &
                & "-----------------------------------------------------------------", &
                & "--------------" 
     write(10,trim(F_ELAP_3)) " AVG-STD = ", OV_F(:,:,j), OV_Ad2(:,j), OV_V_E_E(:,j), OV_pV_E_E(:,j) 
     write(10,*) ""
    end do   !close j

  close(10, iostat=ierror)
  if (ierror/=0) stop "Error in closing BPP output file"
  
  end subroutine write_elastic_parms
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
  !Subroutine used by SerraNA.
  ! Reads directories of trajectory (traj) and topology (top) files.
  ! It returns an integer (str) which tells if the structure
  ! being analysed is double or single stranded, and a similar integer
  ! c_str, which indicates if is a circular or linear structure.
  ! All this information is indicated in the input file "s_NA.in
  subroutine SerraNA_inputs(traj,top,str,c_str)
  implicit none

  character(360), intent(out) :: top, traj
  integer, intent(out) :: str, c_str
  integer :: ierror

  !Closed structure?
  call nextpar(ierror)
  if (ierror /= 0) stop "Error, fill the parameters in input file (.in)"
  read(5,*) str

  !Circular structure?
  call nextpar(ierror)
  if (ierror /= 0) stop "Error, fill the parameters in input file (.in)"
  read(5,*) c_str

  !Topology
  call nextpar(ierror)
  if (ierror /= 0) stop "Error, fill the parameters in input file (.in)"
  read(5,"(A)") top
  top = adjustl(top)

  !Trajectory
  call nextpar(ierror)
  if (ierror /= 0) stop "Error, fill the parameters in input file (.in)"
  read(5,"(A)") traj
  traj = adjustl(traj)

  if (str == 1) then
    write(6,"(A)") "Single stranded structure"
  else if (str == 2) then
    write(6,"(A)") "Double stranded structure"
  else
    stop "Tell me if it is a single or double stranded structure (1 or 2)"
  end if

  if (c_str == 1) then
    write(6,"(A)") "Linear structure"
  else if (c_str == 2) then
    write(6,"(A)") "Circular structure"
  else
    stop "Tell me if it is a circular or linear structure (1 or 2)"
  end if

  write(6,"(2A)") "Topology file =", trim(top)
  write(6,"(2A)") "Trajectory file =", trim(traj)

  end subroutine SerraNA_inputs
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
  !Reads inputs of Analysis (for obtaining overall elastic constants)
  subroutine Analysis_inputs(r_A, r_T, r_S, elas_f, struc_f)
 
  implicit none

  !r_A : Ranges of persistence length dependant parameters
  !r_T : RAnges of twist
  !r_S : Ranges of stretch fitting
  integer, intent(out) :: r_A(2,2), r_T(2,2), r_S(2,2)
  !elas_f : direction of elastic_parameters file
  !struc_f : direction of structural parameters
  character(360), intent(out) :: elas_f, struc_f
  integer :: ierror

  !Ranges A
  call nextpar(ierror)
  if (ierror /= 0) stop "Error, fill the parameters in input file (.in)"
  read(5,*) r_A(:,:)

  !Ranges_T
  call nextpar(ierror)
  if (ierror /= 0) stop "Error, fill the parameters in input file (.in)"
  read(5,*) r_T(:,:)

  !Ranges_S
  call nextpar(ierror)
  if (ierror /= 0) stop "Error, fill the parameters in input file (.in)"
  read(5,*) r_S(:,:)

  !Elastic parameters
  call nextpar(ierror)
  if (ierror /= 0) stop "Error, fill the parameters in input file (.in)"
  read(5,"(A)") elas_f
  elas_f = adjustl(elas_f)

  !Structural parameters
  call nextpar(ierror)
  if (ierror /= 0) stop "Error, fill the parameters in input file (.in)"
  read(5,"(A)") struc_f
  struc_f = adjustl(struc_f)

  end subroutine Analysis_inputs
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
  !Reads inputs for Extract
  subroutine Extract_inputs(file_in, type_ext, sublength, a, b)
 
  implicit none

  !type_ext  : type of extraction ( = 0 or 1)
  !sublength : length to extract ( if type_ext = 0 )
  ![a,b]     : fragment subsection to obtain overalls
  integer, intent(out) :: type_ext, sublength, a, b
  character(360), intent(out) :: file_in
  integer :: ierror

  !Read path to input file
  call nextpar(ierror)
  if (ierror /= 0) stop "Error, fill the parameters in input file (.in)"
  read(5,"(A)",iostat=ierror) file_in
  if (ierror /= 0) stop "Error in reading path to input file"
  file_in = adjustl(file_in)

  !Read the type of extraction
  call nextpar(ierror)
  if (ierror /= 0) stop "Error, fill the parameters in input file (.in)"
  read(5,*,iostat=ierror) type_ext
  if (ierror /= 0) stop "Error in reading type of extraction"

  !Now, read eather the subplength or subfragment
  if (type_ext == 0) then

    !Read the type of extraction
    call nextpar(ierror)
    if (ierror /= 0) stop "Error, fill the parameters in input file (.in)"
    read(5,*,iostat=ierror) sublength
    if (ierror /= 0) stop "Error in reading sublength"

  else if (type_ext == 1) then

    !Read the type of extraction
    call nextpar(ierror)
    if (ierror /= 0) stop "Error, fill the parameters in input file (.in)"
    read(5,*,iostat=ierror) a, b
    if (ierror /= 0) stop "Error in reading subfragment [a,b]"

  else

    write(6,*) "Invalid extraction type, please write eather 0 or 1 in the second input"
    stop
  end if

  end subroutine Extract_inputs
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
  !Subroutine implemented in Extract.
  !Reads extract input file.
  !Depending if it is BSP, elasp, or strucp it allocates the data to
  !the corresponding variables
  subroutine read_extract_input_file(file_in, type_parm, nbp, frames, &
                        & strands, str, n_bsp, &
                        & seq_I, seq_II, BSP, OV_BSP, elasp, OV_elasp, & 
                        & strucp, avstrp, OV_strucp, OV_avstrp, BPP, OV_BPP)
  implicit none

  !file_in : parameters directory
  !nbp : number of base-pairs
  !frames : number of frames
  !strands : 2 if double stranded structure, 1 if single stranded
  !str : 2 if closed structure, 1 if opened
  !n_bsp: number of parameters
  !seq: sequence
  !OV_ : Overall parameters
  !BPP : base-pair parameters
  !BSP: base-step parameters
  !strucp: structural parameters
  !avstrp: average structure parameters
  !elasp: elastic parameters

  character(360), intent(in) :: file_in
  integer, intent(out) :: type_parm

  character(1), allocatable, intent(out) :: seq_I(:), seq_II(:)
  integer, intent(out) :: nbp, frames, strands, str, n_bsp
  real(dp), allocatable, intent(out) :: BSP(:,:,:), OV_BSP(:,:,:), &
                                      & elasp(:,:), OV_elasp(:,:,:), &
                                      & strucp(:,:,:), avstrp(:,:), &
                                      & OV_strucp(:,:,:), OV_avstrp(:,:,:), &
                                      & BPP(:,:,:), OV_BPP(:,:)

  integer :: ierror
  character(360) :: cdum

  open(unit=10, file=trim(file_in), action="read", status="old", iostat=ierror)
  if (ierror /= 0) stop "Error in opening input parameters file"

  !Let's only read the header, in order to identify the type of input
  read(10,"(A)") cdum

  close(10,iostat=ierror)
  if (ierror /= 0) stop "Error in closing input parameters file"


  !Let's check the type of file

  !If BSP
  if (cdum(2:21) .eq. "BASE-STEP PARAMETERS") then

    type_parm = 1
    write(6,*) "Reading base-step parameters"

    !And now the real reading
    call read_BSP(file_in, nbp, frames, strands, str, n_bsp, seq_I, seq_II, &
                      & BSP, OV_BSP)

  !If elastic
  else if (cdum(2:19) .eq. "ELASTIC PARAMETERS") then

    type_parm = 2
    write(6,*) "Reading elastic parameters"

    !And now the real reading
    call read_elastic_parms(file_in, nbp, frames, strands, str, n_bsp, seq_I, seq_II, &
                      & elasp, OV_elasp)

  !If structural
  else if (cdum(2:22) .eq. "STRUCTURAL PARAMETERS") then

    type_parm = 3
    write(6,*) "Reading structural parameters"

    call read_structural_parms(file_in, nbp, frames, strands, str, n_bsp, seq_I, seq_II, &
                               & strucp, avstrp, OV_strucp, OV_avstrp)

  !If BPP
  else if (cdum(2:21) .eq. "BASE-PAIR PARAMETERS") then

    type_parm = 4
    write(6,*) "Reading base-pair parameters"

    !And now the real reading
    call read_BPP(file_in, nbp, frames, strands, str, n_bsp, seq_I, seq_II, &
                      & BPP, OV_BPP)


  !If bad input file
  else 

    stop "Couldn't identify the type of data file, please check the header"

  end if

  end subroutine read_extract_input_file
!-----------------------------------------------------------------------
 

!-----------------------------------------------------------------------
  !Subroutine implemented in Extract.
  !Reads BSP output file

  subroutine read_BSP(BSP_file, nbp, frames, strands, str, n_bsp, &
                               & seq_I, seq_II, BSP, OV_BSP )
  implicit none

  !BSP_file : BSP_file
  !nbp : number of base-pairs
  !frames : number of frames
  !strands : 2 if double stranded structure, 1 if single stranded
  !str : 2 if closed structure, 1 if opened
  !n_bsp: number of parameters
  !seq: sequence
  !BSP: BSP 
  !OV_BSP : Overall BSP
  character(360), intent(in) :: BSP_file

  integer, intent(out) :: nbp, frames, strands, str, n_bsp
  character(1), allocatable, intent(out) :: seq_I(:), seq_II(:)
  real(dp), allocatable, intent(out) :: BSP(:,:,:), OV_BSP(:,:,:)

  integer :: ierror, i, j, l, idum
  character(360) :: cdum, FRMT

  open(unit=10, file=trim(BSP_file), action="read", status="old", iostat=ierror)
  if (ierror /= 0) stop "Error in opening BSP file"

  read(10,"(A)") cdum
  if (cdum(2:21) .ne. "BASE-STEP PARAMETERS") &
  & stop "Wrong BSP input file"
  !READ CLOSED OR OPEN
  read(10,*) cdum
  if (cdum(1:6) .eq. "CLOSED") then
    str = 2
  else
    str = 1
  end if
  !READ STR
  read(10,*) cdum
  if (cdum(1:6) .eq. "DOUBLE") then
    strands = 2
  else
    strands = 1
  end if
  !BP
  read(10,*) cdum, nbp
  !FRAMES
  read(10,*) cdum, frames
  read(10,*) cdum
  !STRAND 1
  allocate(seq_I(nbp), stat=ierror)
  if (ierror/=0) stop "Error in allocating sequence"
  if (nbp < 10) then
    write(FRMT,"(A6,I1,A3)") '(1A13,',nbp,'A1)' 
  else if (nbp < 100) then
    write(FRMT,"(A6,I2,A3)") '(1A13,',nbp,'A1)' 
  else if (nbp < 1000) then
    write(FRMT,"(A6,I3,A3)") '(1A13,',nbp,'A1)' 
  else
    write(FRMT,"(A6,I4,A3)") '(1A13,',nbp,'A1)' 
  end if
  read(10,trim(FRMT)) cdum, seq_I(:)
  !STRAND 2
  if (strands ==2) then
    allocate(seq_II(nbp), stat=ierror)
    if (ierror/=0) stop "Error in allocating sequence"
    read(10,trim(FRMT)) cdum, seq_II(:)
  else
    read(10,*) cdum
  end if
  read(10,*) cdum

  !Retrieve data
  if ( str == 2) then
    n_bsp = nbp! nbp*(nbp-1)   !Circular structure has more parameters
  else
    n_bsp = nbp-1!nbp*(nbp-1)/2
  end if
  allocate(BSP(2,7,n_bsp), OV_BSP(2,7,1), stat=ierror)
  if (ierror/=0) stop "Error in allocating BSP"

  !CLOSED STRUCTURE
  if (str ==2) then

    l = 0
    do j=1,1!nbp-1
      !read(10,*) cdum !bp 
      read(10,*) cdum !base-step, ...
      read(10,*) cdum !---------- ...

      do i=1,nbp
        l=l+1
        read(10,trim(F_BSP_2)) idum, cdum, idum, cdum, cdum, cdum, cdum, cdum, &
               & cdum, BSP(:,:,l)
      end do !close i
      read(10,*) cdum !---------- ...
      read(10,trim(F_BSP_3)) cdum, OV_BSP(:,:,j)
    end do !close j 

   !OPENED STRUCTURE
  else

    l = 0
    do j=1,1!nbp-1
!      read(10,*) cdum !bp 
      read(10,*) cdum !base-step, ...
      read(10,*) cdum !---------- ...
      do i=1,nbp-j
        l=l+1
        read(10,trim(F_BSP_2)) idum, cdum, idum, cdum, cdum, cdum, cdum, cdum, &
               & cdum, BSP(:,:,l)
      end do !close i
      read(10,*) cdum !---------- ...
      read(10,trim(F_BSP_3)) cdum, OV_BSP(:,:,j)
    end do !close j 
 
 end if 
 
  close(10) 
  if (ierror /= 0) stop "Error in closing BSP file"

  end subroutine read_BSP
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
  !Subroutine implemented in Analysis.
  !Reads elastic_parameters output file
  subroutine read_elastic_parms(elasparms, nbp, frames, strands, str, &
                               & n_bsp, seq_I, seq_II, elasp, OV_elasp)
  implicit none

  !elasparms : elastic parameters directory
  !nbp : number of base-pairs
  !frames : number of frames
  !strands : 2 if double stranded structure, 1 if single stranded
  !str : 2 if closed structure, 1 if opened
  !n_bsp: number of parameters
  !seq: sequence
  !elasp: elastic parameters
  !OV_elasp : Overall elastic parameters
  real(dp), allocatable :: elasp(:,:), OV_elasp(:,:,:)
  integer :: nbp, frames, strands, str, n_bsp, ierror, i, j, l, idum
  character(1), allocatable :: seq_I(:), seq_II(:)
  character(360) :: elasparms, cdum, FRMT

  open(unit=10, file=trim(elasparms), action="read", status="old", iostat=ierror)
  if (ierror /= 0) stop "Error in opening elastic parameters file"

  read(10,"(A)") cdum
  if (cdum(2:19) .ne. "ELASTIC PARAMETERS") &
  & stop "Wrong elastic parameters input file"
  !READ CLOSED OR OPEN
  read(10,*) cdum
  if (cdum(1:6) .eq. "CLOSED") then
    str = 2
  else
    str = 1
  end if
  !READ STR
  read(10,*) cdum
  if (cdum(1:6) .eq. "DOUBLE") then
    strands = 2
  else
    strands = 1
  end if
  !BP
  read(10,*) cdum, nbp
  !FRAMES
  read(10,*) cdum, frames
  read(10,*) cdum
  !STRAND 1
  allocate(seq_I(nbp), stat=ierror)
  if (ierror/=0) stop "Error in allocating sequence"
  if (nbp < 10) then
    write(FRMT,"(A6,I1,A3)") '(1A13,',nbp,'A1)' 
  else if (nbp < 100) then
    write(FRMT,"(A6,I2,A3)") '(1A13,',nbp,'A1)' 
  else if (nbp < 1000) then
    write(FRMT,"(A6,I3,A3)") '(1A13,',nbp,'A1)' 
  else
    write(FRMT,"(A6,I4,A3)") '(1A13,',nbp,'A1)' 
  end if
  read(10,trim(FRMT)) cdum, seq_I(:)
  !STRAND 2
  if (strands ==2) then
    allocate(seq_II(nbp), stat=ierror)
    if (ierror/=0) stop "Error in allocating sequence"
    read(10,trim(FRMT)) cdum, seq_II(:)
  else
    read(10,*) cdum
  end if
!  read(10,*) cdum !First column

  !Retrieve data
  if ( str == 2) then
    n_bsp = nbp*(nbp-1)   !Circular structure has more parameters
  else
    n_bsp = nbp*(nbp-1)/2
  end if
  allocate(elasp(13,n_bsp), OV_elasp(2,13,nbp-1), stat=ierror)
  if (ierror/=0) stop "Stop in allocating elastic parameters"

  !CLOSED STRUCTURE
  if (str ==2) then

    l = 0
    do j=1,nbp-1
      read(10,*) cdum !bp 
      read(10,*) cdum !base-step, ...
      read(10,*) cdum !---------- ...

      do i=1,nbp
        l=l+1

        read(10,trim(F_ELAP_2)) idum, cdum, idum, cdum, cdum, cdum, cdum, cdum, &
               & cdum, elasp(:,l)

      end do !close i

      read(10,*) cdum !---------- ...
      read(10,trim(F_ELAP_3)) cdum, OV_elasp(:,:,j)
 
    end do !close j

  !OPENED STRUCTURE
  else

    l = 0
    do j=1,nbp-1
      read(10,*) cdum !bp 
      read(10,*) cdum !base-step, ...
      read(10,*) cdum !---------- ...
      do i=1,nbp-j
        l=l+1
        read(10,trim(F_ELAP_2)) idum, cdum, idum, cdum, cdum, cdum, cdum, cdum, &
               & cdum, elasp(:,l)
      end do !close i
      read(10,*) cdum !---------- ...
      read(10,trim(F_ELAP_3)) cdum, OV_elasp(:,:,j)
    end do !close j 

  end if 
 
  close(10) 
  if (ierror /= 0) stop "Error in closing elastic parameters file"

  end subroutine read_elastic_parms
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
  !Subroutine implemented in Analysis.
  !Reads structural_parameters output file

  subroutine read_structural_parms(strucparms, nbp, frames, strands, str, n_bsp, &
                               & seq_I, seq_II, strucp, avstrp, &
                               & OV_strucp, OV_avstrp)
  implicit none

  !strucparms : sctructural parameters directory
  !nbp : number of base-pairs
  !frames : number of frames
  !strands : 2 if double stranded structure, 1 if single stranded
  !str : 2 if closed structure, 1 if opened
  !n_bsp: number of parameters
  !seq: sequence
  !elasp: elastic parameters
  !OV_elasp : Overall elastic parameters 
  real(dp), allocatable :: strucp(:,:,:), avstrp(:,:), &
                         & OV_avstrp(:,:,:), OV_strucp(:,:,:)
  integer :: nbp, frames, strands, str, n_bsp, ierror, i, j, l, idum
  character(1), allocatable :: seq_I(:), seq_II(:)
  character(360) :: strucparms, cdum, FRMT

  open(unit=10, file=trim(strucparms), action="read", status="old", iostat=ierror)
  if (ierror /= 0) stop "Error in opening elastic parameters file"

  read(10,"(A)") cdum
  if (cdum(2:22) .ne. "STRUCTURAL PARAMETERS") &
  & stop "Wrong structural parameters input file"
  !READ CLOSED OR OPEN
  read(10,*) cdum
  if (cdum(1:6) .eq. "CLOSED") then
    str = 2
  else
    str = 1
  end if
  !READ STR
  read(10,*) cdum
  if (cdum(1:6) .eq. "DOUBLE") then
    strands = 2
  else
    strands = 1
  end if
  !BP
  read(10,*) cdum, nbp
  !FRAMES
  read(10,*) cdum, frames
  read(10,*) cdum
  !STRAND 1
  allocate(seq_I(nbp), stat=ierror)
  if (ierror/=0) stop "Error in allocating sequence"
  if (nbp < 10) then
    write(FRMT,"(A6,I1,A3)") '(1A13,',nbp,'A1)' 
  else if (nbp < 100) then
    write(FRMT,"(A6,I2,A3)") '(1A13,',nbp,'A1)' 
  else if (nbp < 1000) then
    write(FRMT,"(A6,I3,A3)") '(1A13,',nbp,'A1)' 
  else
    write(FRMT,"(A6,I4,A3)") '(1A13,',nbp,'A1)' 
  end if
  read(10,trim(FRMT)) cdum, seq_I(:)
  !STRAND 2
  if (strands ==2) then
    allocate(seq_II(nbp), stat=ierror)
    if (ierror/=0) stop "Error in allocating sequence"
    read(10,trim(FRMT)) cdum, seq_II(:)
  else
    read(10,*) cdum
  end if
  read(10,*) cdum

  !Retrieve data
  if ( str == 2) then
    n_bsp = nbp*(nbp-1)   !Circular structure has more parameters
  else
    n_bsp = nbp*(nbp-1)/2
  end if
  allocate(strucp(2,11,n_bsp), avstrp(3,n_bsp), &
         & OV_strucp(2,11,nbp-1), OV_avstrp(2,3,nbp-1), stat=ierror)
  if (ierror/=0) stop "Stop in allocating structural parameters"

  !CLOSED STRUCTURE
  if (str ==2) then

    l = 0
    do j=1,nbp-1
      read(10,*) cdum !bp 
      read(10,*) cdum !base-step, ...
      read(10,*) cdum !---------- ...

      do i=1,nbp
        l=l+1
        read(10,trim(F_STRP_2)) idum, cdum, idum, cdum, cdum, cdum, cdum, cdum, &
               & cdum, strucp(:,:,l), avstrp(:,l)
      end do !close i
      read(10,*) cdum !---------- ...
      read(10,trim(F_STRP_3)) cdum, OV_strucp(:,:,j), OV_avstrp(:,:,j)
    end do !close j 

   !OPENED STRUCTURE
  else

    l = 0
    do j=1,nbp-1
      read(10,*) cdum !bp 
      read(10,*) cdum !base-step, ...
      read(10,*) cdum !---------- ...
      do i=1,nbp-j
        l=l+1
        read(10,trim(F_STRP_2)) idum, cdum, idum, cdum, cdum, cdum, cdum, cdum, &
               & cdum, strucp(:,:,l), avstrp(:,l)
      end do !close i
      read(10,*) cdum !---------- ...
      read(10,trim(F_STRP_3)) cdum, OV_strucp(:,:,j), OV_avstrp(:,:,j)
    end do !close j 
 
 end if 
 
  close(10) 
  if (ierror /= 0) stop "Error in closing structural parameters file"

  end subroutine read_structural_parms
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
  !Subroutine implemented in Extract.
  !Reads BPP output file

  subroutine read_BPP(BPP_file, nbp, frames, strands, str, n_bsp, &
                               & seq_I, seq_II, BPP, OV_BPP )
  implicit none

  !BPP_file : BPP_file
  !nbp : number of base-pairs
  !frames : number of frames
  !strands : 2 if double stranded structure, 1 if single stranded
  !str : 2 if closed structure, 1 if opened
  !n_bsp: number of parameters
  !seq: sequence
  !BPP: BPP 
  character(360), intent(in) :: BPP_file

  integer, intent(out) :: nbp, frames, strands, str, n_bsp
  character(1), allocatable, intent(out) :: seq_I(:), seq_II(:)
  real(dp), allocatable, intent(out) :: BPP(:,:,:), OV_BPP(:,:)

  integer :: ierror, i, idum
  character(360) :: cdum, FRMT

  open(unit=10, file=trim(BPP_file), action="read", status="old", iostat=ierror)
  if (ierror /= 0) stop "Error in opening BSP file"

  read(10,"(A)") cdum
  if (cdum(2:21) .ne. "BASE-PAIR PARAMETERS") &
  & stop "Wrong BPP input file"
  !READ CLOSED OR OPEN
  read(10,*) cdum
  if (cdum(1:6) .eq. "CLOSED") then
    str = 2
  else
    str = 1
  end if
  !READ STR
  read(10,*) cdum
  if (cdum(1:6) .eq. "DOUBLE") then
    strands = 2
  else
    strands = 1
  end if
  !BP
  read(10,*) cdum, nbp
  !FRAMES
  read(10,*) cdum, frames
  read(10,*) cdum
  !STRAND 1
  allocate(seq_I(nbp), stat=ierror)
  if (ierror/=0) stop "Error in allocating sequence"
  if (nbp < 10) then
    write(FRMT,"(A6,I1,A3)") '(1A13,',nbp,'A1)' 
  else if (nbp < 100) then
    write(FRMT,"(A6,I2,A3)") '(1A13,',nbp,'A1)' 
  else if (nbp < 1000) then
    write(FRMT,"(A6,I3,A3)") '(1A13,',nbp,'A1)' 
  else
    write(FRMT,"(A6,I4,A3)") '(1A13,',nbp,'A1)' 
  end if
  read(10,trim(FRMT)) cdum, seq_I(:)
  !STRAND 2
  if (strands ==2) then
    allocate(seq_II(nbp), stat=ierror)
    if (ierror/=0) stop "Error in allocating sequence"
    read(10,trim(FRMT)) cdum, seq_II(:)
  else
    read(10,*) cdum
  end if
  read(10,*) cdum

  !Retrieve data
  n_bsp = nbp   !Does not matter if circular or linear
  allocate(BPP(2,6,n_bsp), OV_BPP(2,6), stat=ierror)
  if (ierror/=0) stop "Error in allocating BPP"

  read(10,*) cdum !base-pair, ...
  read(10,*) cdum !---------- ...
  !Again, does not matter if opened or closed
  do i = 1, nbp
    read(10,trim(F_BPP_2)) idum, cdum, cdum, cdum, cdum, &
           & BPP(:,:,i)
  end do !close i

  read(10,*) cdum !---------- ...
  read(10,trim(F_BPP_3)) cdum, OV_BPP(:,:)
 
  close(10) 
  if (ierror /= 0) stop "Error in closing BSP file"

  end subroutine read_BPP
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
! This subroutine prints on screen information of elastic constants
! calculated by Analysis.f90
  subroutine print_elastic_constants(Tilt, Roll, Twist, Stretch, A_a, &
                                   & As_a, Ad_a, A_bb, Ad_c, A_d, &
                                   & r_A, r_S, r_T, str, frames, nbp )  

  !Elastic constants
  real(dp), intent(in) :: Tilt(:), Roll(:), Twist(:), Stretch(:), &
                        & A_a(:), As_a(:), Ad_a(:), A_bb(:), Ad_c(:), &
                        & A_d(:)

  !Ranges
  integer, intent(in) :: r_A(:,:), r_S(:,:), r_T(:,:), str, frames, nbp

  if ( str == 2 ) then
    write(6,*) "CLOSED STRUCTURE ANALYSED"
  else
    write(6,*) "LINEAR STRUCTURE ANALYSED"
  end if

  write(6,*) "BASE-PAIRS: ", nbp
  write(6,*) "FRAMES:     ", frames


  write(6,trim(F_OV_1)) "","Elastic cte","Intercept","Slope","Confidence-I","Strd Error"
  write(6,*) '-------------------------------------------------------------------------------------------'
  write(6,trim(F_OV_2)) "Tilt (nm): ", Tilt(1),"#","#","#",Tilt(2)
  write(6,trim(F_OV_2)) "Roll (nm): ", Roll(1),"#","#","#",Roll(2)
  write(6,trim(F_OV_2)) "Twist (nm): ", Twist(1),"#","#","#",Twist(2)
  write(6,trim(F_OV_3)) "Stretch (pN) : ", Stretch(:),"#"
  write(6,trim(F_OV_3)) "A [a] (nm):", A_a(:),"#"
  write(6,trim(F_OV_3)) "As [a] (nm):", As_a(:),"#"
  write(6,trim(F_OV_3)) "Ad [a] (nm):", Ad_a(:),"#"
  write(6,trim(F_OV_4)) "A [b] (nm): ", A_bb(1),"#","#","#","#"
  write(6,trim(F_OV_2)) "Ad [c] (nm): ", Ad_c(1),"#","#","#",Ad_c(2)
  write(6,trim(F_OV_4)) "A [d] (nm): ", A_d(1),"#","#","#","#"

  write(6,*) ''
  if (r_A(1,1) /= r_A(2,1) ) then
    write(6,trim(F_OV_5)) 'Persistence lengths, Tilt and Roll calculated from base: ', r_A(1,1), &
    " to base", r_A(2,1), " and from lengths: ", r_A(1,2), " to ", r_A(2,2)
  else !IF DEFAULT CASE
    write(6,trim(F_OV_6)) "Persistence lengths, Tilt and Roll calculated from whole fragment "//&
                       &  "and from lengths: ", r_A(1,2), " to ", r_A(2,2)
  end if

  if (r_T(1,1) /= r_T(2,1) ) then
    write(6,trim(F_OV_5)) 'Twist calculated from base: ', r_T(1,1), " to base ", r_T(2,1), &
    " and from lengths: ", r_T(1,2), " to ", r_T(2,2)
  else !IF DEFAULT CASE
    write(6,trim(F_OV_6)) "Twist calculated from whole fragment "//&
                       &  "and from lengths: ", r_T(1,2), " to ", r_T(2,2)
  end if

  if (r_S(1,1) /= r_S(2,1) ) then
    write(6,trim(F_OV_5)) 'Stretch calculated from base: ', r_S(1,1), " to base", r_S(2,1), &
    " and from lengths: ", r_S(1,2), " to ", r_S(2,2)
  else !IF DEFAULT CASE
    write(6,trim(F_OV_6)) "Stretch calculated from whole fragment "//&
                       &  "and from lengths: ", r_S(1,2), " to ", r_S(2,2)
  end if

  write(6,*) ''
  write(6,*) 'A = Persistence lentgh, As = Static persistence length, Ad = Dynamic persistence length'
  write(6,*) '(a) => overall fitting'
  write(6,*) '(b) => A = As*Ad/(As+Ad)'
  write(6,*) '(c) => Ad = Tilt*Roll/(Tilt+Roll)'
  write(6,*) '(d) => A = As*Ad/(As+Ad) with Ad obtained by (c) '

  end subroutine print_elastic_constants
!----------------------------------------------------------------------


!-----------------------------------------------------------------------
  !Reads chr, if first character is different from 'c' (comment line), then 
  !it exits and then next line should be readed
  subroutine nextpar(ierror)
  implicit none
    character(len=1) :: chr
    integer, intent(out) :: ierror
    do 
      read(5,'(A1)', iostat=ierror) chr
      if (chr(1:1) /= 'c') exit
      if ( ierror > 0) stop "Error in reading input file"
      if ( ierror < 0) exit
    end do

  end subroutine nextpar
!----------------------------------------------------------------------

  end module IO_mod
