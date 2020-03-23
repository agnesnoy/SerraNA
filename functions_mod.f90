! Module containing all the mathematical functions needed for      31/07/2019
! SerraNA 
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

 module functions_mod

 use parms

 implicit none

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

 contains


!SerraNA FUNCTIONS
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
 !Subroutine used for SerraNA
 !Calculates the orientation matrix "R" and origin vector "O" of a base.
 !"Ex" is the coordinates of the experimental base and "St" are the 
 !coordinates of the standard base. These bases have N atoms
 !N__1 = 1/N, N__1_1 = 1/(N-1)
 subroutine Get_Rotation_R_Origin_O(R, O, Ex, St, N, N__1, N__1_1)
 integer, intent(in) :: N
 real(dp), intent(in) :: Ex(3,N), St(3,N), N__1, N__1_1
 real(dp), intent(out) :: R(3,3), O(3)
 !C covariance matrix
 !I array with entries = 1
 !A, B  auxiliar matrices
 !M real symmetrix matrix
 !S diagonalization of matrix M
 !V eigenvectors of matrix M
 real(dp) :: av_Ex(3), av_St(3), C(3,3), &
           & M(4,4), eigenvec(4), I(N,1), & 
           & S(4,4), V(4,4), A(3,3), B(3,3)
 integer :: k, j

 I = one
 !Calculates averages
 av_Ex = 0.0_dp
 av_St = 0.0_dp
 do k=1,N
   do j=1,3
     av_Ex(j)=av_Ex(j)+Ex(j,k)
     av_St(j)=av_St(j)+St(j,k)
   end do
 end do
 av_Ex=av_Ex/real(N,dp)
 av_St=av_St/real(N,dp)

 !Calculates covariance matrix C
 !** C = [1/(N-1)][S'E-(1/N)S'II'E] this is the original formula
 !** St = S', Ex = E' because of order that the matrices are stored
 !** I = I, N = N
 A = MATMUL(MATMUL(St,I),MATMUL(transpose(I),transpose(Ex))) ![S'II'E]
 B = MATMUL(St,transpose(Ex)) !S'E
 C = N__1_1*( B - N__1*A )

 !Calculates Real Symmetric Matrix M
 call RSmatrixM(C,M)

 !Get eigenvectors of matrix M. 4 dimension of M
 call diagonalization_Jacobi(M,S,V,4) 

 !And the largest eigen vector of V
 call largest_eigenvec(S,V,eigenvec,4)

 !Now to calculate the rotation matrix
 call get_rotation_R(R,eigenvec)

 !And origin O
 O=av_Ex-MATMUL(av_St,transpose(R))

 end subroutine Get_Rotation_R_Origin_O
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! Calculates the 4x4 real symmetric matrix M from 3x3 matrix C
 subroutine RSmatrixM(C,M)
 real(dp), intent(in) :: C(3,3)
 real(dp), intent(out) :: M(:,:)
 integer i,j

 M(1,1)=C(1,1)+C(2,2)+C(3,3)
 M(1,2)=C(2,3)-C(3,2)
 M(1,3)=C(3,1)-C(1,3)
 M(1,4)=C(1,2)-C(2,1)
 M(2,2)=C(1,1)-C(2,2)-C(3,3)
 M(2,3)=C(1,2)+C(2,1)
 M(2,4)=C(3,1)+C(1,3)
 M(3,3)=-C(1,1)+C(2,2)-C(3,3)
 M(3,4)=C(2,3)+C(3,2)
 M(4,4)=-C(1,1)-C(2,2)+C(3,3)
 
 do i=1,4
  do j=1,i
   M(i,j)=M(j,i)
  end do
 end do

 end subroutine RSmatrixM

!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! Only works with real symmetric matrix "M" (ndim x ndim)
! This one returns the Diagonalization matrix S with eigenvalues in the diagonal
! and the eigenvectors V. Eigenvalue S(i,i) corresponds to eigenvector V(:,i)
! It uses the Jacobi diagonalization method

 subroutine diagonalization_Jacobi(M,S,V,ndim)
 integer, intent(in) :: ndim
 real(dp), intent(in) :: M(ndim,ndim)
 real(dp), intent(out) :: S(ndim,ndim),V(ndim,ndim)
 integer :: i,j,l,k, r, cont
 real(dp) :: larg, th, n, cond, co, si,  &
             upS(ndim,ndim), Image(ndim,ndim) ,&
             G(ndim,ndim)
 
 n = real(ndim,dp)
 S = M
 cont = 50
 Image = 0.0_dp
 do k=1,ndim
   Image(k,k) = 1.0_dp
 end do

 V = Image
 cond = 4 !Just because
 r = 0
 do
   if (cond <= eps .or. r==cont) then
     S = upS
     exit
   end if
   if (r/=0) then
     S=upS
   end if
   r=r+1
   larg=-1.0_dp
   do k=1,ndim-1
     do l=k+1,ndim
       if ( larg < abs( S(k,l) ) ) then
         larg = abs( S(k,l) )
         i = k
         j = l
       end if
     end do
   end do
   if ( S(i,i) == S(j,j) ) then
     th= pi/4.0_dp
   else
     th=0.5_dp*datan( 2.0_dp*S(i,j)/( S(j,j) - S(i,i) ) )
   end if
   co = dcos(th)
   si = dsin(th)
   G = Image
   G(i,i) = co
   G(j,j) = co
   G(i,j) = si
   G(j,i) = -si
   V = MATMUL(V,G)
   upS = MATMUL(transpose(G),S)
   upS = MATMUL(upS,G)
   cond = 0.0_dp
   do k=1,ndim-1
     do l=k+1,ndim
       cond = cond + abs( upS(k,l) )
     end do
   end do
   cond=cond*2.0_dp
 end do

 end subroutine diagonalization_Jacobi
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! S is a diagonal matrix obtained from M. The diagonal values are the
! eigenvalues of M, and matrix V contains the eigenvectors
! The eigenvector of eigenvalue S(i,i), is column V(:,i)
! Now, to find the largest eigenvalue and its corresponding eigenvector
 subroutine largest_eigenvec(S,V,eigenvec,n)
 integer, intent(in) :: n
 real(dp), intent(in) :: S(n,n),V(n,n)
 real(dp), intent(out) :: eigenvec(n)
 real(dp) :: larg
 integer :: i,j,k
 
 larg = eps
 k = 0

 do i=1,n
   if( S(i,i) < 0.0_dp) then
     k=k+1
   end if
   if (k == n) stop "Error, all eigenvalues negative? Please check"
   if(larg < S(i,i)) then
     larg = S(i,i)
     j=i
   end if
 end do

 !Eigenvector with largest eigenvalue:
 eigenvec=V(:,j)
 end subroutine largest_eigenvec
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
 !Calculates rotation matrix R from eigenvector q
 subroutine get_rotation_R(R,q)
 real(dp), intent(out) :: R(3,3)
 real(dp), intent(in) :: q(4)

 R(1,1)=q(1)*q(1)+q(2)*q(2)-q(3)*q(3)-q(4)*q(4)
 R(1,2)=2*(q(2)*q(3)-q(1)*q(4))
 R(1,3)=2*(q(2)*q(4)+q(1)*q(3))
 R(2,1)=2*(q(3)*q(2)+q(1)*q(4))
 R(2,2)=q(1)*q(1)-q(2)*q(2)+q(3)*q(3)-q(4)*q(4)
 R(2,3)=2*(q(3)*q(4)-q(1)*q(2))
 R(3,1)=2*(q(4)*q(2)-q(1)*q(3))
 R(3,2)=2*(q(4)*q(3)+q(1)*q(2))
 R(3,3)=q(1)*q(1)-q(2)*q(2)-q(3)*q(3)+q(4)*q(4)
 
 end subroutine get_rotation_R
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
 !Calculates base-pair parameters and store them in BPP
 subroutine basepair_parameters(O1,O2,R1,R2,BPP,Tmbt,Ombt)
 real(dp), intent(out) :: BPP(6),Tmbt(3,3),Ombt(3)
 real(dp), intent(in) :: O1(3),O2(3),R1(3,3),R2(3,3)
 real(dp) :: bp(3), Rbp(3,3), T1(3,3), T2(3,3), h(3), &
             delta, psi
 integer :: i

 !Tmbt-> mid-base-triad
 !Ombt -> position of Tmbt
 !O1, R1 -> position and orientation of base 1
 !For BPP(i)
 !i=1, Shear (Sx)           
 !i=2, Streatch (Sy)       
 !i=3, Stagger (Sz)       
 !i=4, Buckle (kappa)      
 !i=5, Propeller (omega)  
 !i=6, Opening (sigma)    

 !BucklePropeller angle delta = arccos(z1.z2)

 delta=dacos(dot_product(R2(:,3),R1(:,3)))
 
 !BucklePropeller axis bp
 bp = cross_product3(R2(:,3),R1(:,3))
 bp = normalize_vector(bp,3)
 
 !Rotate triads T1, T2 about -+1/2 to bp.
 !Rbp auxiliar rotation matrix
 call general_rotation_matrix(Rbp,bp,-0.5_dp*delta)
 T1= MATMUL(Rbp,R1)

 call general_rotation_matrix(Rbp,bp,0.5_dp*delta)
 T2= MATMUL(Rbp,R2)

 !Construct Tmbt. Is obtained by averaging T1 and T2, then normalazing
 !each column vector

 Tmbt=(T1+T2)
 Tmbt=0.5_dp*Tmbt
 do i=1,3
  Tmbt(:,i) = normalize_vector(Tmbt(:,i), 3)
 end do

 !Origin Ombt is the average position between O1+O2
 Ombt=(O1+O2)
 Ombt=0.5_dp*Ombt

 !Get opening BPP(6) (sigma)

 h = cross_product3(T2(:,2),T1(:,2)) !h auxiliar vector

 !------ !!!!!This is stupid but happens!!!!!
 !So just to be sure
 if (dot_product(T2(:,2),T1(:,2)) > 1) then
   BPP(6)=0.0_dp
 else if (dot_product(T2(:,2),T1(:,2)) < -1) then
   BPP(6)=pi
 else
   BPP(6)=dacos(dot_product(T2(:,2),T1(:,2)))
 end if

 !------ Getting correct sign
 if (dot_product(h,Tmbt(:,3)) < 0.0_dp) then
  BPP(6) = -BPP(6)
 end if
 !--------------------------------------------------
 

 !The angle between the Buckle-Opening axis and the 
 !MBT x-axis is psi
 h = cross_product3(bp,Tmbt(:,2))  !h is auxiliar

 if (dot_product(bp,Tmbt(:,2)) > 1) then
  psi=0.0_dp
 else if (dot_product(bp,Tmbt(:,2)) < -1) then
   psi = pi
 else
   psi = dacos(dot_product(bp,Tmbt(:,2)))
 end if
 !------ Getting correct sign
 if (dot_product(h,Tmbt(:,3)) < 0.0_dp) then
   psi = -psi
 end if

 !Propeller BPP(5)
 BPP(5) = delta*dcos(psi)

 !Buckle BPP(4)    
 BPP(4) = delta*dsin(psi)

 !Displacement parameters
 !BPP(1), Shear (Sx)           
 !BPP(2), Streatch (Sy)       
 !BPP(3), Stagger (Sz)    

 BPP(1:3)=MATMUL(O1-O2,Tmbt)

 !Convert rads to degrees for angular variables
 do i=4,6
  BPP(i)=rad_to_deg*BPP(i)
 end do

 end subroutine basepair_parameters
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
 !Calculates base-step_parameters and store them in BPP
 subroutine basestep_parameters(O1, O2, R1, R2, BSP, totaltwist)
 real(dp), intent(out) :: BSP(:)
 real(dp), intent(in) :: O1(:),O2(:),R1(:,:),R2(:,:), totaltwist
 real(dp) :: rt(3), Rrt(3,3), T1(3,3), T2(3,3), Omst(3), Tmst(3,3), &
          &  h(3), phi
 integer :: i, N
 !BSP(j)
 !Total twist is the sum of previous twists
 !j=1, Shift (Dx)
 !j=2, Slide (Dy)
 !j=3, Rise (Dz)
 !j=4, Tilt (tao)
 !j=5, Roll (rho)
 !j=6, Twist (Omega)
 !j=7, Bending (theta)
 !j=8, Directional decay
 !j=9, Squared bending (theta**2)

 !Directional decay
 BSP(8) = dot_product(R1(:,3),R2(:,3))
 
 !Bending, which is the angle RoltTilt or between z1 and z2
 BSP(7)=dacos(BSP(8))

 !RoltTilt axis rt
 rt = cross_product3(R1(:,3),R2(:,3))
 rt = normalize_vector(rt,3)

 !Rotation matrices Rrt and rotate orientation matrices R1,R2 => T1,T2

 call general_rotation_matrix(Rrt, rt, 0.5_dp*BSP(7))
 T1 = MATMUL(Rrt, R1)

 call general_rotation_matrix(Rrt, rt, -0.5_dp*BSP(7))
 T2 = MATMUL(Rrt, R2)

 !Construct Tmst. Is obtained by averaging T1 and T2, then normalazing
 !each column vector

 Tmst = (T1+T2)
 Tmst = 0.5_dp*Tmst

 do i=1,3
   Tmst(:,i) = normalize_vector(Tmst(:,i),3)
 end do

 !Calculates Twist, BSP(6)

 h = cross_product3(T1(:,2),T2(:,2))
 !------ !!!!!This is stupid but happens!!!!!
 !So just to be sure
 if (dot_product(T1(:,2),T2(:,2)) > one) then
   BSP(6)=0.0_dp
 else if (dot_product(T1(:,2),T2(:,2)) < -one) then
   BSP(6)=pi
 else
   BSP(6)=dacos(dot_product(T1(:,2),T2(:,2)))
 end if
 !------ Getting correct sign
 if (dot_product(h,Tmst(:,3)) < 0.0_dp) then
   BSP(6) = -BSP(6)
 end if

 !-------------------------------------------
 !Obtain the correct twist
 
 N = int(totaltwist/180.0_dp) !Number of half turns

 if (mod(N,2) == 1) then
   BSP(6) = BSP(6)+2.0_dp*pi
 end if

! Origin Ombt is average position
 Omst = (O1+O2)*0.5_dp

 N = int(totaltwist/360.0_dp) !Number of turns
 BSP(6) = BSP(6) + N*2.0_dp*pi

 !Due to problems near pi
 !Note that totaltwist is in degrees and BSP in radians
 !The resulting BSP must be close to totaltwist.
 if (totaltwist-rad_to_deg*BSP(6) >  180.0_dp) BSP(6) = BSP(6)+2.0_dp*pi
 if (totaltwist-rad_to_deg*BSP(6) < -180.0_dp) BSP(6) = BSP(6)-2.0_dp*pi

 !Let's try and correct the x and y axis.
 N = int( ( BSP(6)+pi )/( 2.0_dp*pi ) )

 Tmst(:,1:2)=Tmst(:,1:2)*(-1)**N
 !-------------------------------------------

 !The angle between the Roll-Tilt axis and the 
 !MST y-axis is phi:
 h = cross_product3(rt,Tmst(:,2))
 !------ !!!!!This is stupid but happens!!!!!
 !So just to be sure
 if (dot_product(rt,Tmst(:,2)) > one) then
   phi = 0.0_dp
 else if (dot_product(rt,Tmst(:,2)) < -one) then
   phi = pi
 else
   phi = dacos(dot_product(rt,Tmst(:,2)))
 end if
 !------ Getting correct sign
 if (dot_product(h,Tmst(:,3)) < 0.0_dp) then
   phi = -phi
 end if

 !Tilt BSP(4)
 BSP(4)=BSP(7)*dsin(phi)

 !Roll BSP(5)
 BSP(5)=BSP(7)*dcos(phi)

 !Displacement parameters
 !BSP(1), Shift (Dx)
 !BSP(2), Slide (Dy)
 !BSP(3), Rise (Dz) 
 BSP(1:3) = MATMUL(O2-O1,Tmst)


! Let's convert rads to degrees for angular variables
 do i=4,7!8 !!!!!!NOTE THAT NOW IS 7 BECAUSE COMPONENT 8 IS NOT AN ANGLE 
  BSP(i)=rad_to_deg*BSP(i)
 end do

 BSP(9)=BSP(7)*BSP(7)  !This is squared bending. It will be used in
                       !calculating Persistence length

 end subroutine basestep_parameters
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
 !Calculates covariance matrix V of deformations variables:
 !Roll,Tilt,Twist,Streatch
 !Roll = BSP(5), Tilt = BSP(4), Twist = BSP(6), Streatch=BPP(2)
 !First Covariance matrix V
 subroutine deformation_covariance(V, Roll, Tilt, Twist, Streatch, n)
 real(dp), intent(in) :: Roll(:), Tilt(:), Twist(:), Streatch(:)
 real(dp), intent(out) :: V(:,:)
 real(dp) :: av(4),X(n,4),Y(n,4)
 integer :: i, j, n

 !av is the average vector, which contains the averages of the deformation variables
 X(:,1)=Roll
 X(:,2)=Tilt
 X(:,3)=Twist
 X(:,4)=Streatch
 av=0.0_dp
 av=sum(X,1)
 av=av/real(n,dp)
 do j=1,4
   Y(:,j)=X(:,j)-av(j)
 end do

 do i=1,4
   do j=1,4
     V(i,j)= sum(Y(:,i)*Y(:,j))
   end do
 end do
 V=V/real(n,dp)

 end subroutine deformation_covariance
!-----------------------------------------------------------------

 subroutine inverse_matrix_analytic4x4(V,InverV)
 real(dp), intent(in) :: V(4,4)
 real(dp), intent(out) :: InverV(4,4)
 real(dp) :: b(4,4), det
 integer i

 det= V(1,4)*V(2,3)*V(3,2)*V(4,1) - V(1,3)*V(2,4)*V(3,2)*V(4,1) - &
 V(1,4)*V(2,2)*V(3,3)*V(4,1) + V(1,2)*V(2,4)*V(3,3)*V(4,1)+ &
 V(1,3)*V(2,2)*V(3,4)*V(4,1) - V(1,2)*V(2,3)*V(3,4)*V(4,1) - &
 V(1,4)*V(2,3)*V(3,1)*V(4,2) + V(1,3)*V(2,4)*V(3,1)*V(4,2)+ &
 V(1,4)*V(2,1)*V(3,3)*V(4,2) - V(1,1)*V(2,4)*V(3,3)*V(4,2) - &
 V(1,3)*V(2,1)*V(3,4)*V(4,2) + V(1,1)*V(2,3)*V(3,4)*V(4,2)+ &
 V(1,4)*V(2,2)*V(3,1)*V(4,3) - V(1,2)*V(2,4)*V(3,1)*V(4,3) - &
 V(1,4)*V(2,1)*V(3,2)*V(4,3) + V(1,1)*V(2,4)*V(3,2)*V(4,3)+ &
 V(1,2)*V(2,1)*V(3,4)*V(4,3) - V(1,1)*V(2,2)*V(3,4)*V(4,3) - &
 V(1,3)*V(2,2)*V(3,1)*V(4,4) + V(1,2)*V(2,3)*V(3,1)*V(4,4)+ &
 V(1,3)*V(2,1)*V(3,2)*V(4,4) - V(1,1)*V(2,3)*V(3,2)*V(4,4) - &
 V(1,2)*V(2,1)*V(3,3)*V(4,4) + V(1,1)*V(2,2)*V(3,3)*V(4,4)


 b(1,1)=V(2,3)*V(3,4)*V(4,2)-V(2,4)*V(3,3)*V(4,2)+ &
 V(2,4)*V(3,2)*V(4,3)-V(2,2)*V(3,4)*V(4,3)-V(2,3)*V(3,2)*V(4,4)+ &
 V(2,2)*V(3,3)*V(4,4)

 b(1,2)=V(1,4)*V(3,3)*V(4,2)-V(1,3)*V(3,4)*V(4,2)- &
 V(1,4)*V(3,2)*V(4,3)+V(1,2)*V(3,4)*V(4,3)+V(1,3)*V(3,2)*V(4,4)- &
 V(1,2)*V(3,3)*V(4,4)

 b(1,3)=V(1,3)*V(2,4)*V(4,2)-V(1,4)*V(2,3)*V(4,2)+ &
 V(1,4)*V(2,2)*V(4,3)-V(1,2)*V(2,4)*V(4,3)-V(1,3)*V(2,2)*V(4,4)+ &
 V(1,2)*V(2,3)*V(4,4)

 b(1,4)=V(1,4)*V(2,3)*V(3,2)-V(1,3)*V(2,4)*V(3,2)- &
 V(1,4)*V(2,2)*V(3,3)+V(1,2)*V(2,4)*V(3,3)+V(1,3)*V(2,2)*V(3,4)- &
 V(1,2)*V(2,3)*V(3,4)

 b(2,1)=V(2,4)*V(3,3)*V(4,1)-V(2,3)*V(3,4)*V(4,1)- &
 V(2,4)*V(3,1)*V(4,3)+V(2,1)*V(3,4)*V(4,3)+V(2,3)*V(3,1)*V(4,4)- &
 V(2,1)*V(3,3)*V(4,4)

 b(2,2)=V(1,3)*V(3,4)*V(4,1)-V(1,4)*V(3,3)*V(4,1)+ &
 V(1,4)*V(3,1)*V(4,3)-V(1,1)*V(3,4)*V(4,3)-V(1,3)*V(3,1)*V(4,4)+ &
 V(1,1)*V(3,3)*V(4,4)

 b(2,3)=V(1,4)*V(2,3)*V(4,1)-V(1,3)*V(2,4)*V(4,1)- &
 V(1,4)*V(2,1)*V(4,3)+V(1,1)*V(2,4)*V(4,3)+V(1,3)*V(2,1)*V(4,4)- &
 V(1,1)*V(2,3)*V(4,4)

 b(2,4)=V(1,3)*V(2,4)*V(3,1)-V(1,4)*V(2,3)*V(3,1)+ &
 V(1,4)*V(2,1)*V(3,3)-V(1,1)*V(2,4)*V(3,3)-V(1,3)*V(2,1)*V(3,4)+ &
 V(1,1)*V(2,3)*V(3,4)

 b(3,1)=V(2,2)*V(3,4)*V(4,1)-V(2,4)*V(3,2)*V(4,1)+ &
 V(2,4)*V(3,1)*V(4,2)-V(2,1)*V(3,4)*V(4,2)-V(2,2)*V(3,1)*V(4,4)+ &
 V(2,1)*V(3,2)*V(4,4)

 b(3,2)=V(1,4)*V(3,2)*V(4,1)-V(1,2)*V(3,4)*V(4,1)- &
 V(1,4)*V(3,1)*V(4,2)+V(1,1)*V(3,4)*V(4,2)+V(1,2)*V(3,1)*V(4,4)- &
 V(1,1)*V(3,2)*V(4,4)

 b(3,3)=V(1,2)*V(2,4)*V(4,1)-V(1,4)*V(2,2)*V(4,1)+ &
 V(1,4)*V(2,1)*V(4,2)-V(1,1)*V(2,4)*V(4,2)-V(1,2)*V(2,1)*V(4,4)+ &
 V(1,1)*V(2,2)*V(4,4)

 b(3,4)=V(1,4)*V(2,2)*V(3,1)-V(1,2)*V(2,4)*V(3,1)- &
 V(1,4)*V(2,1)*V(3,2)+V(1,1)*V(2,4)*V(3,2)+V(1,2)*V(2,1)*V(3,4)- &
 V(1,1)*V(2,2)*V(3,4)

 b(4,1)=V(2,3)*V(3,2)*V(4,1)-V(2,2)*V(3,3)*V(4,1)- &
 V(2,3)*V(3,1)*V(4,2)+V(2,1)*V(3,3)*V(4,2)+V(2,2)*V(3,1)*V(4,3)- &
 V(2,1)*V(3,2)*V(4,3)

 b(4,2)=V(1,2)*V(3,3)*V(4,1)-V(1,3)*V(3,2)*V(4,1)+ &
 V(1,3)*V(3,1)*V(4,2)-V(1,1)*V(3,3)*V(4,2)-V(1,2)*V(3,1)*V(4,3)+ &
 V(1,1)*V(3,2)*V(4,3)

 b(4,3)=V(1,3)*V(2,2)*V(4,1)-V(1,2)*V(2,3)*V(4,1)- &
 V(1,3)*V(2,1)*V(4,2)+V(1,1)*V(2,3)*V(4,2)+V(1,2)*V(2,1)*V(4,3)- &
 V(1,1)*V(2,2)*V(4,3)

 b(4,4)=V(1,2)*V(2,3)*V(3,1)-V(1,3)*V(2,2)*V(3,1)+ &
 V(1,3)*V(2,1)*V(3,2)-V(1,1)*V(2,3)*V(3,2)-V(1,2)*V(2,1)*V(3,3)+ &
 V(1,1)*V(2,2)*V(3,3)

 if (det /= 0.0_dp) then
   InverV=(1/det)*b
 else
   write(6,*) 'Error in Inverse of Covariance, DET=0., showing V matrix'
   do i=1, 4
     write(6,*) V(i,:)
   end do
   stop
 end if

 end subroutine inverse_matrix_analytic4x4
!-----------------------------------------------------------------


! Dynamic persistence length through tilt and roll:
! 1/Ad' = 1/2(1/At + 1/Ap ) 
!-----------------------------------------------------------------------
! Calculates the absolute value of vector a and stores it in b
 function dynamic_persistence_length2( tilt, roll )
 implicit none
 real(dp) :: dynamic_persistence_length2, tilt, roll

 dynamic_persistence_length2 = 2.0_dp*(tilt*roll/( tilt + roll ) )

 end function dynamic_persistence_length2

!-----------------------------------------------------------------------


! REBUILDING ALGORITHM
!-----------------------------------------------------------------
! Procedure can be seen at:
! Lu, X. J., El Hassan, M. A. & Hunter, C. A. (1997). Structure
! and conformation of helical nucleic acids: rebuilding

! This subroutine calculates the base triads given a set of 
! base step parameters (BSP). These triads are equivalent to the
! triads which would've been calculated if we had
! the atoms coordinates


 subroutine reverse_algorithm_base_triads(Tg_i,rg_i,BSP,nbp)

 integer, intent(in) :: nbp
 real(dp), intent(in) :: BSP(:,:)
 real(dp), intent(out) :: Tg_i(:,:,:),rg_i(:,:)

 integer :: i
 real(dp) :: an12,ax12(3),phi,an3,h(3), & !auxiliar variables
             R1(3,3),R2(3,3),R3(3,3), &   !auxiliar matrices
             Ti_mst(3,3), &               !matrices used in the process
             x(3),y(3),z(3)               !unitary vectors

  x=(/ 1.0_dp, 0.0_dp, 0.0_dp /)
  y=(/ 0.0_dp, 1.0_dp, 0.0_dp /)
  z=(/ 0.0_dp, 0.0_dp, 1.0_dp /)

! The position of first bp is the origin and its triad
! is the identity matrix
  rg_i(:,1)=0.0_dp
  Tg_i(:,1,1)=x
  Tg_i(:,2,1)=y
  Tg_i(:,3,1)=z

  do i=2,nbp
    !lets get the angles
    an12 = BSP(4,i-1)*BSP(4,i-1)+BSP(5,i-1)*BSP(5,i-1)
    an12 = deg_to_rad*sqrt(an12)

    !buckle-propeller axis
    ax12(1) = BSP(4,i-1)*deg_to_rad/an12
    ax12(2) = BSP(5,i-1)*deg_to_rad/an12
    ax12(3) = 0.0_dp

    !angle phi
    phi = dacos(ax12(2))

    !angle 3 in radians (z axis)
    an3 = BSP(6,i-1)*deg_to_rad
    h = cross_product3(ax12,y)
    if ( dot_product(h,z) > 0.0_dp ) then
      phi = abs(phi)
    else
      phi = -abs(phi)
    end if

    !Orientation of bp i from bp 1
    call general_rotation_matrix(R1, z, an3/2.0_dp - phi)
    call general_rotation_matrix(R2, y, an12)
    call general_rotation_matrix(R3, z, an3/2.0_dp + phi)
    Tg_i(:,:,i)=matmul(R1,R2)          !recycling variables
    Tg_i(:,:,i)=matmul(Tg_i(:,:,i),R3)

    !middle frame from reference II, R1 is going to be the same
    call general_rotation_matrix(R2, y, an12/2.0_dp)  !R1 is the same that we already have  
    call general_rotation_matrix(R3, z, phi)
    Ti_mst = matmul(R1,R2)
    Ti_mst = matmul(Ti_mst,R3)

    !position of bp i from bp 1
    h = (/ BSP(1,i-1), BSP(2,i-1), BSP(3,i-1) /)
    rg_i(:,i) = matmul( h,transpose(Ti_mst) )

    !The ones we need
    Tg_i(:,:,i) = matmul( Tg_i(:,:,i-1), Tg_i(:,:,i) )
    rg_i(:,i) = rg_i(:,i-1) + matmul( rg_i(:,i), transpose(Tg_i(:,:,i-1)) )

  end do

 end subroutine reverse_algorithm_base_triads
!-----------------------------------------------------------------------


! GENERAL FUNCTIONS AND SUBROUTINES
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
 !Average and standard deviations are stored in the same arrays in SerraLINE
 function average_std(X,n)
 implicit none
 real(dp) :: average_std(2)
 integer :: n
 real(dp), intent(in) :: X(n)
 real(dp) :: a !auxiliar variable
 integer :: i

 !Get Average
 average_std(1) = sum(X)/real(n,dp)

 !Now Standard deviation
 average_std(2)=0.0_dp

 do i=1,n
   a = X(i) - average_std(1)
   average_std(2) = average_std(2) + a*a
 end do
 average_std(2) = average_std(2)/real(n,dp)
 average_std(2) = sqrt(average_std(2))

 end function 

!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! Calculates the absolute value of vector a and stores it in b
 function absv(a) result (b)
 implicit none
 real(dp), intent(in) :: a(:)
 real(dp) :: b
 integer :: i
  b=0.
  do i=1,size(a)
   b=b+a(i)**2
  end do
  b=sqrt(b)
 end function absv

!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
 !This function normalizes vector v with dimension n
 function normalize_vector(v,n)
 implicit none
 real(dp) :: normalize_vector(n)
 real(dp) :: v(n)
 integer :: n
 integer :: i
 real(dp) :: v_sum

 v_sum=0._dp

 do i=1,n
  v_sum = v_sum + v(i)*v(i)
 end do

 v_sum=sqrt(v_sum)

 do i=1,n
  normalize_vector(i)=v(i)/v_sum
 end do

 end function normalize_vector
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
 !W= U X V
 function cross_product3(U,V)
 real(dp) :: U(3), V(3)
 real(dp) :: cross_product3(3)

 cross_product3(1) = U(2)*V(3)-U(3)*V(2)
 cross_product3(2) = U(3)*V(1)-U(1)*V(3)
 cross_product3(3) = U(1)*V(2)-U(2)*V(1)

 end function cross_product3
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! Constructs the general rotation matrix 'R' ('a')
! describes a rotation of magnitude 'a' about an
! unit vector 'u'
! (In R3)
 subroutine general_rotation_matrix(R,u,a)
 real(dp), intent(in) :: u(3), a
 real(dp), intent(out) :: R(3,3)

 R(1,1)=dcos(a)+(1-dcos(a))*u(1)**2
 R(1,2)=(1-dcos(a))*u(1)*u(2)-u(3)*dsin(a)
 R(1,3)=(1-dcos(a))*u(1)*u(3)+u(2)*dsin(a)
 R(2,1)=(1-dcos(a))*u(1)*u(2)+u(3)*dsin(a)
 R(2,2)=dcos(a)+(1-dcos(a))*u(2)**2
 R(2,3)=(1-dcos(a))*u(2)*u(3)-u(1)*dsin(a)
 R(3,1)=(1-dcos(a))*u(1)*u(3)-u(2)*dsin(a)
 R(3,2)=(1-dcos(a))*u(2)*u(3)+u(1)*dsin(a)
 R(3,3)=dcos(a)+(1-dcos(a))*u(3)**2

 end subroutine general_rotation_matrix
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
! Constructs the a new central fragment vector, which contains
!  all the averages that will be considered
 subroutine central_fragment(overall, dat, n1, n2, nbp, str)
 real(dp), intent(in) :: dat(:)
 real(dp), intent(out) :: overall(:)
 integer, intent(in) :: n1, n2, nbp, str
 integer i, j, p, l, s

 overall=0.0_dp
 l=1
 p=1

 if ( (n2 < n1) .and. str /= 2) then 
   stop "Second range less than first makes no sense if structure is opened  "
 else if (n2 > nbp .or. n1> nbp) then
   stop "Ranges cannot be greater than number of bases"
 end if
 

 !Opened structure
 if (str /= 2) then
   do j=1,nbp-1
     do i=1,nbp-j
       s = i + j
       ! ( n1 <= i,s <= n2)
       if ((i .ge. n1) .and. (i .le. n2) .and. (s .le. n2)) then
         overall(l)=overall(l)+dat(p)
         if (s == n2) l=l+1
       end if
       p=p+1
     end do
   end do

   !Divide by number of data points
   do i=1,n2-n1
     overall(i)=overall(i)/real(n2-n1-i+1,dp) 
   end do

 !Closed structure
 else

   !if first range less than first one
   if (n1 < n2) then
     do j=1,nbp-1
       do i=1,nbp
         s = i + j  !note that it does not matter if s > nbp since it won't
                    !be in the ranges (n1 <= i,s < = n2)
         if ((n1 .le. i) .and. (i .le. n2) .and. (s .le. n2)) then
           overall(l)=overall(l)+dat(p)
           if (s == n2) l=l+1
         end if
         p=p+1
       end do
     end do

     !Divide by number of data points
     do i=1,n2-n1
       overall(i)=overall(i)/real(n2-n1-i+1,dp) 
     end do

   !If first range is greater than second one. This means that
   !we are on the other side of the circle 
   else
    do j=1,nbp-1
       do i=1,nbp
         s = i + j  !note that it does not matter if s > nbp since it won't
                    !be in the ranges (n1 <= i,s < = n2)
                    !The middle condition might not be needed
         if ( s> nbp) s = s -nbp
         if ((n1 .le. i) .and. (i .ge. n2) .and. (s .le. n2)) then
           overall(l)=overall(l)+dat(p)
           if (s == n2) l=l+1
         end if
         p=p+1
       end do
     end do

     !Divide by number of data points
     s = nbp-n1+n2
     do i=1,nbp-n1+n2
       overall(i)=overall(i)/real(nbp-n1+n2-i+1,dp) 
     end do

   end if

 end if

 end subroutine central_fragment
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! Constructs the a new central fragment vector, which contains
!  all the averages and standard deviations that will be considered
 subroutine central_fragment_mean_std(overall, std, dat, n1, n2, nbp, str)
 real(dp), intent(in) :: dat(:)
 real(dp), intent(out) :: overall(:), std(:)
 integer, intent(in) :: n1, n2, nbp, str
 integer i, j, p, l, s

 overall=0.0_dp
 std = 0.0_dp

 if ( (n2 < n1) .and. str /= 2) then 
   stop "Second range less than first makes no sense if structure is opened  "
 else if (n2 > nbp .or. n1> nbp) then
   stop "Ranges cannot be greater than number of bases"
 end if
 

 !Opened structure
 if (str /= 2) then
   !mean
   l=1
   p=1
   do j=1,nbp-1
     do i=1,nbp-j
       s = i + j
       ! ( n1 <= i,s <= n2)
       if ((i .ge. n1) .and. (i .le. n2) .and. (s .le. n2)) then
         overall(l)=overall(l)+dat(p)
         if (s == n2) l=l+1
       end if
       p=p+1
     end do
   end do

   !Divide by number of data points
   do i=1,n2-n1
     overall(i)=overall(i)/real(n2-n1-i+1,dp) 
   end do

   !std
   l=1
   p=1
   do j=1,nbp-1
     do i=1,nbp-j
       s = i + j
       if ((i .ge. n1) .and. (i .le. n2) .and. (s .le. n2)) then
         std(l)=std(l)+(dat(p)-overall(l))**2
         if (s == n2) l=l+1
       end if
       p=p+1
     end do
   end do

   !Divide by number of data points
   do i=1,n2-n1
     std(i)=sqrt(std(i)/real(n2-n1-i+1,dp))
   end do

 !Closed structure
 else

   !if first range less than first one
   if (n1 < n2) then
 
    !mean
    l=1
    p=1  
    do j=1,nbp-1
       do i=1,nbp
         s = i + j  !note that it does not matter if s > nbp since it won't
                    !be in the ranges (n1 <= i,s < = n2)
         if ((n1 .le. i) .and. (i .le. n2) .and. (s .le. n2)) then
           overall(l)=overall(l)+dat(p)
           if (s == n2) l=l+1
         end if
         p=p+1
       end do
     end do

     !Divide by number of data points
     do i=1,n2-n1
       overall(i)=overall(i)/real(n2-n1-i+1,dp) 
     end do

     !std
     l=1
     p=1
     do j=1,nbp-1
       do i=1,nbp
         s = i + j
         if ((n1 .le. i) .and. (i .le. n2) .and. (s .le. n2)) then
           std(l)=std(l)+(dat(p)-overall(l))**2
           if (s == n2) l=l+1
         end if
         p=p+1
       end do
     end do

     !Divide by number of data points
     do i=1,n2-n1
       std(i)=sqrt(std(i)/real(n2-n1-i+1,dp)) 
     end do

    !If first range is greater than second one. This means that
    !we are on the other side of the circle 
    else

    !mean
    l=1
    p=1   
    do j=1,nbp-1
       do i=1,nbp
         s = i + j  !note that it does not matter if s > nbp since it won't
                    !be in the ranges (n1 <= i,s < = n2)
                    !The middle condition might not be needed
         if ( s> nbp) s = s -nbp
         if ((n1 .le. i) .and. (i .ge. n2) .and. (s .le. n2)) then
           overall(l)=overall(l)+dat(p)
           if (s == n2) l=l+1
         end if
         p=p+1
       end do
     end do

     !Divide by number of data points
     s = nbp-n1+n2
     do i=1,nbp-n1+n2
       overall(i)=overall(i)/real(nbp-n1+n2-i+1,dp) 
     end do

     !std
     l=1
     p=1
     do j=1,nbp-1
       do i=1,nbp
         s = i + j
         if ( s > nbp) s = s -nbp
         if ((n1 .le. i) .and. (i .ge. n2) .and. (s .le. n2)) then
           std(l)=std(l)+(dat(p)-overall(l))**2
           if (s == n2) l=l+1
         end if
         p=p+1
       end do
     end do

     !Divide by number of data points
     s = nbp-n1+n2
     do i=1,nbp-n1+n2
       std(i)=sqrt(std(i)/real(nbp-n1+n2-i+1,dp))
     end do


   end if !n1< n2

 end if   !open or closed

 end subroutine central_fragment_mean_std
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!Extracts sublength from datain and stores it in subdata.
!subdata has to be already allocated.
!The number of lengths depend on the type of structure (str)
!This subroutine handles datas that contain means and standard
!devations, like BSP, and strucp. For other structures with only mean
!values use extract_sublength_2d

 subroutine extract_sublength_3d( datain, subdata, nbp, sublength, str )
 implicit none
 real(dp), intent(in) :: datain(:,:,:)
 real(dp), intent(out) :: subdata(:,:,:)
 integer, intent(in) :: nbp, sublength, str
 integer :: i, j, k, l

 !If linear structure
 if (str == 1) then 

   l = 0
   !To directly select the sublength, will require calculate a factorial,
   !since we don't expect huge amounts of data, it might be easier to
   !use loops
   do j = 1, nbp-1
     
     if ( j == sublength ) then

       do k = 1, nbp-sublength !or nbp-j
         l = l + 1
         subdata(:,:,k) = datain(:,:,l)
       end do !k 

       exit !escape! we have what we wanted

     end if

     do i = 1, nbp-j

       l = l +1
       
     end do !i 

   end do !j

 !If circular
 else if (str == 2) then

   !Here, we can directly know which index we want
   l = nbp*(sublength-1)

   do k = 1, nbp
     l = l + 1
     subdata(:,:,k) = datain(:,:,l)
   end do

 !Something else? That's wrong
 else

   stop "Error in extract_sublength. Unidentified type of structure"

 end if


 end subroutine extract_sublength_3d
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!Extracts sublength from datain and stores it in subdata.
!subdata has to be already allocated.
!The number of lengths depend on the type of structure (str)
!This subroutine handles data that contain means
!like elasp and avstrp. For other structures with mean and stds
!values use extract_sublength_3d

 subroutine extract_sublength_2d( datain, subdata, nbp, sublength, str )
 implicit none
 real(dp), intent(in) :: datain(:,:)
 real(dp), intent(out) :: subdata(:,:)
 integer, intent(in) :: nbp, sublength, str
 integer :: i, j, k, l

 !If linear structure
 if (str == 1) then 

   l = 0
   !To directly select the sublength, will require calculate a factorial,
   !since we don't expect huge amounts of data, it might be easier to
   !use loops
   do j = 1, nbp-1

     if ( j == sublength ) then

       do k = 1, nbp-sublength !or nbp-j
         l = l + 1
         subdata(:,k) = datain(:,l)
       end do !k 

       exit !escape! we have what we wanted

     end if

     do i = 1, nbp-j

       l = l +1
       
     end do !i 

   end do !j

 !If circular
 else if (str == 2) then 

   !Here, we can directly know which index we want
   l = nbp*(sublength-1)

   do k = 1, nbp
     l = l + 1
     subdata(:,k) = datain(:,l)
   end do

 !Something else? That's wrong
 else

   stop "Error in extract_sublength. Unidentified type of structure"

 end if


 end subroutine extract_sublength_2d
!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
! Calculates simple linear regression
! X, Y => axis, a => intercept, b => slope, c => confidence interval of b
 subroutine simple_linear_regression(Y,X,a,b,c)

 real(dp), intent(in) :: X(:),Y(:)
 real(dp), intent(out) :: a,b,c
 real(dp) :: avx,avy,cxy,vx,e2,sb,tn_2
 integer :: i, n

 n = size(X)
 if (size(Y) /= n) stop "Different sizes of X & Y in linear regression (a=0)"

 avx=sum(X(:))/real(n,dp)
 avy=sum(Y(:))/real(n,dp)

!vx & cxy
 vx=0.0_dp
 cxy=0.0_dp
 do i=1,n
  vx=vx+(X(i)-avx)**2
  cxy=cxy+(X(i)-avx)*(Y(i)-avy)
 end do
!Note that vx and cxy weren't divided by n cause they would be cancelled in b 

!b the slope
 b=cxy/vx

!a the intercept
 a=avy-b*avx

!Calculate confidence interval

!e2 sum of squared residual
 e2=0.0_dp
 do i=1,n
  e2 = e2 +  ( Y(i) -a -b*X(i) )**2
 end do

 if (n-2 <= 0) stop "n not large enough in linear regression"
 sb = sqrt( ( e2/(n-2) )/vx )  !standard error of b

 tn_2 = tstudent_th_quantil(n-2) 
 c = tn_2*sb                   !Got the confidence interval

 end subroutine simple_linear_regression
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
! Calculates simple linear regression with alpha=0
! X, Y => axis, a => intercept, b => slope, c => confidence interval
 subroutine simple_linear_regression_a_0(Y,X,a,b,c)

 real(dp), intent(in) :: X(:),Y(:)
 real(dp), intent(out) :: a,b,c
 real(dp) :: avxy,avx2,avx,vx,e2,sb,tn_2
 integer :: i, n

 n = size(X)
 if (size(Y) /= n) stop "Different sizes of X & Y in linear regression (a=0)"
 a=0.0_dp

!avx2 & avxy
 avx2=0.0_dp
 avxy=0.0_dp

 do i=1,n
   avx2=avx2+X(i)**2
   avxy=avxy+X(i)*Y(i)
 end do

 avx2=avx2/real(n,dp)
 avxy=avxy/real(n,dp)

!b the slope
 b=avxy/avx2

!Calculate confidence interval
 avx=sum(X(:))/real(n,dp)

!vx & cxy
 vx=0.0_dp
 do i=1,n
  vx=vx+(X(i)-avx)**2
 end do

!e2 sum of squared residual
 e2=0.0_dp
 do i=1,n
  e2 = e2 + ( Y(i) -a -b*X(i) )**2
 end do

 if (n-2 <= 0) stop "n not large enough in linear regression"
 sb = sqrt( ( e2/(n-2) )/vx )  !standard error of b
 tn_2 = tstudent_th_quantil(n-2)
 c = tn_2*sb

 end subroutine simple_linear_regression_a_0
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! Calculates simple linear regression with alpha=1
! X, Y => axis, a => intercept, b => slope, c => confidence interval
 subroutine simple_linear_regression_a_1(Y,X,a,b,c)

 real(dp), intent(in) :: X(:),Y(:)
 real(dp), intent(out) :: a,b,c
 real(dp) :: avx, avxy,avx2,e2,vx,sb,tn_2
 integer :: i, n

 n = size(X)
 if (size(Y) /= n) stop "Different sizes of X & Y in linear regression (a=0)"
 a=1.0_dp

 !average x
 avx = sum(X(:))/real(n,dp)

 !avx2 & avxy
 avx2=0.0_dp
 avxy=0.0_dp

 do i=1,n
   avx2=avx2+X(i)**2
   avxy=avxy+X(i)*Y(i)
 end do

 avx2=avx2/real(n,dp)
 avxy=avxy/real(n,dp)

!b the slope
 b=(avxy-avx)/avx2

!Calculate confidence interval
!vx 
 vx=0.0_dp
 do i=1,n
  vx=vx+(X(i)-avx)**2
 end do

!e2 sum of squared residual
 e2=0.0_dp
 do i=1,n
  e2 = e2 + ( Y(i) -a -b*X(i) )**2
 end do

 if (n-2 <= 0) stop "n not large enough in linear regression"
 sb = sqrt( ( e2/(n-2) )/vx )  !standard error of b

 tn_2 = tstudent_th_quantil(n-2)
 c =  tn_2*sb

 end subroutine simple_linear_regression_a_1
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! Assigns the corresponding th quantil of the t-student distribution.
 function tstudent_th_quantil(df)
 real(dp) :: tstudent_th_quantil
 integer :: df, n, area(1)

 if (df < 1) then
   write(6,*) "Error in tstudent_th_quantil function."
   write(6,*) "-> Degrees of freedom less than 1?"
 end if

 !Assign equivalent degree of freedom to corresponding index in array
 ! t_quantil.
 if (df > size(t_quantil,2) ) then
   n = size(t_quantil,2)
 else
   n = df
 end if

 area = minloc( abs(t_percent - confidence_l) ) !Find closest area 

 tstudent_th_quantil = t_quantil(area(1),n) !assign quantil

 end function tstudent_th_quantil
!-----------------------------------------------------------------------

 end module functions_mod

