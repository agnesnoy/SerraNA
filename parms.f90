!------------------------------------------------------------- 25/08/2021
! This module contains parameters needed for SerraNA
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

  module parms
  implicit none

  !GENERAL PARAMETERS
!-----------------------------------------------------------------------

  integer, parameter :: sp= kind(0.0), dp = kind(0.d0), &
                        n_R=9, n_Y=6    !purine and pyrimidine ring atoms

  real(dp), parameter :: pi = datan(1.0_dp)*2.0_dp*2.0_dp, &
                       & rad_to_deg = 180.0_dp/pi, deg_to_rad = pi/180.0_dp, &
                       & one = 1.0_dp, eps = 1.0E-12_dp, &  !^-12 is enough
                       & w_0 = 2.0_dp*pi/34.0_dp, &         !factor w_0 in angstroms
                       & bnm = 0.34_dp, &                   !bp rise in nanometers
                       & bKTpN=4.14E-21_dp, &               !Boltzmann constant?
                       & confidence_l = 0.70_dp             !confidence level for linear regression (do not change it)

  !BASE PARAMETERs
!-----------------------------------------------------------------------
  !PURINES
  real (dp), dimension(3,n_R), parameter :: &

  ! Guanine
  & G_b = reshape( (/ -1.289_dp,  4.551_dp, -0.000_dp , &    !N9
  &  0.023_dp,  4.962_dp,  0.000_dp, &                       !C8
  &  0.870_dp,  3.969_dp, -0.000_dp, &                       !N7
  &  0.071_dp,  2.833_dp,  0.000_dp, &                       !C5
  &  0.424_dp,  1.460_dp,  0.000_dp, &                       !C6
  & -0.700_dp,  0.641_dp, -0.000_dp, &                       !N1
  & -1.999_dp,  1.087_dp, -0.000_dp, &                       !C2
  & -2.342_dp,  2.364_dp,  0.001_dp, &                       !N3
  & -1.265_dp,  3.177_dp, -0.000_dp /), &                    !C4
  & (/ 3,n_R /) ), &

  ! Adenine
  & A_b = reshape( (/ -1.291_dp,  4.498_dp, -0.000_dp , &    !N9
  &  0.024_dp,  4.897_dp,  0.000_dp, &                       !C8
  &  0.877_dp,  3.902_dp,  0.000_dp, &                       !N7
  &  0.071_dp,  2.771_dp,  0.000_dp, &                       !C5
  &  0.369_dp,  1.398_dp,  0.000_dp, &                       !C6
  & -0.668_dp,  0.532_dp, -0.000_dp, &                       !N1
  & -1.912_dp,  1.023_dp,  0.000_dp, &                       !C2
  & -2.320_dp,  2.290_dp,  0.000_dp, &                       !N3
  & -1.267_dp,  3.124_dp,  0.000_dp /), &                    !C4
  & (/ 3,n_R /) )

  !PYRIMIDINES
  real (dp), dimension(3,n_Y), parameter :: &

  ! Uracil
  & U_b = reshape( (/ -1.284_dp,  4.500_dp, 0.000_dp , &      !N1
  & -1.462_dp,   3.131_dp, -0.000_dp, &                       !C2
  & -0.302_dp,   2.397_dp,  0.000_dp, &                       !N3
  &  0.989_dp,   2.884_dp, -0.000_dp, &                       !C4
  &  1.089_dp,   4.311_dp,  0.000_dp, &                       !C5
  &  -0.024_dp,  5.053_dp, -0.000_dp /), &                    !C6
  & (/ 3,n_Y /) ), &

  ! Cytosine
  & C_b = reshape( (/ -1.285_dp,  4.542_dp, -0.000_dp , &    !N1
  & -1.472_dp,  3.158_dp,  0.000_dp, &                       !C2
  & -0.391_dp,  2.344_dp, -0.000_dp, &                       !N3
  &  0.837_dp,  2.868_dp, -0.000_dp, &                       !C4
  &  1.056_dp,  4.275_dp,  0.000_dp, &                       !C5
  & -0.023_dp,  5.068_dp, -0.000_dp /), &                    !C6
  & (/ 3,n_Y /) ), &

  ! Thymine
  & T_b = reshape( (/ -1.284_dp,  4.500_dp,  0.000_dp , &    !N1
  & -1.462_dp,  3.135_dp,  0.000_dp, &                       !C2
  & -0.298_dp,  2.407_dp, -0.000_dp, &                       !N3
  &  0.994_dp,  2.897_dp, -0.000_dp, &                       !C4
  &  1.106_dp,  4.338_dp,  0.000_dp, &                       !C5
  & -0.024_dp,  5.057_dp, -0.000_dp /), &                    !C6
  & (/ 3,n_Y /) )

  !BASE LABELS
!-----------------------------------------------------------------------

  character(len=4), dimension(6), parameter :: &

  !ADENINE
  & A_l = ["A   ", "A5  ", "A3  ", "DA  ", "DA5 ", "DA3 "], &

  !GUANINE
  & G_l = ["G   ", "G5  ", "G3  ", "DG  ", "DG5 ", "DG3 "], &

  !CYTOSINE
  & C_l = ["C   ", "C5  ", "C3  ", "DC  ", "DC5 ", "DC3 "], &

  !THYMINE
  & T_l = ["T   ", "T5  ", "T3  ", "DT  ", "DT5 ", "DT3 "], &

  !URACIL
  & U_l = ["U   ", "U5  ", "U3  ", "DU  ", "DU5 ", "DU3 "]
!-----------------------------------------------------------------------


  !RING ATOMS LABELS
!-----------------------------------------------------------------------
  character(len=4), parameter :: &

  N1 = "N1  ", C2 = "C2  ", N3 = "N3  ", C4 = "C4  ", C5 = "C5  ", &
  C6 = "C6  ", N7 = "N7  ", C8 = "C8  ", N9 = "N9  "

  !INPUT/OUTPUT FORMATS
!-----------------------------------------------------------------------
  character(50), parameter :: &

  !SERRANA OUTPUTS
  !---------------
  !BPP
  & F_BPP_1 = "(1A10,6A20)", F_BPP_2 = "(I6,4A1,12F10.3)", &
  & F_BPP_3 = "(1A10,12F10.3)", &

  !BSP 
  & F_BSP_1 = "(1A15,7A20)", F_BSP_2 = "(1I4,1A1,1I4,6A1,14F10.3)", &
  & F_BSP_3 = "(1A15,14F10.3)", &

  !STRUCTURAL PARAMETERS 
  & F_STRP_1 = "(1A15,11A20,3A20)", &
  & F_STRP_2 = "(1I4,1A1,1I4,6A1,20F10.3,2F10.6,2F20.3,1F20.6)", &
  & F_STRP_3 = "(1A15,20F10.3,2F10.6,4F10.3,2F10.6)", &

!  & F_STRP_1 = "(1A15,11A40,3A20)", &
!  & F_STRP_2 = "(1I4,1A1,1I4,6A1,20F20.3,2F20.6,2F20.3,1F20.6)", &
!  & F_STRP_3 = "(1A15,20F20.3,2F20.6,4F10.3,2F10.6)", &

  !ELASTIC PARAMETERS
  & F_ELAP_1 = "(1A15,10A20,1A20, 2A20)", &
  & F_ELAP_2 = "(1I4,1A1,1I4,6A1,10F20.3,1F20.3,1F20.5,1F20.5)", &
  & F_ELAP_3 = "(1A15,20F10.3,2F10.3,2F10.5,2F10.5)", &

  !OVERALL ELASTIC CONSTANTS
  & F_OV_1 =  "(6A15)", &
  & F_OV_2 =  "(A15,1F15.3,3A15,1F15.3)", &
  & F_OV_3 =  "(A15,1F15.3,2F15.5,1F15.3,A15)", &
  & F_OV_4 =  "(A15,1F15.3,4A15)", &
  & F_OV_5 =  "(A,1I7,A,1I7,A,1I7,A,1I7)", &
  & F_OV_6 =  "(A,1I7,A,1I7)"

!----------------------------------------------------------------------- 
!----------------------------------------------------------------------- 
!!!!!!!!!!!!!!!!!!DO NOT EDIT THE FOLLOWING PARAMETERS!!!!!!!!!!!!!!!!!!

  !T STUDENT QUANTILS: This is used for calculating confidence 
  !intervals when using linear regression
!-----------------------------------------------------------------------

  !PERCENTAGES. This corresponds to the y-axis in the t_quantil 
  !             parameter. It will help to get the approximate value.
  real(dp), dimension(1), parameter :: &

  t_percent = 0.7_dp

  !QUANTILS. And the quantils. This will be huge...
  !          The x-axis corresponds to the area (determined by
  !          t_percent parameter array) and y-axis to the degrees
  !          of freedom
  !          In this case, we will only use confidence interval of 70%
  real (dp), dimension(1,100), parameter :: &

  t_quantil = reshape( (/  1.96261050553_dp, &
 & 1.38620656017_dp, &
 & 1.24977810368_dp, &
 & 1.18956685235_dp, &
 & 1.15576734288_dp, &
 & 1.13415693067_dp, &
 & 1.11915912836_dp, &
 & 1.10814544456_dp, &
 & 1.09971619639_dp, &
 & 1.09305807359_dp, &
 & 1.08766638035_dp, &
 & 1.08321142046_dp, &
 & 1.07946873704_dp, &
 & 1.07628024458_dp, &
 & 1.07353139558_dp, &
 & 1.07113716328_dp, &
 & 1.06903311062_dp, &
 & 1.06716951554_dp, &
 & 1.06550739859_dp, &
 & 1.06401577116_dp, &
 & 1.06266968707_dp, &
 & 1.06144884144_dp, &
 & 1.06033654227_dp, &
 & 1.05931892292_dp, &
 & 1.05838439553_dp, &
 & 1.05752318312_dp, &
 & 1.05672697701_dp, &
 & 1.05598869993_dp, &
 & 1.05530224628_dp, &
 & 1.05466234517_dp, &
 & 1.054064417_dp, &
 & 1.05350446369_dp, &
 & 1.05297897928_dp, &
 & 1.05248487644_dp, &
 & 1.05201942582_dp, &
 & 1.05158020564_dp, &
 & 1.05116505954_dp, &
 & 1.05077206119_dp, &
 & 1.05039948443_dp, &
 & 1.05004577793_dp, &
 & 1.04970954359_dp, &
 & 1.04938951808_dp, &
 & 1.049084557_dp, &
 & 1.0487936212_dp, &
 & 1.04851576493_dp, &
 & 1.04825012562_dp, &
 & 1.04799591491_dp, &
 & 1.04775241082_dp, &
 & 1.04751895094_dp, &
 & 1.0472949264_dp, &
 & 1.04707977655_dp, &
 & 1.04687298431_dp, &
 & 1.04667407199_dp, &
 & 1.04648259763_dp, &
 & 1.04629815171_dp, &
 & 1.04612035424_dp, &
 & 1.04594885211_dp, &
 & 1.04578331682_dp, &
 & 1.04562344228_dp, &
 & 1.04546894304_dp, &
 & 1.04531955249_dp, &
 & 1.04517502137_dp, &
 & 1.04503511638_dp, &
 & 1.04489961891_dp, &
 & 1.04476832393_dp, &
 & 1.04464103891_dp, &
 & 1.04451758291_dp, &
 & 1.04439778572_dp, &
 & 1.04428148707_dp, &
 & 1.04416853591_dp, &
 & 1.04405878977_dp, &
 & 1.04395211414_dp, &
 & 1.04384838194_dp, &
 & 1.04374747302_dp, &
 & 1.04364927364_dp, &
 & 1.04355367614_dp, &
 & 1.04346057846_dp, &
 & 1.04336988381_dp, &
 & 1.04328150034_dp, &
 & 1.04319534082_dp, &
 & 1.04311132236_dp, &
 & 1.04302936611_dp, &
 & 1.04294939709_dp, &
 & 1.04287134389_dp, &
 & 1.04279513846_dp, &
 & 1.04272071598_dp, &
 & 1.04264801458_dp, &
 & 1.04257697525_dp, &
 & 1.04250754164_dp, &
 & 1.04243965992_dp, &
 & 1.04237327864_dp, &
 & 1.04230834858_dp, &
 & 1.04224482268_dp, &
 & 1.04218265587_dp, &
 & 1.04212180501_dp, &
 & 1.04206222874_dp, &
 & 1.04200388744_dp, &
 & 1.04194674309_dp, &
 & 1.04189075925_dp, &
 & 1.0418359009_dp /), &
 & (/ 1,100 /) )

!----------------------------------------------------------------------- 
!----------------------------------------------------------------------- 

  contains

!-----------------------------------------------------------------------

  end module parms
