___________________________________________________________________
                        \       _
			'       '
                \       ---------        /
		|--------SERRANA--------|
                        ---------
                        |       |
                        -       -
     by Victor Velasco with the guidance of Agnes Noy
___________________________________________________________________
___________________________________________________________________
                       OVERVIEW
                       --------

SerraNA is composed by three executables: SerraNA, Analysis and Extract

The whole program calculates flexibility parameters from simulations of 
nucleic acid structures (circular or linear) at different length scales, 
by following the  Length-Dependent Elastic Model (LDEM).

The methodology implemented, allows to study how global elastic
properties arise from the base-pair (bp) level, and to infer elastic
constants that describe the fragment's overall flexibility.
These elastic constants correspond to torsional, stretch modulus,
and the persistence length.

For more information look the manual and the paper:

Velasco, Burman, Shepherd, Leake, Golestanian and Noy (2020), biorxiv:
"SerraNA: a program to infer elastic constants from local to global nucleic acids simulation data"
This is the paper to cite! :). For more information, please visit:
https://agnesnoylab.wordpress.com/

___________________________________________________________________
                       REQUIREMENTS
                       ------------

You only need a fortran compiler (gfortran), that's it.

___________________________________________________________________

                       COMPILATION
                       -----------

Type make SerraNA for compiling SerraNA.
     make Analysis to compile Analysis section.
     make Extract for compiling the extraction program.
     or 
     make all for all of them.

___________________________________________________________________

                       SerraNA
                       -------

First main process that carries the calculation of structural and 
flexible parameters flexibility for every possible pair of bp.
For running SerraNA type:

./SerraNA < s_NA.in

A trajectory file and a topology file in AMBER style format is needed 
(10F8.3 for the trajectory). The files can contain ions or other residues 
and SerraNA will ignore them.

s_NA.in is the input file that indicates: 
  1.- The path for topology and trajectory
  2.- If the structure is double-stranded ("2") or if it is 
      single-stranded ("1").
  3.- If the structure is linear ("1") or closed ("2"). 
      For linear NA, SerraNA ignores the two base-pairs at each end for
      avoiding end effects.

SerraNA creates 4 ouputs:

- BPP.out which have the base-pair parameters as they are calculated in 3DNA
  
- BSP.out which have the base-step parameters as they are calculated in 3DNA
  plus the total bending

- structural_parameters.out which have variables describing the geometry of the
  DNA molecule for all possible sub-fragments:
  1.- Added-rise, added-slide, added-rise, roll, tilt and twist as an extension
      of CEHS algorithm to sequences > 2 bp.
  2.- End-to-end distance and contour length
  3.- Bending, bending**2, and directional correlation (D correlation) 
  4.- From averaged structure: Bending (AVSTR B), bending**2 (AVSTR B**2) 
      and directional correlation (AVSTR D C) 

- elastic_parameters.out have elastic constants for Stretch (pN), Twist (nm), 
  Roll (nm), Tilt (nm), as well as their couplings (nm), together with the
  variance and partial variance  of the end-to-end distance (angstroms).

___________________________________________________________________

                       Analysis
                       --------

Second main process which analyze the parameters calculated by SerraNA 
and estimates global elastic constants. For running Analysis type:

./Analysis < ov_NA.in

ov_NA.in is the input file that indicates:
   1.- The path to elastic_parameters.out and structural_parameters.out 
   2.- The part of the molecule that will be used to calculate global elastic
       constants. Two ranges should be provided:
       (i) The first one defines the region of the molecule used (from bp "a" 
           to bp "b"). 
       (ii) The second indicates the sub-lengths considered, being from number 
            of bp-steps "c" to number of bp-steps "d". Note that c > 0 and d <= b-a.
       DEFAULT OPTION, which is specified by a=b=c=d=0, uses the whole fragment
       and applies the recommended methodology described on the paper and the manual.

Analysis prints information on screen regarding the global elastic constants,
stretch modulus, torsion an persistence lenght.

___________________________________________________________________

                       Extract
                       -------

The Extract program process SerraNA ouputs, elastic_parameters.out and 
structural_parameters.out, creating simple files ready to plot. 
You can filter a particular sublength or you can extract averages 
and standard deviations as a function of length from a particular region 
 to obtain plots similar to Figure 3 of Velaco et al. Extract also can 
process BPP.out and BSP.out for easier plotting. For running Analysis type:

./Extract < ex_NA.in

ex_NA.in is the input file that indicates:
   1.- Path to either BPP, BSP, structural or elastic parameters file.
       If you selected to extract BPP or BSP, then all other inputs will be ignored.
   2.- Type "0" for extracting a sub-length or "1" for getting avg+-sd as 
       a function of length.
   3.- (i) If you type 0, then indicate the length (l) you want to process,
           which should be 0 < l < N, where N is the number of bp-steps.
       (ii) If you typed 1, then indicate the region (a,b) from which you want
            to extract avg+-sd as a function of length.
            If it is linear DNA, then 0 < a < b < N
            If it is circular DNA, then both a < b or b < a, are valid
            DEFAULT OPTION, a=b=0, consider all possible lengths in the whole 
            fragment.

The program can create different types of outputs:

- BPP_plot.out, if BPP.out is processed
- BSP_plot.out, if BSP.out is processed
- structural_lmer.out, if structural_parameters.out is processed to extract a 
  particular sub-length l
- structural_[a:b].out, if structural_parameters.out is processed to extract a 
  length-dependence for a particular sub-fragment
- structural_plot.out, if a=b=0
- elastic_lmer.out, if elastic_parameters.out  is processed to extract a 
  particular sub-length l
- elastic_[a:b].out, if elastic_parameters.out is processed to extract a
  length-dependence for a particular sub-fragment
- elastic_plot.out, if a=b=0


All formats can be modified in parms.f90, in section INPUT/OUTPUT FORMATS,
where for BPP, BSP, structural and elastic, _2 corresponds to parameters
by sublength and _3 correspond to overalls.

Finally, we provide an example of the python script used to make the plots at
Velasco et al.
___________________________________________________________________

                         TEST
                       -------

We provide an example/test with topology and trajectory files corresponding
to the 32-mer from Velasco et al.

The first step is to compile the programs by typing "make" and then
run the main code (SerraNA.f90) by typing:

./SerraNA < s_NA.in

This will produce the following files: BPP.out BSP.out 
structural_parameters.out and elastic_parameters.out

The input files *.in have the options set for processing the 32-mer trajectory.

You can edit ov_NA.in and run the Analysis for experimenting different results.

Once the structural and elastic parameters are calculated, one can use the
Extract process to extract the data that one is interested in or to filter
a particular length (l) or a region [a,b].

The bash script run-example.sh runs and edits ex_NA.in multiple times to extract
structural_parameters.out at the lengths of 2, 12 and 22 bp, and the avg-std
of elastic_parameters.out at regions [1,6], [7,12], [13,18] and [19,24], and
finally, the overalls for the elastic parameters. Then it executes a python script
plot-example.py to plot 1) bending angles for 2mer, 12mer and 22mer.
                        2) Twist elastic constant for the 4 regions as a 
                           function of length.
                        3) Stretch modulus as a function of length for the 
                           whole fragment.

This produces a pdf image (my_first_result.pdf) with the 3 subplots.
___________________________________________________________________

                        LICENSE
                       ---------

    SerraNA is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published 
    by the Free Software Foundation, either version 3 of the License, or
    any later version.

    SerraNA is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.

___________________________________________________________________
