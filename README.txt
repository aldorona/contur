# contur
CONTUR code for contoured nozzle design by J.C. Sivells ADEC TR 78 63 December 1978

This respository contains the CONTUR code in 
Appendix D of J.C. Sivells, A computer program for the aerodynamic
design of axisymmetric and planar nozzles for supersonic and hypersonic
wind tunnels, ARO Inc., a Sverdrup Corporation Company, ADEC TR 78 63,
December 1978.

The code uses a combination of analytical solutions, the method of 
characteristics, and centerline distributions in order to calculate the
divergent section of a convergent-divergent de Laval nozzle.
 
Folder sivells contains a FORTRAN77 source code which runs through 
a series of 16 subroutines and uses 7 user-defined input cards
describing the flow conditions of the desired nozzle profile.

Folder src contains the FORTRAN90 source code version of folder sivells

Sivells reports the input cards and output file of a MACH 4 
axisymmetric nozzle. The version of the code in this repository was
run in May 2019 on this test case and the same output was obtained.

Florentina-Luiza Zavalan, under the guide of project
supervisor Aldo Rona, typeset Appendix D from the freely available
source code listing at 
https://apps.dtic.mil/dtic/tr/fulltext/u2/a062944.pdf
Significant effort was put in to interpret the low-resolution scan and
to disambigue similar typographical symbols (e.g. * from +, 5 from S).

Aldo Rona translated the code in FORTRAN90. If you find any bug
while testing, please report it to the authors: 
aldo.rona@le.ac.uk, flz1@leicester.ac.uk

License: See LICENSE.txt. Code users must acknowledge the
provenance of the listing by using the following acknowledgement in
their published work: "This work
used the CONTUR source code by L.F. Zavalan and A. Rona, based on the
computer program by J.C. Sivells".

Compilation instructions

Linux/Unix users:
1. Create a new directory, suggested name: sivells
2. Download in the SAME directory the bundle of 20 *.f source files,
   one input.txt file, and makefile.
3. You are advised to also download Sivells.pdf, which is article by
   Sivells, and output.txt, which is the sample output from the MACH 4 
   test case. These files are in folders sivells/ and docs/.
4. Compile the code by just typing: make. This will create the object
   files *.o and the executable MAIN.exe
5. Run the executable MAIN.exe

Windows users:
1. Download and install Microsoft Visual Studio
2. Open Microsoft Visual Studio and create a new Project, suggested name: contur
3. Move input.txt and input2d.txt from the /src directory into the
   'resources' directory of the Project.
4. Move all .f files from the /src directory to the 'source' directory of the Project.
   Make sure you only move the .f files from /src, which have all small case names, do
   not add any of the .f files from the /sivells directory, which have CAPITAL CASE names,
   as this will over-define symbols in the compilation.
5. Right click on the makefile you have downloaded and open it with a text editor.
   Do not double click on it. Check the list of files you now have in 'source' 
   of your Project. You should only have files with root names <root>.f matching 
   the object list root names <root>.o in makefile, e.g. axial.f -> axial.o. You should
   also have main.f in your 'source'. You should not 
   have any .txt, makefile, or .pdf files in your 'source'.
6. Remove makefile. You do not need this file for compiling under Windows.
7. Select Compile and Run from the taskbar of Microsoft Visual Studio

Frequently Asked Questions
1. You may need to let makefile know what compiler you are using. To
do so, use a text editor to edit makefile and change FC=ifort to your
own compiler, e.g. gfort. Run the command "man -k fortran" to find out
what compiler is installed on your system.

2. You may need to change the compiler flags according to what is
available from your compiler. In unix/linux, you can type
"man mycompiler" to find out what flags are available.

08 May 2019 - A. Rona, F.L. Zavalan

Card notes:

Card 5: Use an integer value for XJ as JX=XJ and JX is integer

Code updates:
14May2019: AR moved label 5 from line BOU 72 to BOU 69 to avoid jumping
into the DO IV=1,IW loop.
16May2019: AR variables ns and nc removed from AXIAL.f as not used
16May2019: AR removed M from common block of NEO.f
11Jun2019: AR code changed to f90. Added input2d.txt as sample input
file to generate a two-dimensional contoured nozzle. To use rename
input2d.txt as input.txt and rund the code

-----
end of notes
