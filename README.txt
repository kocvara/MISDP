1. CODES RELATED TO THE ARTICLE

"Truss topology design with integer variables made easy"

by Michal Kocvara, University of Birmingham, kocvara@maths.bham.ac.uk
see www.optimization-online.org

2. THIS DIRECTORY CONTAINS

- this README.txt file

- subdirectory GEO with many examples of truss ground structures
written in ascii format

- kobum.m ...code for reading the truss data and creating the "par" 
             structure

- minvol.m ...code for minimum volume problem with continuous variables
- minvolbb.m  ...minimum volume, single load, binary variables
- minvolbb_m.m ... minimum volume, multiple load, binary variables
- minvolbb_v.m ...minimum volume, single load, vibration constraints, 
                  binary variables
- minvolbb_v_m ...minimum volume, multiple-mass problem, binary variables

- pic_ini.m, pic.m, topo.m ...auxiliary files

3. HOW TO SOLVE A PROBLEM

in Matlab:

>> par = kobum('GEO/t3x3.geo');
>> minvol  (or minvolbb, minvolbb_m, ...)

4. WHAT SHOULD BE CHANGED MANUALLY

Some parameters, such as \gamma or \lambdamin are different for each 
example (see the article). These must be changed manually in the respective
minvol*.m file

The loads for the multiple load case are not included in the *.geo file and 
must be set manually in minvol_m.m

NOTICE that with incorrectly chosen parameters, the problem may be 
infeasible or badly scaled

5. WHAT IS NEEDED

Matlab, YALMIP and SeDuMi (or any other linear SDP solver)

6. NAMING CONVENTION IN GEO

Most of the data files in GEO are named

t<n1>x<n2>.geo or t<n1>x<n2>ff.geo 

Here n1 and n2 are the numbers of nodes of the ground structure
in the horizontal and vertical direction, respectively. So the truss
"t11x3.geo" will have 33 nodes.

"ff" in the name means that all the nodes are connected by potential bars.
If "ff" is not in the name, then only neighbouring nodes may be connected.


