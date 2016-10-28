%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 Compute closed null-geodesics in 2D               %
%                                                                   %
%                      Mattia Serra (ETH Zurich)                    %                               
%                   http://www.zfm.ethz.ch/~serra/                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

License:

This software is made public for research use only. It may be modified and redistributed
under the terms of the GNU General Public License. 

Algorithm:

This code implements theoretical results developed by Serra & Haller (2016) [1] at ETH Zurich. 
Specifically, the code developed here can be used to compute closed null-geodesic in 2D flows 
in a fully automated fashion. Examples include elliptic Lagrangian Coherent Structures 
(LCSs) [3], and elliptic Objective Eulerian Coherent Structures (OECSs) [2].

Citation:

Please cite [1] if you use the code in your own work.

----------------------------------------------------------------------------- 
References:

[1] M. Serra, and G. Haller, Efficient Computation of Null-Geodesic 
    with Applications to Coherent Vortex Detection, submitted, (2016). 
	
[2] M. Serra, and G. Haller, Objective Eulerian Coherent Structures,  
    Chaos, 26 (2016) 053110.

[3] G. Haller, FJ. Beron-Vera, Coherent Lagrangian vortices: the black
    holes of turbulence.  J. Fluid Mech. 731 (2013) R4
-----------------------------------------------------------------------------

Tested on Matlab R2015b.

Installation notes :

1) After you unzipped the files to mydir, 
   put the Current Directory in Matlab to mydir

2) In the Matlab command prompt,
   type "add_path" to add the necessary folders to the top of the search path

3) You can run the Main file called: ClosedNullGeodesics.m contained in the folder \Main.

The altimeter products used in this work are produced by SSALTO/DUACS and distributed by AVISO, 
with support from CNES (http://www.aviso.oceanobs.com/duacs). 

NOTE: This code may be inporved and subjected to several changes. Therefore, we suggest to visit this 
page and check if the version you downloaded is up to date.  


Maintained by Mattia Serra,
serram at ethz dot ch
http://www.zfm.ethz.ch/~serra/
October 26, 2016.
