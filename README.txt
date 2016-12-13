This code provides a minimal working example (MWE) of the paper
"Inverting RANSAC: Global Model Detection via Inlier Rate Estimation"
presented at CVPR 2015

---	Version 0.1 ---

Please note:
- The code is provided "as-is" and my not work on your machine - please let us know if this happens. The code was tested on matlab 2013b.
- Some functions are implemented in c++. We provide compiled mex files for windows 64-bit. 
- Some matlab functions were compiled using "matlab coder". We provide compiled mex files for windows 64-bit. Aditionally we added the "prj" file and the "codegen" folder.
- Run the 'compile-mex' script
- To run the MWE - try "script_MWE_ComparisonToUSAC.m" in the "code" folder.
- The MWE code relies on the following 3rd party packages (not provided, will be added later):
	1) USAC
	2) Mikolajczyk's data-set
	3) vl_feat
	
	
Good luck!
The authors