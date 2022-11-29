Date: 
Nov. 28 2022

Description:
This project achieve an open-loop control of flexible needle base on FEA method. This repository include all the source code and experiment data used in the project. The FBG is not included in our project, and the source code related to it is not completed. Simply set the "FBG_switch" to 0 then the codejust ignore all FBG functions. 

Content in Folders:
#-- Data Postprocess --#
"experiment_data_postprocess": include experiment data, photo, analysis code and outputs which generated from those data and code. 
"image_analysis": include source code used in data postprocess.

#-- FBG Related --#
"write_to_file": include source code to save data to files, this is prepared for FBG validation. In our project this is not used.
"rawdata_process": include source code and data that could generate calibration matrix. Also have the function to convert FBG reading into curvatures.
"sm130_interrogator_matlab": include method to connect and get signal from sm130 interrogator.

#-- FEM Related --#
"needle_FEM_realtime": include FEA method and funtion integrating FBG and FEM together to achieve realtime shape sensing. The function has "FBG_switch" which can switch on/off the FBG.

#-- flexible needle control --#
"needle_control": include souce code to control the needle tip to desire position.

#-- visualization --#
"visualization": include function that can achieve parallel visualization of data (position, orientation, curvature and shape). This is complished by "memmapfile", run the function in another matlab process. If the function give error about "memmapfile", the possible solution is delete the "communication.dat" at /tmp.

#-- implement script --#
"test_script": include files below

""" files """: 
"test_data_process.m": test the function of process of raw data that convert FBG singal into curvatures.
"test_FBG_FEM_control": this script test the function that integrating FBG, FEM and control.
"test_FBG_FEM_realtime": this script test the function that integrating FBG and FEM.
"test_FBG_integrated_FEM_Ogden_Utru.m", "test_FEM_returned_curvature.m": these script test the functions that use FBG and FEM.
"test_FBG_FEM_unintegrated.m": this script seperatly test the functions of FBG and FEM.
"test_interrogator.m": this function test the use of interrogator.



