* Main path - excluding drive letter - to project folder
global name	 "\Users\Student\OneDrive - Nexus365\Replication Project\Practical\Specification Curve\Working Folder"

* Network letter 
global drive 		"C:"			


*** Working Drive ***
global workingDir	"${drive}/${name}"




*==============================================================================*
*                   Load and Save HRS Datasets as Version 11                   *
*==============================================================================*

use "$workingDir\H06A_H.dta", clear 
saveold H06A_H.dta, version (11) replace
   
use "$workingDir\H06A_R.dta", clear 
saveold H06A_R.dta, version (11) replace
   
use "$workingDir\H06B_R.dta", clear 
saveold H06B_R.dta, version (11) replace
   
use "$workingDir\H06C_R.dta", clear 
saveold H06C_R.dta, version (11) replace
   
use "$workingDir\H06D_R.dta", clear 
saveold H06D_R.dta, version (11) replace
   
use "$workingDir\H08D_R.dta", clear 
saveold H08D_R.dta, version (11) replace
   
use "$workingDir\H08PR_R.dta", clear 
saveold H08PR_R.dta, version (11) replace
   
use "$workingDir\H10G_R.dta", clear 
saveold H10G_R.dta, version (11) replace
   
use "$workingDir\H12G_R.dta", clear 
saveold H12G_R.dta, version (11) replace
   
use "$workingDir\H12I_R.dta", clear 
saveold H12I_R.dta, version (11) replace
   
use "$workingDir\H14A_R.dta", clear 
saveold H14A_R.dta, version (11) replace
   
use "$workingDir\H14LB_R.dta", clear 
saveold H14LB_R.dta, version (11) replace
   
use "$workingDir\H16C_R.dta", clear 
saveold H16C_R.dta, version (11) replace
   
use "$workingDir\H16PR_H.dta", clear 
saveold H16PR_H.dta, version (11) replace
   
use "$workingDir\H18D_R.dta", clear 
saveold H18D_R.dta, version (11) replace
   
use "$workingDir\H18PR_R.dta", clear 
saveold H18PR_R.dta, version (11) replace
   
use "$workingDir\H06G_R.dta", clear 
saveold H06G_R.dta, version (11) replace
   
use "$workingDir\H06PR_SB.dta", clear 
saveold H06PR_SB.dta, version (11) replace
   
use "$workingDir\H08G_R.dta", clear 
saveold H08G_R.dta, version (11) replace
   
use "$workingDir\H08Q_H.dta", clear 
saveold H08Q_H.dta, version (11) replace
   
use "$workingDir\H10I_R.dta", clear 
saveold H10I_R.dta, version (11) replace
   
use "$workingDir\H10Q_H.dta", clear 
saveold H10Q_H.dta, version (11) replace
   
use "$workingDir\H12IO_R.dta", clear 
saveold H12IO_R.dta, version (11) replace
   
use "$workingDir\H14B_R.dta", clear 
saveold H14B_R.dta, version (11) replace
   
use "$workingDir\H14PR_H.dta", clear 
saveold H14PR_H.dta, version (11) replace
   
use "$workingDir\H16D_R.dta", clear 
saveold H16D_R.dta, version (11) replace
   
use "$workingDir\H16PR_R.dta", clear 
saveold H16PR_R.dta, version (11) replace
   
use "$workingDir\H18G_R.dta", clear 
saveold H18G_R.dta, version (11) replace
   
use "$workingDir\H18Q_H.dta", clear 
saveold H18Q_H.dta, version (11) replace
   
use "$workingDir\X10C_R.dta", clear 
saveold X10C_R.dta, version (11) replace
   
use "$workingDir\H06I_R.dta", clear 
saveold H06I_R.dta, version (11) replace
   
use "$workingDir\H06Q_H.dta", clear 
saveold H06Q_H.dta, version (11) replace
   
use "$workingDir\H08I_R.dta", clear 
saveold H08I_R.dta, version (11) replace
   
use "$workingDir\H10A_R.dta", clear 
saveold H10A_R.dta, version (11) replace
   
use "$workingDir\H10IO_R.dta", clear 
saveold H10IO_R.dta, version (11) replace
   
use "$workingDir\H12A_R.dta", clear 
saveold H12A_R.dta, version (11) replace
   
use "$workingDir\H12LB_R.dta", clear 
saveold H12LB_R.dta, version (11) replace
   
use "$workingDir\H14C_R.dta", clear 
saveold H14C_R.dta, version (11) replace
   
use "$workingDir\H14PR_R.dta", clear 
saveold H14PR_R.dta, version (11) replace
   
use "$workingDir\H16G_R.dta", clear 
saveold H16G_R.dta, version (11) replace
   
use "$workingDir\H16Q_H.dta", clear 
saveold H16Q_H.dta, version (11) replace
   
use "$workingDir\H18I_R.dta", clear 
saveold H18I_R.dta, version (11) replace
   
use "$workingDir\H06LB_R.dta", clear 
saveold H06LB_R.dta, version (11) replace
   
use "$workingDir\X12C_R.dta", clear 
saveold X12C_R.dta, version (11) replace
   
use "$workingDir\H06IO_R.dta", clear 
saveold H06IO_R.dta, version (11) replace
   
use "$workingDir\H06W_R.dta", clear 
saveold H06W_R.dta, version (11) replace
   
use "$workingDir\H08IO_R.dta", clear 
saveold H08IO_R.dta, version (11) replace
   
use "$workingDir\H10B_R.dta", clear 
saveold H10B_R.dta, version (11) replace
   
use "$workingDir\H10LB_R.dta", clear 
saveold H10LB_R.dta, version (11) replace
   
use "$workingDir\H12B_R.dta", clear 
saveold H12B_R.dta, version (11) replace
   
use "$workingDir\H12PR_H.dta", clear 
saveold H12PR_H.dta, version (11) replace
   
use "$workingDir\H14D_R.dta", clear 
saveold H14D_R.dta, version (11) replace
   
use "$workingDir\H14Q_H.dta", clear 
saveold H14Q_H.dta, version (11) replace
   
use "$workingDir\H16I_R.dta", clear 
saveold H16I_R.dta, version (11) replace
   
use "$workingDir\H18A_R.dta", clear 
saveold H18A_R.dta, version (11) replace
   
use "$workingDir\H18IO_R.dta", clear 
saveold H18IO_R.dta, version (11) replace
   
use "$workingDir\X06C_R.dta", clear 
saveold X06C_R.dta, version (11) replace
   
use "$workingDir\X14C_R.dta", clear 
saveold X14C_R.dta, version (11) replace
   
use "$workingDir\H06LB_R.dta", clear 
saveold H06LB_R.dta, version (11) replace
   
use "$workingDir\H08A_R.dta", clear 
saveold H08A_R.dta, version (11) replace
   
use "$workingDir\H08LB_R.dta", clear 
saveold H08LB_R.dta, version (11) replace
   
use "$workingDir\H10C_R.dta", clear 
saveold H10C_R.dta, version (11) replace
   
use "$workingDir\H10PR_H.dta", clear 
saveold H10PR_H.dta, version (11) replace
   
use "$workingDir\H12C_R.dta", clear 
saveold H12C_R.dta, version (11) replace
   
use "$workingDir\H12PR_R.dta", clear 
saveold H12PR_R.dta, version (11) replace
   
use "$workingDir\H14G_R.dta", clear 
saveold H14G_R.dta, version (11) replace
   
use "$workingDir\H16A_R.dta", clear 
saveold H16A_R.dta, version (11) replace
   
use "$workingDir\H16IO_R", clear 
saveold H16IO_R.dta, version (11) replace
   
use "$workingDir\H18B_R.dta", clear 
saveold H18B_R.dta, version (11) replace
   
use "$workingDir\H18LB_R.dta", clear 
saveold H18LB_R.dta, version (11) replace
   
use "$workingDir\X06PR_R.dta", clear 
saveold X06PR_R.dta, version (11) replace
   
use "$workingDir\X16C_R.dta", clear 
saveold X16C_R.dta, version (11) replace
   
use "$workingDir\H06PR_H.dta", clear 
saveold H06PR_H.dta, version (11) replace
   
use "$workingDir\H08B_R.dta", clear 
saveold H08B_R.dta, version (11) replace
   
use "$workingDir\H08PR_H.dta", clear 
saveold H08PR_H.dta, version (11) replace
   
use "$workingDir\H10D_R.dta", clear 
saveold H10D_R.dta, version (11) replace
   
use "$workingDir\H10PR_R.dta", clear 
saveold H10PR_R.dta, version (11) replace
   
use "$workingDir\H12D_R.dta", clear 
saveold H12D_R.dta, version (11) replace
   
use "$workingDir\H12Q_H.dta", clear 
saveold H12Q_H.dta, version (11) replace
   
use "$workingDir\H14I_R.dta", clear 
saveold H14I_R.dta, version (11) replace
   
use "$workingDir\H16B_R.dta", clear 
saveold H16B_R.dta, version (11) replace
   
use "$workingDir\H16LB_R.dta", clear 
saveold H16LB_R.dta, version (11) replace
   
use "$workingDir\H18C_R.dta", clear 
saveold H18C_R.dta, version (11) replace
   
use "$workingDir\H18PR_H.dta", clear 
saveold H18PR_H.dta, version (11) replace
   
use "$workingDir\X08C_R.dta", clear 
saveold X08C_R.dta, version (11) replace
   
use "$workingDir\H06PR_MC.dta", clear 
saveold H06PR_MC.dta, version (11) replace
   
use "$workingDir\H08C_R.dta", clear 
saveold H08C_R.dta, version (11) replace
   
use "$workingDir\H06PR_R.dta", clear 
saveold H06PR_R.dta, version (11) replace



