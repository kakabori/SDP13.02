###############################################################################
# Eventually, all the global variables should be initialized in this file  
# This file should be sourced at the end of bilayer.tcl                
###############################################################################

# For link between the tcl x array and C++ x array, see G2model::linkvar() 
# Here, show the variable names with their index
# 0 XCG, 1 CCG, 2 SCG, 3 XPh, 4 CPh, 5 SPh, 6 XCh, 7 CCh, 8 SCh 
# 9 XC, 10 CC, 11 SC, 12 Xc1, 13 Cc1, 14 Sc1, 15 Xc3, 16 Cc3, 17 Sc3
# 18 r, 19 r12, 20 RCG, 21 RPh, 22 Rm, 23 sigR,
# Among the above 24 parameters, 6 are not fitting parameters, so we have
# total of 18 fitting parameters => Pnum = 18 (see bilayer.cpp)

array set upperBounds {}
array set lowerBounds {}
for {set i 0} {$i < $Tnum} {
  set hasLowerBound($i) false
  set hasUpperBound($i) false
}

array set carbGlyc "target_c 14 target_a 1 target_s 3 tol_c -1 tol_a -1 tol_s -1"
array set phosphate "target_c 18 target_a 1 target_s 3 tol_c -1 tol_a -1 tol_s -1"
array set choline {
  target_c 20
  target_a 1
  target_s 3
  tol_c -1
  tol_a -1
  tol_s -1
}
array set methine {
  target_c 5
  target_a 1
  target_s 3
  tol_c -1
  tol_a -1
  tol_s -1
}
array set methyl {
  target_c 0
  target_a 1
  target_s 3
  tol_c -1
  tol_a -1
  tol_s -1
}
