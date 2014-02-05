###############################################################################
# Eventually, all the global variables should be initialized in this file  
# This file should be sourced at the end of bilayer.tcl                
###############################################################################
array set carbGlyc {
  #center 1
  #ampl 1
  #sigma 1
  target_c 14
  target_a 1
  target_s 3
  tol_c -1
  tol_a -1
  tol_s -1
}

array set phosphate {
  target_c 18
  target_a 1
  target_s 3
  tol_c -1
  tol_a -1
  tol_s -1
}

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
