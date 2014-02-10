proc export2g {filename} {
	global x smpfilename fract Tnum chifactorN DBgibbs xpin SNUM status max_peak
  global choline phosphate carbGlyc methine methyl
	global lowerBounds upperBounds hasLowerBound hasUpperBound
	set fid [open $filename w]
	puts $fid "# [clock format [clock second]]"
	puts $fid "# smp_filename: $smpfilename"
	puts $fid "# [spinfo]"
	puts $fid "# [spinfoX]"
	puts $fid "# [spinfoN]"
	puts $fid "# A             =$x(24)"
	puts $fid "# CG(X,C,S)	   =$x(0),	$x(1),	$x(2)"
	puts $fid "# Ph            =$x(3),	$x(4),	$x(5)"
	puts $fid "# Ch            =$x(6),	$x(7),	$x(8)"
	puts $fid "# hydroCarbon   =$x(9),	$x(10),	$x(11)"
	puts $fid "# SigmaC        =$x(11),	$x(56),	$x(57)"
	puts $fid "# DoubleBond    =$x(12),	$x(13),	$x(14)"
	puts $fid "# Methyl        =$x(15),	$x(16),	$x(17)"
	puts $fid "# r(r,t,sr)     =$x(18),	$x(48),	$x(49)"
	puts $fid "# r12(r,t,sr)   =$x(19),	$x(50),	$x(51)"
	puts $fid "# RCG(R,t,sR)   =$x(20),	$x(52),	$x(53)"
	puts $fid "# RPh(R,t,sR)   =$x(21),	$x(54),	$x(55)"
	puts $fid "# Rm		   =$x(22)"
	puts $fid "# sigR		   =$x(23)"
	puts $fid "# n(C2,C1,C3)   =$x(25),	$x(26),	$x(27)"
	puts $fid "# V(C2,C1,C3)   =$x(35),	$x(36),	$x(37)"
	puts $fid "# VL		   =$x(28)"
	puts $fid "# VHL=VHR       =$x(29)"
	puts $fid "# VCG           =$x(32)"
	puts $fid "# VPh, VCh      =$x(33),	$x(34)"
	puts $fid "# VH=VHL+VHpep  =$x(31)"
	puts $fid "# Vpep          =$x(30)"
	puts $fid "# D             =$x(38)"
	puts $fid "# DPP           =$x(39)"
	puts $fid "# DB            =$x(40)"
	puts $fid "# DH1           =$x(42),	$x(43),	$x(44)"
	puts $fid "# DC            =$x(41),	$x(58),	$x(59)"
	puts $fid "# dXH           =$x(45),	$x(46),	$x(47)"
	puts $fid "# dXH2          =$x(60),	$x(61),	$x(62)"
	puts $fid "# dXH3          =$x(67),	$x(68),	$x(69)"
	
#	puts $fid "# nw'*Vw        =[expr (($x(33)+9)*$x(32)-$x(21)-$x(24))]"
#	puts $fid "# DW            =[expr $x(25)-2.0*($x(21)+$x(24))/$x(32)]"
#	puts $fid "# nW*Vw         =[expr ($x(32)*$x(25)/2.0-$x(21)-$x(24))]"
	puts $fid "# fract(2,0,1)  =$fract(2),	$fract(0),	$fract(1)"
  puts $fid "set max_peak $max_peak"
	puts $fid "set chifactorN $chifactorN"
	for {set i 0} {$i<$SNUM} {incr i} {
		if {$status($i)==1} {puts $fid "active $i"}
	}
	
	# Save the x array, which holds the values of model parameters
	for {set i 0} {$i<$Tnum} {incr i} {
		puts $fid "set x($i) $x($i); [expr $xpin($i)?"fix":"free"] $i"
	}

  # Save the additional soft constraints; see optionwin.tcl
  foreach key [array names carbGlyc] {
    puts $fid "set carbGlyc($key) $carbGlyc($key)"
  }
  foreach key [array names phosphate] {
    puts $fid "set phosphate($key) $phosphate($key)"
  }
  foreach key [array names choline] {
    puts $fid "set choline($key) $choline($key)"
  }
  foreach key [array names methine] {
    puts $fid "set methine($key) $methine($key)"
  }
  foreach key [array names methyl] {
    puts $fid "set methyl($key) $methyl($key)"
  }
  
  # Save the elements of bound-related arrays; see optionwin.tcl
  # These are the indices of parameters with bound constraints
  foreach i {0 2 3 5 6 8 12 14 15 17} {
    puts $fid "set upperBounds($i) $upperBounds($i)"
  }
  foreach i {0 2 3 5 6 8 12 14 15 17} {
    puts $fid "set lowerBounds($i) $lowerBounds($i)"
  }	
  foreach i {0 2 3 5 6 8 12 14 15 17} {
    puts $fid "set hasUpperBound($i) $hasUpperBound($i)"
  }
  foreach i {0 2 3 5 6 8 12 14 15 17} {
    puts $fid "set hasLowerBound($i) $hasLowerBound($i)"
  }
  
  close $fid
}
