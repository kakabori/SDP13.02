###############################################################################
set Xpep 16
set Cpep 1
set Spep 3
set Vpep 48.1
set Vlipid 1313.5
set VCG 153.7
set VPh 180.5
set VMe 27.2
set NePep 21.487
###############################################################################
# Set up a window to set the peptide Gaussian parameters.
proc pepOptWin {} {
	global head Xpep Cpep Spep Vpep Vlipid VCG VPh VMe NePep
	
	toplevel .fitopts
	wm title .fitopts "Peptide Gaussian"

	#############
	set optframe [ frame .fitopts.optframe ]

	label $optframe.head -text "peptide in head region (uncheck if in tail region)"
	checkbutton $optframe.headbutton -variable head
	grid $optframe.headbutton -row 1 -column 1
	grid $optframe.head -row 1 -column 2

	grid $optframe -row 1
	
	#############
	set tframe [ frame .fitopts.tframe ]
	
	#Set up an entry field to specify Xpep, the center of a peptide Gaussian
	label $tframe.xlabel -text "Xpep"
	entry $tframe.xentry -textvariable Xpep -width 6 -borderwidth 2 
	grid $tframe.xlabel -row 1 -column 1
	grid $tframe.xentry -row 1 -column 2
	
	#Set up an entry field to specify Cpep, the amplitude of a peptide Gaussian
	label $tframe.clabel -text "Cpep"
	entry $tframe.centry -textvariable Cpep -width 6 -borderwidth 2
	grid $tframe.clabel -row 2 -column 1
	grid $tframe.centry -row 2 -column 2
	
	#Set up an entry field to specify Spep, the width of a peptide Gaussian
	label $tframe.slabel -text "Spep"
	entry $tframe.sentry -textvariable Spep -width 6 -borderwidth 2 
	grid $tframe.slabel -row 3 -column 1
	grid $tframe.sentry -row 3 -column 2
	
	#Set up an entry field to specify Vpep, the volume of the peptide
	label $tframe.vplabel -text "Vpep"
	entry $tframe.vpentry -textvariable Vpep -width 6 -borderwidth 2 
	grid $tframe.vplabel -row 4 -column 1
	grid $tframe.vpentry -row 4 -column 2
	
	#Set up an entry field to specify Vlipid, the volume of the lipid
	label $tframe.vllabel -text "Vlipid"
	entry $tframe.vlentry -textvariable Vlipid -width 6 -borderwidth 2 
	grid $tframe.vllabel -row 5 -column 1
	grid $tframe.vlentry -row 5 -column 2
	
	#Set up an entry field to specify VCG, the volume of the CG group
	label $tframe.vcglabel -text "VCG"
	entry $tframe.vcgentry -textvariable VCG -width 6 -borderwidth 2 
	grid $tframe.vcglabel -row 6 -column 1
	grid $tframe.vcgentry -row 6 -column 2
	
	#Set up an entry field to specify VPh, the volume of the Ph group
	label $tframe.vphlabel -text "VPh"
	entry $tframe.vphentry -textvariable VPh -width 6 -borderwidth 2 
	grid $tframe.vphlabel -row 7 -column 1
	grid $tframe.vphentry -row 7 -column 2
	
	#Set up an entry field to specify VMe, the volume of the methylene group
	label $tframe.vmelabel -text "VMe"
	entry $tframe.vmeentry -textvariable VMe -width 6 -borderwidth 2 
	grid $tframe.vmelabel -row 8 -column 1
	grid $tframe.vmeentry -row 8 -column 2
	
	grid $tframe -row 2
	#############
	set nframe [ frame .fitopts.nframe ]
	
	#Set up an entry field to specify NePep, the number of electrons in the peptide
	label $nframe.neplabel -text "Number of electrons in peptide"
	entry $nframe.nepentry -textvariable NePep -width 6 -borderwidth 2 
	grid $nframe.neplabel -row 1 -column 1
	grid $nframe.nepentry -row 1 -column 2
	
	grid $nframe -row 3
	#############
	
	#Set up an Apply button to reflect the parameter setting
	set bframe [ frame .fitopts.bframe ]
	button $bframe.keep -text Apply -borderwidth 4 -command { setPepGauss $head \
	  $Xpep $Cpep $Spep $Vpep $Vlipid $VCG $VPh $VMe $NePep }
	grid $bframe.keep -row 1 -column 1
	grid $bframe -row 4
}

#Set peptide's Gaussian parameters, relevant volumes and volume ratio in appropriate entries.
#Also, set the number of electrons in the choline or double-bond group to that in the peptide.
proc setPepGauss {head Xpep Cpep Spep Vpep Vlipid VCG VPh VMe NePep} {
	global x
	if {$head == 1} {
		#peptide is in head region; using choline group as a peptide.
		#Fix the double-bond parameters to zero. Vhead is the head group volume.
		set Vhead [expr $Vpep + $VCG + $VPh]
		#choline Gaussian parameters
		set x(6) $Xpep
		set x(7) $Cpep
		set x(8) $Spep
		#Methine (double bond) Gaussian parameters
		set x(12) 0
		set x(13) 0
		set x(14) 0
		#r12 parameter, the ratio, methine/methylene
		set x(19) 0
		set x(51) 0
		#RCG parameter
		set x(20) [expr double($VCG) / $Vhead]
		set x(53) [expr double($VCG) / $Vhead]
		#RPh parameter
		set x(21) [expr double($VPh) / $Vhead]
		set x(55) [expr double($VPh) / $Vhead]
		#number of methine needs to be set to 0
		set x(26) 0
		#Total volume = lipid volume + peptide volume
		set x(28) [expr $Vpep + $Vlipid]
		#Head group volume
		set x(29) $Vhead
		free 6
		free 8
		#Methine related parameters are fixed, not free
		fix 12
		fix 14
		fix 19
		#This command sets the number of electrons to the user input value
		setNePep $head $NePep
	} else {
		#peptide is in tail region
		set Vhead [expr $VCG + $VPh]
		#In this case, peptide is assumed to be in both chains. So, the peptide volume 
		#per chain is half of its value when peptide is in the head region.
		set VpepPerChain [expr $Vpep / 2.0]
		set x(6) 0
		set x(7) 0
		set x(8) 0
		set x(12) $Xpep
		set x(13) $Cpep
		set x(14) $Spep
		set x(19) [expr double($VpepPerChain) / $VMe]
		set x(51) [expr double($VpepPerChain) / $VMe]
		set x(20) [expr double($VCG) / $Vhead]
		set x(53) [expr double($VCG) / $Vhead]
		set x(21) [expr double($VPh) / $Vhead]
		set x(55) [expr double($VPh) / $Vhead]
		#number of methine needs to be set to 1
		set x(26) 1
		set x(28) [expr $Vpep + $Vlipid]
		set x(29) $Vhead
		#Choline Gaussian parameters are fixed
		fix 6
		fix 8
		free 12
		free 14
		free 19
		#Also, number of electrons per peptide per chain is cut in half
		set NePep [expr $NePep / 2.0]
		setNePep $head $NePep
	}
	
}
###############################################################################
set rXCG 0.1
set rSCG 0.1
set rXPh 0.1
set rSPh 0.1
set rXCh 0.1
set rSCh 0.1
set rXC 0.1
set rSC 0.1
set rXc1 0.1
set rSc1 0.1
set rSc3 0.1
###############################################################################
# Set up a window to set randomization parameters.
proc randOptWin {} {
	global rXCG rSCG rXPh rSPh rXCh rSCh rXC rSC rXc1 rSc1 rSc3 rr rr12 rRCG rRPh
	global random repeat
	
	toplevel .randopts
	wm title .randopts "Randomization Parameters"

	set tframe [ frame .randopts.tframe ]
	
	#Set up an entry field to specify Xpep, the center of a peptide Gaussian
	label $tframe.rlabel -text "random"
	entry $tframe.rentry -textvariable random -width 4
	pack $tframe.rlabel $tframe.rentry
	
	#Set up an entry field to specify Cpep, the amplitude of a peptide Gaussian
	label $tframe.rplabel -text "repeat"
	entry $tframe.rpentry -textvariable repeat -width 4
	pack $tframe.rplabel $tframe.rpentry
	
	#Set up an entry field to specify Xpep, the center of a peptide Gaussian
	label $tframe.r0label -text "rXCG"
	entry $tframe.r0entry -textvariable rXCG -width 4
	pack $tframe.r0label $tframe.r0entry
	
	#Set up an entry field to specify Cpep, the amplitude of a peptide Gaussian
	label $tframe.r2label -text "rSCG"
	entry $tframe.r2entry -textvariable rSCG -width 4
	pack $tframe.r2label $tframe.r2entry
	
	#Set up an entry field to specify Spep, the width of a peptide Gaussian
	label $tframe.r3label -text "rXPh"
	entry $tframe.r3entry -textvariable rXPh -width 4
	pack $tframe.r3label $tframe.r3entry
	
	#Set up an entry field to specify Vpep, the volume of the peptide
	label $tframe.r5label -text "rSPh"
	entry $tframe.r5entry -textvariable rSPh -width 4
	pack $tframe.r5label $tframe.r5entry
	
	#Set up an entry field to specify Vlipid, the volume of the lipid
	label $tframe.r6label -text "rXCh"
	entry $tframe.r6entry -textvariable rXCh -width 4
	pack $tframe.r6label $tframe.r6entry
	
	#Set up an entry field to specify VCG, the volume of the CG group
	label $tframe.r8label -text "rSCh"
	entry $tframe.r8entry -textvariable rSCh -width 4
	pack $tframe.r8label $tframe.r8entry
	
	#Set up an entry field to specify VPh, the volume of the Ph group
	label $tframe.r9label -text "rXC"
	entry $tframe.r9entry -textvariable rXC -width 4
	pack $tframe.r9label $tframe.r9entry
	
	#Set up an entry field to specify VMe, the volume of the methylene group
	label $tframe.r11label -text "rSC"
	entry $tframe.r11entry -textvariable rSC -width 4
	pack $tframe.r11label $tframe.r11entry
	
	#Set up an entry field to specify NePep, the number of electrons in the peptide
	label $tframe.r12label -text "rXc1"
	entry $tframe.r12entry -textvariable rXc1 -width 4
	pack $tframe.r12label $tframe.r12entry
	
	#Set up an entry field to specify NePep, the number of electrons in the peptide
	label $tframe.r14label -text "rSc1"
	entry $tframe.r14entry -textvariable rSc1 -width 4
	pack $tframe.r14label $tframe.r14entry
	
	#Set up an entry field to specify NePep, the number of electrons in the peptide
	label $tframe.r17label -text "rSc3"
	entry $tframe.r17entry -textvariable rSc3 -width 4
	pack $tframe.r17label $tframe.r17entry
	
	grid $tframe -row 0
}
#########################################################################################
proc exportOpts {filename} {
	global head Xpep Cpep Spep Vpep Vlipid VCG VPh VMe NePep
	global rXCG rSCG rXPh rSPh rXCh rSCh rXC rSC rXc1 rSc1 rSc3 rr rr12 rRCG rRPh
	global random repeat
	
	set fid [open $filename w]

	#peptide options
	puts $fid "set head $head"
	puts $fid "set Xpep $Xpep"
	puts $fid "set Cpep $Cpep"
	puts $fid "set Spep $Spep"
	puts $fid "set Vpep $Vpep"
	puts $fid "set Vlipid $Vlipid"
	puts $fid "set VCG $VCG"
	puts $fid "set VPh $VPh"
	puts $fid "set VMe $VMe"
	puts $fid "set NePep $NePep"
	
	#randomization options
	puts $fid "set rXCG $rXCG"
	puts $fid "set rSCG $rSCG"
	puts $fid "set rXPh $rXPh"
	puts $fid "set rSPh $rSPh"
	puts $fid "set rXCh $rXCh"
	puts $fid "set rSCh $rSCh"
	puts $fid "set rXC $rXC"
	puts $fid "set rSC $rSC"
	puts $fid "set rXc1 $rXc1"
	puts $fid "set rSc1 $rSc1"
	puts $fid "set rSc3 $rSc3"
	
	close $fid
}

#########################################################################################
proc go {direc} {
	switch $direc {
		"dopc" {
			cd ~/WinE/chess12/tat_dopc
		}
		"dope" {
			cd ~/WinE/chess12/tat_dopc_dope
		}
		"dops" {
			cd ~/WinE/chess12/tat_dopc_dops
		}
		"nuclear" {
			cd ~/WinE/chess12/tat_nuclear
		}
		"LB" {
			cd ~/WinE/chess12/LB
		}
		default {
			puts "'go' followed by 'dopc', 'dope', 'dops', or 'nuclear'"
		}
	}
}


###############################################################################
# Set up a window for additional soft constraints and upper/lower bounds
###############################################################################
proc constraints_window {} {
	global choline phosphate carbGlyc methine methyl
	global lowerBounds upperBounds hasLowerBound hasUpperBound
	
	# The width of entry fields
	set w 7 
	
	toplevel .c
	wm title .c "Additional Soft Constraints"
	set f [ frame .c.f]
	
	label $f.l_1 -text "target\nvalue"
	label $f.l_2 -text "tolerance"
	label $f.l_3 -text "lower\nbound"
	label $f.l_4 -text "upper\nbound"
	label $f.l_5 -text "lower\nbound on"
	label $f.l_6 -text "upper\nbound on"
	
  set n XCG
  label $f.l_${n}_0 -text "XCG"
	entry $f.e_${n}_1 -textvariable carbGlyc(target_c) -width $w
	entry $f.e_${n}_2 -textvariable carbGlyc(tol_c) -width $w
	entry $f.e_${n}_3 -textvariable lowerBounds(0) -width $w
	entry $f.e_${n}_4 -textvariable upperBounds(0) -width $w 
	checkbutton $f.c_${n}_5 -variable hasLowerBound(0) -width $w
	checkbutton $f.c_${n}_6 -variable hasUpperBound(0) -width $w
	
	set n SCG
	label $f.l_${n}_0 -text "SCG"
	entry $f.e_${n}_1 -textvariable carbGlyc(target_s) -width $w
	entry $f.e_SCG_2 -textvariable carbGlyc(tol_s) -width $w
	entry $f.e_SCG_3 -textvariable lowerBounds(2) -width $w 
	entry $f.e_SCG_4 -textvariable upperBounds(2) -width $w
  checkbutton $f.c_SCG_5 -variable hasLowerBound(2) -width $w
	checkbutton $f.c_SCG_6 -variable hasUpperBound(2) -width $w	
	
	set n XPh
	set i 3
	label $f.l_XPh_0 -text "XPh"
	entry $f.e_XPh_1 -textvariable phosphate(target_c) -width $w
	entry $f.e_XPh_2 -textvariable phosphate(tol_c) -width $w	
	entry $f.e_XPh_3 -textvariable lowerBounds(3) -width $w 
	entry $f.e_XPh_4 -textvariable upperBounds(3) -width $w
  checkbutton $f.c_${n}_5 -variable hasLowerBound($i) -width $w
	checkbutton $f.c_${n}_6 -variable hasUpperBound($i) -width $w		
	
	set n SPh
	set i 5
	label $f.l_SPh_0 -text "SPh"
	entry $f.e_SPh_1 -textvariable phosphate(target_s) -width $w
	entry $f.e_SPh_2 -textvariable phosphate(tol_s) -width $w	
	entry $f.e_SPh_3 -textvariable lowerBounds(5) -width $w 
	entry $f.e_SPh_4 -textvariable upperBounds(5) -width $w
  checkbutton $f.c_${n}_5 -variable hasLowerBound($i) -width $w
	checkbutton $f.c_${n}_6 -variable hasUpperBound($i) -width $w		

	set n XCh
	set i 6
	label $f.l_XCh_0 -text "XCh"	
	entry $f.e_XCh_1 -textvariable choline(target_c) -width $w
	entry $f.e_XCh_2 -textvariable choline(tol_c) -width $w
	entry $f.e_XCh_3 -textvariable lowerBounds(6) -width $w
	entry $f.e_XCh_4 -textvariable upperBounds(6) -width $w
  checkbutton $f.c_${n}_5 -variable hasLowerBound($i) -width $w
	checkbutton $f.c_${n}_6 -variable hasUpperBound($i) -width $w			

	set n SCh
	set i 8
	label $f.l_SCh_0 -text "SCh"	
	entry $f.e_SCh_1 -textvariable choline(target_s) -width $w
	entry $f.e_SCh_2 -textvariable choline(tol_s) -width $w
	entry $f.e_SCh_3 -textvariable lowerBounds(8) -width $w 
	entry $f.e_SCh_4 -textvariable upperBounds(8) -width $w	
  checkbutton $f.c_${n}_5 -variable hasLowerBound($i) -width $w
	checkbutton $f.c_${n}_6 -variable hasUpperBound($i) -width $w
	
	set n Xc1
	set i 12		
	label $f.l_Xc1_0 -text "Xc1"	
	entry $f.e_Xc1_1 -textvariable methine(target_c) -width $w
	entry $f.e_Xc1_2 -textvariable methine(tol_c) -width $w
	entry $f.e_Xc1_3 -textvariable lowerBounds(12) -width $w 
	entry $f.e_Xc1_4 -textvariable upperBounds(12) -width $w 
  checkbutton $f.c_${n}_5 -variable hasLowerBound($i) -width $w
	checkbutton $f.c_${n}_6 -variable hasUpperBound($i) -width $w

	set n Sc1
	set i 14
	label $f.l_Sc1_0 -text "Sc1"
	entry $f.e_Sc1_1 -textvariable methine(target_s) -width $w	
	entry $f.e_Sc1_2 -textvariable methine(tol_s) -width $w
	entry $f.e_Sc1_3 -textvariable lowerBounds(14) -width $w 
	entry $f.e_Sc1_4 -textvariable upperBounds(14) -width $w
  checkbutton $f.c_${n}_5 -variable hasLowerBound($i) -width $w
	checkbutton $f.c_${n}_6 -variable hasUpperBound($i) -width $w 

	set n Xc3
	set i 15
	label $f.l_Xc3_0 -text "Xc3"
	entry $f.e_Xc3_1 -textvariable methyl(target_c) -width $w
	entry $f.e_Xc3_2 -textvariable methyl(tol_c) -width $w
	entry $f.e_Xc3_3 -textvariable lowerBounds(15) -width $w
	entry $f.e_Xc3_4 -textvariable upperBounds(15) -width $w		
  checkbutton $f.c_${n}_5 -variable hasLowerBound($i) -width $w
	checkbutton $f.c_${n}_6 -variable hasUpperBound($i) -width $w

	set n Sc3
	set i 17
	label $f.l_Sc3_0 -text "Sc3"
	entry $f.e_Sc3_1 -textvariable methyl(target_s) -width $w
	entry $f.e_Sc3_2 -textvariable methyl(tol_s) -width $w
	entry $f.e_Sc3_3 -textvariable lowerBounds(17) -width $w
	entry $f.e_Sc3_4 -textvariable upperBounds(17) -width $w
  checkbutton $f.c_${n}_5 -variable hasLowerBound($i) -width $w
	checkbutton $f.c_${n}_6 -variable hasUpperBound($i) -width $w

	grid $f -column 0 -row 0
	
	grid $f.l_1 -column 1 -row 0
	grid $f.l_2 -column 2 -row 0
	grid $f.l_3 -column 3 -row 0
	grid $f.l_4 -column 4 -row 0
	grid $f.l_5 -column 5 -row 0
	grid $f.l_6 -column 6 -row 0
  
  set i 1
  foreach name {XCG SCG XPh SPh XCh SCh Xc1 Sc1 Xc3 Sc3} {
    grid $f.l_${name}_0 -column 0 -row $i
    for {set j 1} {$j < 5} {incr j} {
	    grid $f.e_${name}_$j -column $j -row $i
    }
    puts $name
    grid $f.c_${name}_5 -column 5 -row $i
    grid $f.c_${name}_6 -column 6 -row $i
    incr i
  }
  
}
