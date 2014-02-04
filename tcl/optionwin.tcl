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
		#This command sets the appropriate number of electrons to the user input value
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
# Set up a window for additional soft constraints
###############################################################################
proc soft_constraint_window {} {
	global choline phosphate carbGlyc methine methyl
	
	set w 10
	
	toplevel .sfcnst
	wm title .sfcnst "Additional Soft Constraints"

	set f [ frame .sfcnst.f]
	
	label $f.l_Ch -text "Ch"
	label $f.l_Ph -text "Ph"
	label $f.l_CG -text "CG"
	label $f.l_CH -text "c1"
	label $f.l_CH3 -text "c3"
	
	label $f.l_tar_c -text "X"
	label $f.l_tol_c -text "t"
	label $f.l_tar_a -text "C"
	label $f.l_tol_a -text "t"
	label $f.l_tar_s -text "S"
	label $f.l_tol_s -text "t"
	
	entry $f.e_Ch_0 -textvariable choline(target_c) -width $w
	entry $f.e_Ch_1 -textvariable choline(tol_c) -width $w
	entry $f.e_Ch_2 -textvariable choline(target_a) -width $w
	entry $f.e_Ch_3 -textvariable choline(tol_a) -width $w
	entry $f.e_Ch_4 -textvariable choline(target_s) -width $w
	entry $f.e_Ch_5 -textvariable choline(tol_s) -width $w
	
	entry $f.e_Ph_0 -textvariable phosphate(target_c) -width $w
	entry $f.e_Ph_1 -textvariable phosphate(target_a) -width $w
	entry $f.e_Ph_2 -textvariable phosphate(target_s) -width $w
	entry $f.e_Ph_3 -textvariable phosphate(tol_c) -width $w
	entry $f.e_Ph_4 -textvariable phosphate(tol_a) -width $w
	entry $f.e_Ph_5 -textvariable phosphate(tol_s) -width $w

	entry $f.e_CG_0 -textvariable carbGlyc(target_c) -width $w
	entry $f.e_CG_1 -textvariable carbGlyc(target_a) -width $w
	entry $f.e_CG_2 -textvariable carbGlyc(target_s) -width $w
	entry $f.e_CG_3 -textvariable carbGlyc(tol_c) -width $w
	entry $f.e_CG_4 -textvariable carbGlyc(tol_a) -width $w
	entry $f.e_CG_5 -textvariable carbGlyc(tol_s) -width $w

	entry $f.e_CH_0 -textvariable methine(target_c) -width $w
	entry $f.e_CH_1 -textvariable methine(target_a) -width $w
	entry $f.e_CH_2 -textvariable methine(target_s) -width $w
	entry $f.e_CH_3 -textvariable methine(tol_c) -width $w
	entry $f.e_CH_4 -textvariable methine(tol_a) -width $w
	entry $f.e_CH_5 -textvariable methine(tol_s) -width $w

	entry $f.e_CH3_0 -textvariable methyl(target_c) -width $w
	entry $f.e_CH3_1 -textvariable methyl(target_a) -width $w
	entry $f.e_CH3_2 -textvariable methyl(target_s) -width $w
	entry $f.e_CH3_3 -textvariable methyl(tol_c) -width $w
	entry $f.e_CH3_4 -textvariable methyl(tol_a) -width $w
	entry $f.e_CH3_5 -textvariable methyl(tol_s) -width $w

  grid $f -column 0 -row 0
  grid $f.l_tar_c -column 0 -row 1 
  grid $f.l_tol_c -column 0 -row 2 
  grid $f.l_tar_a -column 0 -row 3
  grid $f.l_tol_a -column 0 -row 4
  grid $f.l_tar_s -column 0 -row 5
  grid $f.l_tol_s -column 0 -row 6
  
	grid $f.l_Ch -column 1 -row 0
	grid $f.e_Ch_0 -column 1 -row 1
	grid $f.e_Ch_1 -column 1 -row 2
	grid $f.e_Ch_2 -column 1 -row 3
	grid $f.e_Ch_3 -column 1 -row 4
	grid $f.e_Ch_4 -column 1 -row 5
	grid $f.e_Ch_5 -column 1 -row 6
	
	grid $f.l_Ph -column 2 -row 0
	grid $f.e_Ph_0 -column 2 -row 1
	grid $f.e_Ph_1 -column 2 -row 2
	grid $f.e_Ph_2 -column 2 -row 3
	grid $f.e_Ph_3 -column 2 -row 4
	grid $f.e_Ph_4 -column 2 -row 5
	grid $f.e_Ph_5 -column 2 -row 6

	grid $f.l_CG -column 3 -row 0
	grid $f.e_CG_0 -column 3 -row 1
	grid $f.e_CG_1 -column 3 -row 2
	grid $f.e_CG_2 -column 3 -row 3
	grid $f.e_CG_3 -column 3 -row 4
	grid $f.e_CG_4 -column 3 -row 5
	grid $f.e_CG_5 -column 3 -row 6  
}
