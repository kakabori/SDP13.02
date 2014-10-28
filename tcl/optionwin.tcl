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
	
	toplevel .c
	wm title .c "Peptide Gaussian"
	set f [ frame .c.f ]
	grid $f -column 0 -row 0

	label $f.head -text "peptide in head region (uncheck if in tail region)"
	checkbutton $f.headbutton -variable head 
	grid $f.head -row 0 -column 1 -sticky "w"
	grid $f.headbutton -row 0 -column 0 -sticky "e"
	
	#Set up an entry field to specify Xpep, the center of a peptide Gaussian
	label $f.xlabel -text "Gaussian center"
	entry $f.xentry -textvariable Xpep -width 6 -borderwidth 2 
	grid $f.xlabel -row 1 -column 1 -sticky "w"
	grid $f.xentry -row 1 -column 0
	
	# Set up an entry field to specify Cpep, the amplitude of a peptide Gaussian
	#label $f.clabel -text "Gaussian amplitude (constrained by volume)"
	#entry $f.centry -textvariable Cpep -width 6 -borderwidth 2
	#grid $f.clabel -row 2 -column 0 -sticky "w"
	#grid $f.centry -row 2 -column 1
	
	#Set up an entry field to specify Spep, the width of a peptide Gaussian
	label $f.slabel -text "Gaussian width"
	entry $f.sentry -textvariable Spep -width 6 -borderwidth 2 
	grid $f.slabel -row 3 -column 1 -sticky "w"
	grid $f.sentry -row 3 -column 0
	
	#Set up an entry field to specify Vpep, the volume of the peptide
	label $f.vplabel -text "volume per (fractional) peptide"
	entry $f.vpentry -textvariable Vpep -width 6 -borderwidth 2 
	grid $f.vplabel -row 4 -column 1 -sticky "w"
	grid $f.vpentry -row 4 -column 0
	
	#Set up an entry field to specify Vlipid, the volume of the lipid
	label $f.vllabel -text "volume per lipid"
	entry $f.vlentry -textvariable Vlipid -width 6 -borderwidth 2 
	grid $f.vllabel -row 5 -column 1 -sticky "w"
	grid $f.vlentry -row 5 -column 0
	
	#Set up an entry field to specify VCG, the volume of the CG group
	label $f.vcglabel -text "volume per CG group"
	entry $f.vcgentry -textvariable VCG -width 6 -borderwidth 2 
	grid $f.vcglabel -row 6 -column 1 -sticky "w"
	grid $f.vcgentry -row 6 -column 0
	
	#Set up an entry field to specify VPh, the volume of the Ph group
	label $f.vphlabel -text "volume per Ph group"
	entry $f.vphentry -textvariable VPh -width 6 -borderwidth 2 
	grid $f.vphlabel -row 7 -column 1 -sticky "w"
	grid $f.vphentry -row 7 -column 0
	
	#Set up an entry field to specify VMe, the volume of the methylene group
	label $f.vmelabel -text "volume per methylene group"
	entry $f.vmeentry -textvariable VMe -width 6 -borderwidth 2 
	grid $f.vmelabel -row 8 -column 1 -sticky "w"
	grid $f.vmeentry -row 8 -column 0
	
	#Set up an entry field to specify NePep, the number of electrons in the peptide
	label $f.neplabel -text "number of electrons per (fractional) peptide"
	entry $f.nepentry -textvariable NePep -width 6 -borderwidth 2 
	grid $f.neplabel -row 9 -column 1 -sticky "w"
	grid $f.nepentry -row 9 -column 0
	
	#Set up an Apply button to reflect the parameter setting
	button $f.keep -text Apply -borderwidth 4 -command { setPepGauss $head \
	    $Xpep $Cpep $Spep $Vpep $Vlipid $VCG $VPh $VMe $NePep }
	grid $f.keep -row 10 -column 0 -columnspan 2 -sticky "n"
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
	wm title .c "Additional Constraints"
	set f [ frame .c.f]
	
	frame $f.f_1 -width 10
	frame $f.f_2 -width 20
	
	label $f.l_left -text "soft constraint\nparameters"
	label $f.l_right -text "upper/lower bound parameters"
	
	# Set up labels
	label $f.l_1 -text "target\nvalue"
	label $f.l_2 -text "tolerance"
	label $f.l_3 -text "lower\nbound"
	label $f.l_4 -text "upper\nbound"
	label $f.l_5 -text ""
	label $f.l_6 -text ""
	
    set n XCG
    set i 0
    label $f.l_${n}_0 -text "XCG"
	entry $f.e_${n}_2 -textvariable carbGlyc(target_c) -width $w
	entry $f.e_${n}_3 -textvariable carbGlyc(tol_c) -width $w
	entry $f.e_${n}_5 -textvariable lowerBounds($i) -width $w
	entry $f.e_${n}_6 -textvariable upperBounds($i) -width $w 
	checkbutton $f.c_${n}_7 -variable hasLowerBound($i)
	checkbutton $f.c_${n}_8 -variable hasUpperBound($i)
	
	set n SCG
	set i 2
	label $f.l_${n}_0 -text "SCG"
	entry $f.e_${n}_2 -textvariable carbGlyc(target_s) -width $w
	entry $f.e_${n}_3 -textvariable carbGlyc(tol_s) -width $w
	entry $f.e_${n}_5 -textvariable lowerBounds($i) -width $w 
	entry $f.e_${n}_6 -textvariable upperBounds($i) -width $w
    checkbutton $f.c_${n}_7 -variable hasLowerBound($i) 
	checkbutton $f.c_${n}_8 -variable hasUpperBound($i) 
	
	set n XPh
	set i 3
	label $f.l_XPh_0 -text "XPh"
	entry $f.e_XPh_2 -textvariable phosphate(target_c) -width $w
	entry $f.e_XPh_3 -textvariable phosphate(tol_c) -width $w	
	entry $f.e_${n}_5 -textvariable lowerBounds($i) -width $w 
	entry $f.e_${n}_6 -textvariable upperBounds($i) -width $w
    checkbutton $f.c_${n}_7 -variable hasLowerBound($i) 
	checkbutton $f.c_${n}_8 -variable hasUpperBound($i) 
	
	set n SPh
	set i 5
	label $f.l_SPh_0 -text "SPh"
	entry $f.e_SPh_2 -textvariable phosphate(target_s) -width $w
	entry $f.e_SPh_3 -textvariable phosphate(tol_s) -width $w	
	entry $f.e_${n}_5 -textvariable lowerBounds($i) -width $w 
	entry $f.e_${n}_6 -textvariable upperBounds($i) -width $w
    checkbutton $f.c_${n}_7 -variable hasLowerBound($i)
	checkbutton $f.c_${n}_8 -variable hasUpperBound($i)	

	set n XCh
	set i 6
	label $f.l_XCh_0 -text "XCh"	
	entry $f.e_XCh_2 -textvariable choline(target_c) -width $w
	entry $f.e_XCh_3 -textvariable choline(tol_c) -width $w
	entry $f.e_${n}_5 -textvariable lowerBounds($i) -width $w 
	entry $f.e_${n}_6 -textvariable upperBounds($i) -width $w
    checkbutton $f.c_${n}_7 -variable hasLowerBound($i)
	checkbutton $f.c_${n}_8 -variable hasUpperBound($i)			

	set n SCh
	set i 8
	label $f.l_SCh_0 -text "SCh"	
	entry $f.e_SCh_2 -textvariable choline(target_s) -width $w
	entry $f.e_SCh_3 -textvariable choline(tol_s) -width $w
	entry $f.e_${n}_5 -textvariable lowerBounds($i) -width $w 
	entry $f.e_${n}_6 -textvariable upperBounds($i) -width $w
    checkbutton $f.c_${n}_7 -variable hasLowerBound($i)
	checkbutton $f.c_${n}_8 -variable hasUpperBound($i)
	
	set n Xc1
	set i 12		
	label $f.l_Xc1_0 -text "Xc1"	
	entry $f.e_Xc1_2 -textvariable methine(target_c) -width $w
	entry $f.e_Xc1_3 -textvariable methine(tol_c) -width $w
	entry $f.e_${n}_5 -textvariable lowerBounds($i) -width $w 
	entry $f.e_${n}_6 -textvariable upperBounds($i) -width $w
    checkbutton $f.c_${n}_7 -variable hasLowerBound($i)
	checkbutton $f.c_${n}_8 -variable hasUpperBound($i)

	set n Sc1
	set i 14
	label $f.l_Sc1_0 -text "Sc1"
	entry $f.e_Sc1_2 -textvariable methine(target_s) -width $w	
	entry $f.e_Sc1_3 -textvariable methine(tol_s) -width $w
	entry $f.e_${n}_5 -textvariable lowerBounds($i) -width $w 
	entry $f.e_${n}_6 -textvariable upperBounds($i) -width $w
    checkbutton $f.c_${n}_7 -variable hasLowerBound($i)
	checkbutton $f.c_${n}_8 -variable hasUpperBound($i)	

	set n Xc3
	set i 15
	label $f.l_Xc3_0 -text "Xc3"
	entry $f.e_Xc3_2 -textvariable methyl(target_c) -width $w
	entry $f.e_Xc3_3 -textvariable methyl(tol_c) -width $w
	entry $f.e_${n}_5 -textvariable lowerBounds($i) -width $w 
	entry $f.e_${n}_6 -textvariable upperBounds($i) -width $w
    checkbutton $f.c_${n}_7 -variable hasLowerBound($i)
	checkbutton $f.c_${n}_8 -variable hasUpperBound($i)

	set n Sc3
	set i 17
	label $f.l_Sc3_0 -text "Sc3"
	entry $f.e_Sc3_2 -textvariable methyl(target_s) -width $w
	entry $f.e_Sc3_3 -textvariable methyl(tol_s) -width $w
	entry $f.e_${n}_5 -textvariable lowerBounds($i) -width $w 
	entry $f.e_${n}_6 -textvariable upperBounds($i) -width $w
    checkbutton $f.c_${n}_7 -variable hasLowerBound($i)
	checkbutton $f.c_${n}_8 -variable hasUpperBound($i)

    # Pack widgets
	grid $f -column 0 -row 0
    
    grid $f.f_1 -column 1 -row 0 -rowspan 12
	grid $f.l_left -column 2 -row 0 -columnspan 2 -sticky "n"
	grid $f.f_2 -column 4 -row 0 -rowspan 12
	grid $f.l_right -column 5 -row 0 -columnspan 4 -sticky "n"
	
	grid $f.l_1 -column 2 -row 1
	grid $f.l_2 -column 3 -row 1
	grid $f.l_3 -column 5 -row 1
	grid $f.l_5 -column 6 -row 1 -sticky "w"
	grid $f.l_4 -column 7 -row 1
	grid $f.l_6 -column 8 -row 1 -sticky "w"
    
    # starting row value
    set i 2
    foreach name {XCG SCG XPh SPh XCh SCh Xc1 Sc1 Xc3 Sc3} {
        grid $f.l_${name}_0 -column 0 -row $i
        grid $f.e_${name}_2 -column 2 -row $i
        grid $f.e_${name}_3 -column 3 -row $i
        grid $f.e_${name}_5 -column 5 -row $i
        grid $f.c_${name}_7 -column 6 -row $i -sticky "w"
        grid $f.e_${name}_6 -column 7 -row $i      
        grid $f.c_${name}_8 -column 8 -row $i -sticky "w"
        incr i
    }
  
}
