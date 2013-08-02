# This script takes care of header for SIMtoEXP form factor file (*.nff)

set StE_header "##################################################################################
## experimental scattering form factors obtained for DOPC+cholesterol bilayers	##
## Kucerka et al., The European Physical Journal E 23 (2007) 247-254		##
##################################################################################"

set StE_header "##################################################################################
## experimental neutron scattering form factors obtained for ULV data		##
## Kucerka et al., Biophysical Journal 97 (2009) 1926-1932			##
##################################################################################"

set StE_header "##################################################################################
## experimental scattering form factors obtained for ULV data			##
## Kucerka et al., Biochim Biophys Acta - Biomembranes 1808 (2011) 2761-2771	##
##################################################################################"

set StE_header "##################################################################################
## experimental scattering form factors obtained for ULV data			##
## Pan et al., Biochim Biophys Acta - Biomembranes (2012)			##
##################################################################################"

set StE_header ""

set StE_contrast ""
proc SIMtoEXP_contrast {mode contrast rw} {
	global StE_contrast
	set bH -0.374E-4
	set bD  0.667E-4
	set b_contrast [expr $contrast/100.*$bD+(1-$contrast/100.)*$bH]
	set eH 1.
	set eD 1.
	set e_contrast [expr $contrast/100.*$eD+(1-$contrast/100.)*$eH]
	set StE_contrast "##
"
	append StE_contrast "## redefinition of scattering power for D such that water corresponds to "
	append StE_contrast $contrast
	append StE_contrast "% D2O
"
	append StE_contrast "#SLwin
"
	if {$mode=="x"} {
		append StE_contrast "#set e(1) "
		append StE_contrast $e_contrast
		append StE_contrast "
#setRHO_wat X "
	} elseif {$mode=="n"} {
		append StE_contrast "#set b(1) "
		append StE_contrast $b_contrast
		append StE_contrast "
#setRHO_wat N "
	}
	append StE_contrast $rw
	append StE_contrast "
#updateSIM
"
	append StE_contrast "##"

#	set StE_contrast ""
}

proc EXPtoSMP {fin fout} {
	set fidin [open $fin r]
	set fidout [open $fout w]
	set l ""
	while { [gets $fidin l] > 0 } {
		if {[string is double [lindex $l 0]]==1} {
			set q [lindex $l 0]
			set F [lindex $l 1]
			set F2 [expr $F*$F]
			set err [lindex $l 2]
			puts $fidout "$F2 $q $err \\"
		}
	}
	close $fidout
	close $fidin
}
