# This script takes file with scaling factors information (frm.dat)
# and converts it to the Df input file (*.smp).
# 3rd column is error input by user.
# It writes some default values of parameters,
# that should be changed to the actual ones.
# It is necessary to erase backslash in the last line!

#fin can take a string containing multiple .dat file names, each separated by a space. 
#Use: convert {file1.dat file2.dat file3.dat} files.smp 0.5
#Use: convert {file1.dat file2.dat file3.dat} files.smp ka
proc convert {fin fout err} {
    set fins [split $fin " "]
    ##set fidin [open $fin r]
    set fidout [open $fout w]
    set l ""
    set counter 1
    ## Begin Header ##
    puts $fidout "set direct_err 1"
    puts $fidout "set stepsize_integral 0.05"
    puts $fidout "\n"
    ## End Header
    foreach fi $fins {
		set fidin [open $fi r]
		set name [file rootname $fi]
		puts $fidout "\n"   
		puts $fidout "samplist $counter $name"
		puts $fidout "parameter $counter nobeam \\"
		puts $fidout "1.1775 2.3 5 0 10 1 -63.2 65 0.0\\"
		puts $fidout "x 0.333 67.0 97.0 21.487 7.875 0.0 9.0 0.0\\"
		while { [gets $fidin l] >= 0 } {
			set F [lindex $l 1]
			set q [lindex $l 3]
			if { $err == "ka" } {
				set err_a [lindex $l 4]
				puts $fidout "$F $q $err_a \\"
			} else {puts $fidout "$F $q $err \\"}
			#if { $F<0 } { set F 0 }
    	}
   	 	set counter [expr $counter + 1]
    	close $fidin
 	}
    close $fidout
}

proc moderr {fin sampl qlow qhigh newerr} {
    set fidin [open $fin r]
    set ftmp [open tmp w]
    set l ""
  while { [gets $fidin l] > 0 } {
  puts $ftmp "$l"
  while { [gets $fidin l] > 0 } {
    if { [lindex $l 0]=="parameter" && [lindex $l 1]==$sampl } {
      puts $ftmp "$l"
      gets $fidin l
      puts $ftmp "$l"
      gets $fidin l
      puts $ftmp "$l"
      while { [gets $fidin l] > 0 } {
	 set F [lindex $l 0]
	 set q [lindex $l 1]
	 set err [lindex $l 2]
	 if {$q>=$qlow && $q<$qhigh} {set err $newerr}
	 puts $ftmp "$F $q $err \\"
     }
    }
    puts $ftmp "$l"
  }
  puts $ftmp "$l"
  }
  close $ftmp
  close $fidin
  file rename -force tmp $fin   
}

proc multierr {fin sampl coef} {
    set fidin [open $fin r]
    set ftmp [open tmp w]
    set l ""
  while { [gets $fidin l] > 0 } {
  puts $ftmp "$l"
  while { [gets $fidin l] > 0 } {
    if { [lindex $l 0]=="parameter" && [lindex $l 1]==$sampl } {
      puts $ftmp "$l"
      gets $fidin l
      puts $ftmp "$l"
      gets $fidin l
      puts $ftmp "$l"
      while { [gets $fidin l] > 0 } {
	 set F [lindex $l 0]
	 set q [lindex $l 1]
	 set err [expr [lindex $l 2]*$coef]
	 puts $ftmp "$F $q $err \\"
     }
    }
    puts $ftmp "$l"
  }
  puts $ftmp "$l"
  }
  close $ftmp
  close $fidin
  file rename -force tmp $fin   
}

proc modint {fin sampl qlow qhigh newint} {
    set fidin [open $fin r]
    set ftmp [open tmp w]
    set l ""
  while { [gets $fidin l] > 0 } {
  puts $ftmp "$l"
  while { [gets $fidin l] > 0 } {
    if { [lindex $l 0]=="parameter" && [lindex $l 1]==$sampl } {
      puts $ftmp "$l"
      gets $fidin l
      puts $ftmp "$l"
      gets $fidin l
      puts $ftmp "$l"
      while { [gets $fidin l] > 0 } {
	 set F [lindex $l 0]
	 set q [lindex $l 1]
	 set err [lindex $l 2]
         if {$q>=$qlow && $q<$qhigh} {set F $newint}
         puts $ftmp "$F $q $err \\"
      }
    }
    puts $ftmp "$l"
  }
  puts $ftmp "$l"
  }
  close $ftmp
  close $fidin
  file rename -force tmp $fin   
}

proc delint {fin sampl qlow qhigh} {
    set fidin [open $fin r]
    set ftmp [open tmp w]
    set l ""
  while { [gets $fidin l] > 0 } {
  puts $ftmp "$l"
  while { [gets $fidin l] > 0 } {
    if { [lindex $l 0]=="parameter" && [lindex $l 1]==$sampl } {
      puts $ftmp "$l"
      gets $fidin l
      puts $ftmp "$l"
      gets $fidin l
      puts $ftmp "$l"
      while { [gets $fidin l] > 0 } {
	 set F [lindex $l 0]
	 set q [lindex $l 1]
	 set err [lindex $l 2]
         if {$q<$qlow || $q>=$qhigh}  {puts $ftmp "$F $q $err \\"}
      }
    }
    puts $ftmp "$l"
  }
  puts $ftmp "$l"
  }  
  close $ftmp
  close $fidin
  file rename -force tmp $fin   
}

proc bin {fin sampl n} {
    set fidin [open $fin r]
    set ftmp [open tmp w]
    set l ""
    set F 0.1
    set q 0.1
    set err 0.1
    set end 0
  while { [gets $fidin l] > 0 } {
  puts $ftmp "$l"
  while { [gets $fidin l] > 0 } {
    if { [lindex $l 0]=="parameter" && [lindex $l 1]==$sampl } {
      puts $ftmp "$l"
      gets $fidin l
      puts $ftmp "$l"
      gets $fidin l
      puts $ftmp "$l"
      while { [gets $fidin l] > 0 && $end==0 } {
	    set j 1
	    set F [lindex $l 0]
	    set q [lindex $l 1]
	    set err [lindex $l 2]
	    set sumq $q
	    set sumF $F
	    set sumF2 [expr $F*$F]
	    for {set i 1} {$i<$n} {incr i} {
		if { [gets $fidin l] > 0 } {
		    set F [lindex $l 0]
		    set q [lindex $l 1]
		    set sumq [expr $sumq+$q]
		    set sumF [expr $sumF+$F]
		    set sumF2 [expr $sumF2+$F*$F]
		    incr j
		} else {
		     set i $n
		     set end 1
		}
	    }
	    if {$j==1} {puts $ftmp "$F $q $err \\"} else {
		set F [expr ($sumF+0.0)/$j]
		set q [expr ($sumq+0.0)/$j]
		if {[expr $sumF2/($j+0.0)]>[expr $F*$F]} {set err [expr sqrt(($sumF2/($j+0.0)-$F*$F)/($j-1.0))]} else {set err 0.000001}
		puts $ftmp "$F $q $err \\"
	    }
	    if {$end==1} {puts $ftmp ""}
      }
    }
    puts $ftmp "$l"
  }
  puts $ftmp "$l"
  }
    close $ftmp
    close $fidin
    file rename -force tmp $fin   
}

proc convbin {fin fout n} {
    set fidin [open $fin r]
    set fidout [open $fout w]
    set l ""
    puts $fidout "# This is header for dopcb46 as was used in YL's thesiss" 
    puts $fidout "# Set your own parameters and erase last backslash!"  
    puts $fidout "set direct_err 1"
    puts $fidout "set stepsize_integral 0.05"
    puts $fidout "samplist 0 void"
    puts $fidout "samplist 1 dopcb46"
    puts $fidout "parameter 1 nobeam \\"
    puts $fidout "1.1808 2.5 4 1 10 1 -63.2 48 1.075 \\"
    puts $fidout "x 0.333 67.0 97.0 7.875 9.0 0.0 \\"
    set F 0.1
    set q 0.1
    set err 0.1
    while { [gets $fidin l] >= 0 } {
	    set j 1
	    set F [lindex $l 1]
	    set q [lindex $l 3]
	    set sumq $q
	    set sumF $F
	    set sumF2 [expr $F*$F]
	    for {set i 1} {$i<$n} {incr i} {
		if { [gets $fidin l] > 0 } {
		    set F [lindex $l 1]
		    set q [lindex $l 3]
		    set sumq [expr $sumq+$q]
		    set sumF [expr $sumF+$F]
		    set sumF2 [expr $sumF2+$F*$F]
		    incr j
		}
	    }
	    if {$j==1} {
		#if { $F<0 } { set F 0 }
		puts $fidout "$F $q $err \\"
	    } else {
		set F [expr $sumF/($j+0.0)]
		set q [expr $sumq/($j+0.0)]
		if {[expr $sumF2/($j+0.0)]>[expr $F*$F]} {set err [expr sqrt(($sumF2/($j+0.0)-$F*$F)/($j-1.0))]} else {set err 0.000001}
		#if { $F<0 } { set F 0 }
		puts $fidout "$F $q $err \\"
	    }
	}    
    close $fidin
    close $fidout
}

proc correctq {fin sampl alpha} {
# does correction of q
    set fidin [open $fin r]
    set ftmp [open tmp w]
    set l ""  
  while { [gets $fidin l] > 0 } {
  puts $ftmp "$l"
  while { [gets $fidin l] > 0 } {
    if { [lindex $l 0]=="parameter" && [lindex $l 1]==$sampl } {
      puts $ftmp "$l"
      gets $fidin l
      puts $ftmp "$l"
      gets $fidin l
      puts $ftmp "$l"
      while { [gets $fidin l] > 0 } {
	 set F [lindex $l 0]
	 set tmpq [lindex $l 1]
	 set err [lindex $l 2]
	 set q [expr $tmpq*(1.0+$alpha)]
	 puts $ftmp "$F $q $err \\"
      }
    }
    puts $ftmp "$l"
  }
  puts $ftmp "$l"
  }
  close $ftmp
  close $fidin
  file rename -force tmp $fin   
}
