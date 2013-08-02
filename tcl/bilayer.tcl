package require BLT
namespace import blt::*
set zoomMod "Control-"
set tcl_precision 5
wm geometry . =900x650
########################
#general setup         #
########################
set actv {}
set commonusefont {Times 10}
set method 2
proc samplist {num name} {
	global listw ctfwN ctfwX rhow actv commonusefont
	set i $listw.$num
	frame $i
	$listw window create end -window $i
#	radiobutton $i.1 -width 2 -text $num -borderwidth 0  -font $commonusefont -variable actv -value $num -command "plot ctf $num; plot rho $num; colorp"
	label $i.1 -width 4 -text $num -borderwidth 0  -font "$commonusefont bold"
	#-width value determines the width of texts of sample names
	label $i.2 -text $name -width 25 -bg white -font $commonusefont
	bind $i.2 <1> "plot frm $num"
	pack $i.1 -side left
	pack $i.2 -side left -fill y
	newsample $num $name
}
proc deleteall {} {
	global listw ctfwN ctfwX rhow rhowN rhowX Tnum
	deletesamples
	for {set i 0} {$i<$Tnum} {incr i} {
		if {[winfo exists $listw.$i]} {
			$listw delete $listw.$i; 
			catch {$ctfwN element delete $i}
			catch {$ctfwN element delete m$i; $rhowN element delete s$i}
			catch {$ctfwX element delete $i}
			catch {$ctfwX element delete m$i; $rhowX element delete s$i}
		}
	}
}	
proc export2g {filename} {
	global x smpfilename fract Tnum chifactorN DBgibbs xpin SNUM status max_peak
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
	for {set i 0} {$i<$Tnum} {incr i} {
		puts $fid "set x($i) $x($i); [expr $xpin($i)?"fix":"free"] $i"
	}
	close $fid
}
proc active {num} {
	global status listw
	if {[winfo exists $listw.$num]} {
		if {$status($num)==0} {plot frm $num}
	}
}

proc report {num} {
	set w .rpt$num
	catch {destroy $w}
	toplevel $w
	wm title $w "$num report"
#	scrollbar $w.sbar -orient horizontal -command "$w.t xview"
#	text $w.t -xscrollcommand "$w.sbar set" -height 20 -width 14 -cursor arrow
#	pack $w.sbar $w.t -side bottom -fill x
	text $w.t -height 12 -width 70 -cursor arrow -tabs {1c 3c} 
	pack $w.t -expand yes -fill both
	getreport $num $w.t
}
########################
#setup for .t BLT graph#
########################
set w .t
set ctfwN $w.cN
set ctfwX $w.cX
set rhow $w.r
set rhowN $w.rN
set rhowX $w.rX
tabset $w -side top -borderwidth 0 -selectpad 0 -tabborderwidth 1 -highlightthickness 0
$w configure -font "Times 7"
graph $ctfwN -width 1604 -height 1504 -plotrelief flat -plotpadx {0 0} -plotpady {0 0} -plotborderwidth 0
graph $ctfwX -width 1604 -height 1504 -plotrelief flat -plotpadx {0 0} -plotpady {0 0} -plotborderwidth 0
graph $rhow -width 1604 -height 1504 -plotrelief flat -plotpadx {0 0} -plotpady {0 0} -plotborderwidth 0
graph $rhowN -width 1604 -height 1504 -plotrelief flat -plotpadx {0 0} -plotpady {0 0} -plotborderwidth 0
graph $rhowX -width 1604 -height 1504 -plotrelief flat -plotpadx {0 0} -plotpady {0 0} -plotborderwidth 0
$ctfwN axis configure x -min "" -max "" -loose 0 -subdivisions 10 -stepsize 0.1
$ctfwN axis configure y -min "" -max "" -loose 0 -subdivisions 10 -stepsize 1
$ctfwX axis configure x -min "" -max "" -loose 0 -subdivisions 10 -stepsize 0.1
$ctfwX axis configure y -min "" -max "" -loose 0 -subdivisions 10 -stepsize 1
$rhow axis configure x -subdivisions 5 -stepsize 5
$rhow axis configure y -min "" -max "" -subdivisions 5 -stepsize 0.05
$rhowN axis configure x -subdivisions 5 -stepsize 5
$rhowN axis configure y -min "" -max "" -subdivisions 5 -stepsize 0.05
$rhowX axis configure x -subdivisions 5 -stepsize 5
$rhowX axis configure y -min "" -max "" -subdivisions 5 -stepsize 0.05

$ctfwN legend configure -hide yes
$ctfwX legend configure -hide yes
$rhow legend configure -hide yes
$rhowN legend configure -hide yes
$rhowX legend configure -hide yes
$w insert 0 "CTFn" -window $ctfwN -padx 0 -pady 0
$w insert 1 "CTFx" -window $ctfwX -padx 0 -pady 0
$w insert 2 "PROFILE" -window $rhow -padx 0 -pady 0
$w insert 3 "SLD" -window $rhowN -padx 0 -pady 0
$w insert 4 "EDP" -window $rhowX -padx 0 -pady 0

.t.r element create mCG -xdata xmr -ydata ymCG -pixels 0 -linewidth 1 -color black
.t.r element create mPh -xdata xmr -ydata ymPh -pixels 0 -linewidth 1 -color red
.t.r element create mCh -xdata xmr -ydata ymCh -pixels 0 -linewidth 1 -color green
.t.r element create mc2 -xdata xmr -ydata ymc2 -pixels 0 -linewidth 1 -color blue
.t.r element create mc1 -xdata xmr -ydata ymc1 -pixels 0 -linewidth 1 -color cyan
.t.r element create mc3 -xdata xmr -ydata ymc3 -pixels 0 -linewidth 1 -color magenta
.t.r element create mw -xdata xmr -ydata ymw -pixels 0 -linewidth 1 -color yellow
.t.r element create mpep -xdata xmr -ydata ympep -pixels 0 -linewidth 1 -color #888800

bind $ctfwN	<Motion>	{set ycord [%W axis invtransform y %y]; set xcord [%W axis invtransform x %x]}
bind $ctfwX	<Motion>	{set ycord [%W axis invtransform y %y]; set xcord [%W axis invtransform x %x]}
bind $rhow	<Motion>	{set ycord [%W axis invtransform y %y]; set xcord [%W axis invtransform x %x]}
bind $rhowN	<Motion>	{set ycord [%W axis invtransform y %y]; set xcord [%W axis invtransform x %x]}
bind $rhowX	<Motion>	{set ycord [%W axis invtransform y %y]; set xcord [%W axis invtransform x %x]}



########################
#setup for .f          #
#menu, text widget     #
#flip sign, xy coords  #
########################
frame .f
set w .f.m
frame $w
set m $w.file.m
menubutton $w.file -text File -menu $m
menu $m -tearoff 0
$m add command -label "Open" -command {
   set types {
	{{smp Files}       {.smp}        }
	{{All Files}        *             }
    }
    set smpfilename [tk_getOpenFile -initialdir "./" -filetypes $types]
    cd [file dirname $smpfilename]
    update
    if {$smpfilename != ""} {
		deleteall
		source $smpfilename
    }
}
$m add command -label "reload" -command {deleteall; source $smpfilename}
$m add command -label "print CTFn" -command "Blt_PostScriptDialog $ctfwN"
$m add command -label "print CTFx" -command "Blt_PostScriptDialog $ctfwX"
$m add command -label "print PROFILE" -command "Blt_PostScriptDialog $rhow"
$m add command -label "print SLD" -command "Blt_PostScriptDialog $rhowN"
$m add command -label "print EDP" -command "Blt_PostScriptDialog $rhowX"
set m $w.tls.m
menubutton $w.tls -text Tool -menu $m
menu $m -tearoff 0
$m add command -label "export" -command {export all [file rootname $smpfilename]}
$m add command -label "peptide" -command {pepOptWin}
$m add command -label "randomization" -command {randOptWin}

set m $w.opt.m
menubutton $w.opt -text Opt -menu $m
menu $m -tearoff 0
$m add command -label "scale bar" -command {cfgbars}
#$m add check -label "aline on model" -variable aline_on_model
#$m add check -label "use frm as is" -variable as_is
$m add check -label "noscale" -variable noscale
$m add check -label "use sign" -variable usesign
#$m add check -label "HG constraint" -variable cnst(0)
#$m add check -label "MT constraint" -variable cnst(1)
#$m add cascade -label "Method" -menu $w.opt.m.radio
$m add cascade -label "hide/show" -menu $w.opt.m.hide
#set m $w.opt.m.radio
#menu $m -tearoff 0
#$m add radio -label "strict" -variable method -value 1
#$m add radio -label "baysian" -variable method -value 2
set m $w.opt.m.hide
menu $m -tearoff 0
$m add check -label "model cft" -variable hide_mcft -command {hide_model $hide_mcft}
#$m add check -label "model rho" -variable hide_mdr -command {$rhow element configure mdr -hide [expr $hide_mdr?"yes":"no"]}
$m add check -label "gaussian 1" -variable hide_mg1 -command {$rhow element configure mg1 -hide [expr $hide_mg1?"yes":"no"]}
$m add check -label "gaussian 2" -variable hide_mg2 -command {$rhow element configure mg2 -hide [expr $hide_mg2?"yes":"no"]}
$m add check -label "errorf c" -variable hide_mfc -command {$rhow element configure mfc -hide [expr $hide_mfc?"yes":"no"]}
$m add check -label "water" -variable hide_mfw -command {$rhow element configure mfw -hide [expr $hide_mfw?"yes":"no"]}
$m add check -label "methyl" -variable hide_mgm -command {$rhow element configure mgm -hide [expr $hide_mgm?"yes":"no"]}
$m add check -label "peptide" -variable hide_mgp -command {$rhow element configure mgp -hide [expr $hide_mgp?"yes":"no"]}

proc hide_model {hide_mcft} {
	global listw ctfwN ctfwX rhow rhowN rhowX
	for {set i 0} {$i<$SNUM} {incr i} {
		if {[winfo exists $listw.$i]} {
			catch { $ctfwN element configure m$i -hide [expr $hide_mcft?"yes":"no"] }
			catch { $ctfwX element configure m$i -hide [expr $hide_mcft?"yes":"no"] }
		}
	}
}


pack $w.file $w.tls $w.opt -side left

set w .f.t
frame $w
scrollbar $w.sbar -orient vertical -command "$w.t yview"
#-width value determines the size of the window that displays sample names
text $w.t -yscrollcommand "$w.sbar set" -height 20 -width 30 -cursor arrow -font $commonusefont
pack $w.sbar $w.t -side left -fill y

set listw $w.t

#disabled
#set w .f.p
#frame $w -borderwidth 3 -relief flat
#for {set i 1} {$i<16} {incr i} {
#	label $w.$i -text $i -bg white -width 2 -font $commonusefont
#	bind $w.$i <1> "flip $i"
#}
#grid $w.1 $w.2 $w.3 $w.4 $w.5 -padx 1 -pady 1
#grid $w.6 $w.7 $w.8 $w.9 $w.10 -pady 1
#grid $w.11 $w.12 $w.13 $w.14 $w.15 -pady 1

set w .f.n
frame $w
label .f.n.x -textvariable xcord -width 10 -anchor w -borderwidth 0 -padx 0 -pady 0 -font $commonusefont
label .f.n.y -textvariable ycord -width 10 -anchor w -borderwidth 0 -padx 0 -pady 0 -font $commonusefont
pack .f.n.x -side right
pack .f.n.y -side left

foreach i { 1 2 3 4 5} {
	entry .f.msg$i -textvariable msg$i -borderwidth 1 -width 16 -font $commonusefont
}
bind .f.msg1 <Return> {eval $msg1}
bind .f.msg2 <Return> {eval $msg2}
bind .f.msg3 <Return> {eval $msg3}
bind .f.msg4 <Return> {eval $msg4}
bind .f.msg5 <Return> {eval $msg5}

frame .f.par
#pack .f.msg5 .f.msg4 .f.msg3 .f.msg2 .f.msg1 .f.par .f.n .f.p -side bottom -fill x
pack .f.msg5 .f.msg4 .f.msg3 .f.msg2 .f.msg1 .f.par .f.n -side bottom -fill x
pack .f.m .f.t

foreach graph {.t.cN .t.cX .t.r .t.rN .t.rX} {
	Blt_ZoomStack $graph
	Blt_Crosshairs $graph
	Blt_PrintKey $graph
#shift 3
}

########################
#setup for .p          #
#entris scrollbar      #
########################
frame .p
#Notation here is that cfg(00), cfg(10), and cfg(20) are for XCG scroll bar setting.
#cfg(00) is the low limit, cfg(10) is the high limit, and cfg(20) is the increment.
set cfg(00)  0    ; set cfg(10)   30; set cfg(20)  0.0001;
#set cfg(01)  0    ; set cfg(11)   1 ; set cfg(21)  0.0001;
set cfg(02)  0.1  ; set cfg(12)   10; set cfg(22)  0.0001;
set cfg(03)  0    ; set cfg(13)   30; set cfg(23)  0.0001;
#set cfg(04)  0    ; set cfg(14)   1 ; set cfg(24)  0.0001;
set cfg(05)  0.1  ; set cfg(15)   10; set cfg(25)  0.0001;
set cfg(06)  0    ; set cfg(16)   30; set cfg(26)  0.0001;
#set cfg(07)  0    ; set cfg(17)   1 ; set cfg(27)  0.0001;
set cfg(08)  0.1  ; set cfg(18)   10; set cfg(28)  0.0001;
set cfg(09)  0    ; set cfg(19)   30; set cfg(29)  0.0001;
#set cfg(010) 0    ; set cfg(110)  1 ; set cfg(210) 0.0001;
set cfg(011) 0.1  ; set cfg(111)  10; set cfg(211) 0.0001;
set cfg(012) 0    ; set cfg(112)  30; set cfg(212) 0.0001;
#set cfg(013) 0    ; set cfg(113)  1 ; set cfg(213) 0.0001;
set cfg(014) 0.1  ; set cfg(114)  10; set cfg(214) 0.0001;
set cfg(015) 0    ; set cfg(115)  30; set cfg(215) 0.0001;
#set cfg(016) 0    ; set cfg(116)  1; set cfg(216) 0.0001;
set cfg(017) 0    ; set cfg(117)  10; set cfg(217) 0.0001;
set cfg(018) 0    ; set cfg(118)  200; set cfg(218) 0.1;
set cfg(019) 0    ; set cfg(119)  20; set cfg(219) 0.0001;
set cfg(020) 0    ; set cfg(120)  200; set cfg(220) 0.1;
set cfg(021) 0    ; set cfg(121)  20; set cfg(221) 0.0001;
set cfg(022) 0    ; set cfg(122)  200; set cfg(222) 0.1;
set cfg(023) 0    ; set cfg(123)  20; set cfg(223) 0.0001;

set pname(0) XCG
set pname(1) CCG
set pname(2) SCG
set pname(3) XPh
set pname(4) CPh
set pname(5) SPh
set pname(6) XCh
set pname(7) CCh
set pname(8) SCh
set pname(9) XC
set pname(10) CC
set pname(11) SC
set pname(12) Xc1
set pname(13) Cc1
set pname(14) Sc1
set pname(15) Xc3
set pname(16) Cc3
set pname(17) Sc3

set pname(18) r
set pname(19) r12
set pname(20) RCG
set pname(21) RPh
set pname(22) Rm
set pname(23) sigR

set pname(24) A
set pname(25) nC2
set pname(26) nC1
set pname(27) nC3
set pname(28) VL
set pname(29) VHL
set pname(30) Vpep
set pname(31) VH
set pname(32) VCG
set pname(33) VPh
set pname(34) VCh
set pname(35) Vc2
set pname(36) Vc1
set pname(37) Vc3
set pname(38) D
set pname(39) DPP
set pname(40) DB
set pname(41) DC

set pname(42) DH1
set pname(43) t
set pname(44) DH1
set pname(45) dXH
set pname(46) t
set pname(47) dXH
set pname(48) t
set pname(49) r
set pname(50) t
set pname(51) r12
set pname(52) t
set pname(53) RCG
set pname(54) t
set pname(55) RPh
set pname(56) t
set pname(57) SC
set pname(58) t
set pname(59) DC
set pname(60) dXH2
set pname(61) t
set pname(62) dXH2

set pname(63) X2
set pname(64) RSD

set pname(65) negP
set pname(66) StEP

set pname(67) dXH3
set pname(68) t
set pname(69) dXH3

proc sclconfig {w ext} {
	global pname cfg
	$w.l configure -text $pname($ext)
	$w.s configure -from $cfg(0$ext) -to $cfg(1$ext) -variable x($ext) -resolution $cfg(2$ext)
}

proc pin {ext w} {
	global xpin
	set xpin($ext) [expr ($xpin($ext)+1)%2]
	$w configure -fg [expr $xpin($ext)?"red":"cyan"]
}
proc fix {ext} {
	global xpin
	if ($xpin($ext)==0) {.p.s.1.l$ext configure -fg "red"}
	set xpin($ext) 1
}
proc free {ext} {
	global xpin
	set xpin($ext) 0
	.p.s.1.l$ext configure -fg "cyan"
}

proc parscale {w ext} {
	global x pname
	label $w.l$ext -text $pname($ext) -font {Times 10 bold} -borderwidth 0 -fg cyan
	bind $w.l$ext <3> "pin $ext $w.l$ext"
	bind $w.l$ext <1> "sclconfig .p.s.0 $ext"
	entry $w.e$ext -textvariable x($ext) -font {Times 10} -width 7 -borderwidth 0 -background #dfdfdf
}
proc parentry {name ext label parname} {
	global cfg
	label $name.l$ext -text $label
	foreach i {0 1 2 3} {
		entry $name.e$i$ext -width 6 -textvariable [lindex $parname $i]
	}
	grid $name.l$ext $name.e0$ext $name.e1$ext $name.e2$ext
}
proc valueE {w ext} {
	global pname x
	label $w.l$ext -text $pname($ext) -borderwidth 0 -font {Times 10 bold} -fg blue
	entry $w.e$ext -width 7 -textvariable x($ext) -borderwidth 0 -background #dfdfdf -font {Times 10}
}
proc valueF {w ext} {
	global fract
	label $w.lfract$ext -text fract$ext -borderwidth 0 -font {Times 10 bold} -fg blue
	entry $w.efract$ext -width 7 -textvariable fract($ext) -borderwidth 0 -background #dfdfdf -font {Times 10}
}
proc valueE2 {w ext var} {
	global pname $var
	label $w.l$ext -text $pname($ext) -borderwidth 0 -font {Times 10 bold} -fg blue
	entry $w.e$ext -width 6 -textvariable $var -borderwidth 0 -background #dfdfdf -font {Times 8 bold}
}
proc getbars2 {} {
	global cfg errf0 range av chisq x info pname
	bind $w.2.l9 <3> {}
}
proc cfgbars {} {
	global pname
	toplevel .cfgbars
	wm title .cfgbars "Configure parameter bars"
	foreach i {parameters minimum maximum resolution} {label .cfgbars.$i -text $i}
	grid .cfgbars.parameters .cfgbars.minimum .cfgbars.maximum .cfgbars.resolution
	foreach i {0 2 3 5 6 8 9 11 12 14 15 17 18 19 20 21 22 23} {parentry .cfgbars $i $pname($i) "cfg(0$i) cfg(1$i) cfg(2$i)"}
}


set w .p.i
frame $w
pack $w -fill x
label $w.am -text amoeba -pady 0 -relief raised -width 7
label $w.ft -text Fourier -pady 0 -relief raised -width 7
set repeat 10
proc amb1 {} {
	global repeat NMAX
	.p.i.am configure -relief sunken -fg red
	update
	set NMAX 10000
	for {set i 0} {$i<$repeat} {incr i} {amoeba}
	calcmodel
	display
        calculateEDP
	.p.i.am configure -relief raised -fg black
}
proc amb3 {} {
	calcmodel; display
}
proc scalesim {samplist} {
	scale2sim $samplist; display
}
bind $w.am <1> { amb1 }
bind $w.am <3> {
calcmodel
display
calculateEDP
}
proc ft {} {
#	catch {.t.r element delete mg1}
#	catch {.t.r element delete mg2}
#	catch {.t.r element delete mfc}
#	catch {.t.r element delete mfw}
#	catch {.t.r element delete mgm}
#	catch {.t.r element delete mgp}
	calcFourier
}
bind $w.ft <1> { ft }
pack $w.ft $w.am -side right

entry $w.msga -textvariable msga -borderwidth 1 -width 22 -font "Times 10"
entry $w.msgb -textvariable msgb -borderwidth 1 -width 22 -font "Times 10"
entry $w.msgc -textvariable msgc -borderwidth 1 -width 22 -font "Times 10"
#msgi added by KA
entry $w.msgi -textvariable msgi -borderwidth 1 -width 22 -font "Times 10"
entry $w.msgd -textvariable msgd -borderwidth 1 -width 22 -font "Times 10"
entry $w.msg6 -textvariable msg6 -borderwidth 1 -font "Times 10"
pack $w.msga $w.msgb $w.msgc $w.msgi $w.msgd -side right
pack $w.msg6 -side left -fill x -expand yes

set w .p.ii
frame $w
pack $w -fill x
label $w.am -text "Neutrons" -pady 0 -relief flat -width 15
pack $w.am -side right
entry $w.msge -textvariable msge -borderwidth 1 -width 22 -font "Times 10"
entry $w.chin -textvariable chifactorN -borderwidth 1 -width 22 -font "Times 10"
entry $w.msgf -textvariable msgf -borderwidth 1 -width 45 -font "Times 10"
entry $w.msg7 -textvariable msg7 -borderwidth 1 -font "Times 10"
pack $w.msge $w.chin $w.msgf -side right
pack $w.msg7 -side left -fill x -expand yes

set w .p.iii
frame $w
pack $w -fill x
label $w.am -text "X-rays" -pady 0 -relief flat -width 15
pack $w.am -side right
entry $w.msgg -textvariable msgg -borderwidth 1 -width 22 -font "Times 10"
entry $w.chix -textvariable chifactorX -borderwidth 1 -width 22 -font "Times 10"
entry $w.msgh -textvariable msgh -borderwidth 1 -width 45 -font "Times 10"
entry $w.msg8 -textvariable msg8 -borderwidth 1 -font "Times 10"
pack $w.msgg $w.chix $w.msgh -side right
pack $w.msg8 -side left -fill x -expand yes

bind .p.i.msg6 <Return> {eval $msg6}
bind .p.ii.msg7 <Return> {eval $msg7}
bind .p.iii.msg8 <Return> {eval $msg8}

set w .p.s
frame $w
pack $w -fill x
frame $w.0;	pack $w.0 -side bottom -fill x
frame $w.1;	pack $w.1 -side left
frame $w.2;	pack $w.2 -side left
label $w.0.l -text $pname(0) -font {Times 8 bold} -foreground #007700 -borderwidth 0
scale $w.0.s -showvalue 0 -borderwidth 1 -from $cfg(00) -to $cfg(10) -variable x(0) -orient horizontal -resolution $cfg(20) -width 12
bind $w.0.s  <ButtonRelease-1> {mdplot}
bind $w.0.s  <B1-Motion> {mdplot}
pack $w.0.l -side left
pack $w.0.s -side right -fill x -expand yes
foreach i {0 2 3 5 6 8 9 11 12 14 15 17 18 19 20 21 22 23} {
	parscale $w.1 $i
}
foreach i {1 4 7 10 13 16} {
	valueE $w.1 $i
}
for {set i 24} {$i<$Tnum} {incr i} {
	valueE $w.1 $i
}
for {set i 0} {$i<3} {incr i} {
	valueF $w.1 $i
}

label $w.1.div
label $w.1.div2
#     0                  1                  2                  3                  4                  5                  6                  7                  8                  9                  10                 11
grid  $w.1.l0  $w.1.e0   $w.1.l3  $w.1.e3   $w.1.l6  $w.1.e6   $w.1.l9  $w.1.e9   $w.1.l57 $w.1.e57  $w.1.l12 $w.1.e12  $w.1.l15 $w.1.e15  $w.1.l49 $w.1.e49  $w.1.l51 $w.1.e51  $w.1.l53 $w.1.e53  $w.1.l55 $w.1.e55  $w.1.l22 $w.1.e22
grid  $w.1.l1  $w.1.e1   $w.1.l4  $w.1.e4   $w.1.l7  $w.1.e7   $w.1.l10 $w.1.e10  $w.1.l56 $w.1.e56  $w.1.l13 $w.1.e13  $w.1.l16 $w.1.e16  $w.1.l48 $w.1.e48  $w.1.l50 $w.1.e50  $w.1.l52 $w.1.e52  $w.1.l54 $w.1.e54  $w.1.l23 $w.1.e23
grid  $w.1.l2  $w.1.e2   $w.1.l5  $w.1.e5   $w.1.l8  $w.1.e8   $w.1.l11 $w.1.e11  $w.1.div $w.1.div  $w.1.l14 $w.1.e14  $w.1.l17 $w.1.e17  $w.1.l18 $w.1.e18  $w.1.l19 $w.1.e19  $w.1.l20 $w.1.e20  $w.1.l21 $w.1.e21

grid $w.1.div2
grid  $w.1.l38 $w.1.e38  $w.1.l28 $w.1.e28  $w.1.l59 $w.1.e59  $w.1.l44 $w.1.e44  $w.1.l47 $w.1.e47  $w.1.l62 $w.1.e62  $w.1.l69 $w.1.e69  $w.1.l32 $w.1.e32  $w.1.l35 $w.1.e35  $w.1.l24 $w.1.e24  $w.1.l65 $w.1.e65  $w.1.lfract2 $w.1.efract2
grid  $w.1.l25 $w.1.e25  $w.1.l29 $w.1.e29  $w.1.l58 $w.1.e58  $w.1.l43 $w.1.e43  $w.1.l46 $w.1.e46  $w.1.l61 $w.1.e61  $w.1.l68 $w.1.e68  $w.1.l33 $w.1.e33  $w.1.l36 $w.1.e36  $w.1.l39 $w.1.e39  $w.1.l66 $w.1.e66  $w.1.lfract0 $w.1.efract0
grid  $w.1.l26 $w.1.e26  $w.1.l30 $w.1.e30  $w.1.l41 $w.1.e41  $w.1.l42 $w.1.e42  $w.1.l45 $w.1.e45  $w.1.l60 $w.1.e60  $w.1.l67 $w.1.e67  $w.1.l34 $w.1.e34  $w.1.l37 $w.1.e37  $w.1.l40 $w.1.e40  $w.1.div $w.1.div  $w.1.lfract1 $w.1.efract1

foreach i {25 26 27 28 29 30 38 43 44 46 47 48 49 50 51 52 53 54 55 56 57 58 59 61 62 66 68 69} {
	$w.1.l$i configure -fg #008800
}

######################
#samll calculations  #
######################

proc correct {alpha} {
# does final correction
	global x
#shift the peak positions
	set x(0) [expr $x(0)/(1.0+$alpha)]
	set x(3) [expr $x(3)/(1.0+$alpha)]
	set x(6) [expr $x(6)/(1.0+$alpha)]
	set x(9) [expr $x(9)/(1.0+$alpha)]
	set x(12) [expr $x(12)/(1.0+$alpha)]
#scale the widths
	set x(2) [expr $x(2)-$x(0)*$alpha]
	set x(5) [expr $x(5)-$x(3)*$alpha]
	set x(8) [expr $x(8)-$x(6)*$alpha]
	set x(11) [expr $x(11)-$x(9)*$alpha]
	set x(14) [expr $x(14)-$x(12)*$alpha]
#calculate a new area
#and scale the amplitudes
	amb3
}

proc raiseq {iw i} {
	global hide_mcft
	set t [lindex [split $iw {}] 4]

	set w .t.c$t
	set ql [$w element show]
	set rq [lsearch $ql $i]
	set ql [lreplace $ql $rq $rq]
	set ql [linsert $ql end $i]
	catch {$w element show $ql}

	set j m$i
	set ql [$w element show]
	set rq [lsearch $ql $j]
	set ql [lreplace $ql $rq $rq]
	set ql [linsert $ql end $j]
	catch {$w element show $ql}
	hide_model $hide_mcft

	set w .t.r$t
	set j s$i
	set ql [$w element show]
	set rq [lsearch $ql $j]
	set ql [lreplace $ql $rq $rq]
	set ql [linsert $ql end $j]
	catch {$w element show $ql}
}

pack .p -side bottom -fill x
pack .f .t -side left -fill y -fill x

proc myfit {s x1 x2 iter} {
	for {set i 0} {$i<$iter} {incr i} {
		puts "iteration $i"
		amb3
		borrow1 $s $x1 $x2 0
		update
		amb1
	}
}

set showsign 0
set usesign 1
set program_mode 3
source frm2smp.tcl
source SIMtoEXP.tcl
set redraw 1

foreach ext {15 22 23} {
    pin $ext .p.s.1.l$ext
}

#cd /mnt/E/work/chess04b
source optionwin.tcl
