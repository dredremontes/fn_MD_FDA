# copyright Bogdan Costescu 2010-2012

# colors atoms for each timestep based on data from external file
# run as:
#
# package require pf_loaduser
# pf_loaduser "peratom-sum.pfa"
#
# the file contains a line per timestep, each lines contains N floating point
# numbers, where N is the number of atoms, separated by spaces
#
# the color_min and color_max represent the range of colors to scale;
# if they are -1, it means to use the value found in the file, otherwise
# use the value given

package provide pf_loaduser 20120829

proc pf_min_ignore_zero { numbers } {
	set numberslen [llength $numbers]
	if {$numberslen == 0} {
		return 0.0
	}
	set idx 0
	set min 0.0
	while {$idx < $numberslen} {
		set x [lindex $numbers $idx]
		if {$x != 0.0} {
			set min $x
			break
		}
		incr idx
	}
	if {$min == 0.0} {
		return 0.0
	}
	while {$idx < $numberslen} {
		set x [lindex $numbers $idx]
		if {$x != 0.0 && $x < $min} {
			set min $x
		}
		incr idx
	}
	return $min
}

proc pf_loaduser { fname {residues_renumbered false} {ignore_zeros 0} {color_scale "BGR"} {color_min -1} {color_max -1} } {
	# check whether it's a file containing residue data
	set fnamelen [string length $fname]
	if {[string compare [string tolower [string range $fname [expr $fnamelen - 4] $fnamelen]] ".psr"] == 0} {
		set pf_residues 1
	} else {
		set pf_residues 0
	}
	if {$pf_residues == 1} then {
		puts "now loading per-residue data"
	} else {
		puts "now loading per-atom data"
	}

	# get rid of the first frame which is the PDB structure; this way the
	# framecount below is the nr. of frames in the trajectory file
	#animate delete beg 0 end 0 skip 0 0

	# make sure the representation is set to coloring by User and color is		# updated for every frame
	mol modcolor 0 top User
	mol colupdate 0 top 1
	color scale method $color_scale
	# find out the nr. of frames loaded and make a selection composed of all atoms
	set framecount [molinfo top get numframes]
	set sel [atomselect top all]

	# min is set to a big number, so actual numbers are smaller than it
	set color_scale_min 9.9e20
	# data should always be positive, so it's safe to set max to 0.0
	set color_scale_max 0.0

	set dat [open $fname r]
        gets $dat file_type

	for {set frame_no 0} {$frame_no<$framecount} {incr frame_no} {
		# color by data read from file
		# read data from file; one line per frame, atom values separated by space
		#set data($i) [gets $dat]
		if {[gets $dat one_frame] < 0} {
			puts "No user data for frame $frame_no"
			continue
		}
		if {$ignore_zeros != 0} {
			set frame_min [pf_min_ignore_zero $one_frame]
		} else {
			set frame_min [::tcl::mathfunc::min {*}$one_frame]
		}
		set color_scale_min [::tcl::mathfunc::min $color_scale_min $frame_min]
		set frame_max [::tcl::mathfunc::max {*}$one_frame]
		set color_scale_max [::tcl::mathfunc::max $color_scale_max $frame_max]
		#puts "per frame extremes: $frame_min $frame_max"
		if {$pf_residues == 1} then {
			for {set idx 0} {$idx < [llength $one_frame]} {incr idx} {
				if {$residues_renumbered} {
					set ressel [atomselect top "residue $idx"]
				} else {
					set ressel [atomselect top "resid $idx"]
				}
				$ressel frame $frame_no
				$ressel set user [lindex $one_frame $idx]
				$ressel delete
			}
		} else {
			$sel frame $frame_no
			$sel set user $one_frame
		}
		#$sel set user [gets $dat]
		# color by distance between atom and COM of selection
#		set com [measure center $sel]
#		set dlist ""
#		foreach c [$sel get {x y z}] {
#			lappend dlist [veclength [vecsub $c $com]]
#		}
#		$sel set user $dlist
		# color by time percentage
		# this works because the min/max points of the scale are set to 0 and 1
		# respectively by the line 'mol scaleminmax top 0 0.000000 1.000000'
#		$sel set user [expr 1.0*$frame_no/$framecount]
	}

	if { $color_min == -1} {
		set color_min $color_scale_min
	}
	if { $color_max == -1} {
		set color_max $color_scale_max
	}
	puts "setting $color_scale color scale to: $color_min $color_max"
	mol scaleminmax top 0 $color_min $color_max
	close $dat
}

