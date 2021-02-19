#mol pdbload 2F4K

mol load psf ../../common/TrpCage_mut_wb_ion.psf pdb ../../common/TrpCage_mut_wb_ion.pdb

set fil [open conditional_pairs.dat w]
set foo [open conditional_coordinate.dat w]


##### all the protein's atoms in atomselection
set sel [atomselect top protein]                     

##### all carbon alpha in the protein
set sel1 [atomselect top "protein not hydrogen"]         

##### distance between two atoms less than cutoff point
set sll [measure contacts 4.5 $sel1 $sel1]


set atomsA [lindex $sll 0]
set atomsB [lindex $sll 1]
#lindex [$sel1 get {x y z}] 1
set length [llength $atomsA]
for {set r 0} {$r <= $length } {incr r} {


set atomA [lindex $sll 0 $r]
set atomB [lindex $sll 1 $r] 

##### Add one because it starts at 0
set num1 [expr $atomA+0] 
set num2 [expr $atomB+0]


set s1 [lindex [$sel get resid] $num1]
set s2 [lindex [$sel get resid] $num2]



if {(abs($s1-$s2))>3} {

#####get distance
set sd1 [lindex [$sel get {x y z}] $num1]
set sd2 [lindex [$sel get {x y z}] $num2]

set dis [veclength [vecsub $sd1 $sd2]]

set s1 [expr $s1-1]
set s2 [expr $s2-1]

#set noatmA [expr $num1-1]
#set noatmB [expr $num2-1]
puts -nonewline  $foo " $dis \t"
puts   $fil "$num1	$num2	\t"}

}
puts -nonewline $foo  "\r"

close $fil
set no [$sel molid]
mol delete $no
# ######################################################################################
# # next lines find those pairs that have the condition of cutoff point for each frame #
# ######################################################################################
# mol load psf ../../common/BBA_wb_ion.psf dcd ../../box/BBA_wb_eq.dcd
# set nf [molinfo top get numframes]
# 
# for {set x 0} {$x < $nf} {incr x} {
# set frames [atomselect top protein frame $x] 
# #######################################
# for {set r 0} {$r <= $length } {incr r} {
# set atomA [lindex $sll 0 $r]
# set atomB [lindex $sll 1 $r] 
# set num1 [expr $atomA+0] 
# set num2 [expr $atomB+0]
# set s1 [lindex [$frames get resid] $num1]
# set s2 [lindex [$frames get resid] $num2]
# if {(abs($s1-$s2))>3} {
# 
# set dist1 [lindex [$frames get {x y z}] $num1]
# set dist2 [lindex [$frames get {x y z}] $num2]
# 
# set dist [veclength [vecsub $dist1 $dist2]]
# 
# puts -nonewline  $foo " $dist \t"
# 
# }
#  
# }
# 
# #################
# puts -nonewline $foo  "\r"
# 
# }

close $foo
exit
