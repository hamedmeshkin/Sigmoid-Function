
mol new ../common/pbt_wb.psf
mol addfile pbt_wb_eq.dcd waitfor all

set fram {400 796 1198 1596 1981 2400 2795 3198 3595 3994 7394 4799 5196 5589 5983 6388 6794 7198 7594 7994 8391 8797 9196 9596 10000 10400 10800 11200 11600 12000 12400 12790}
for {set fr 0} {$fr < [llength $fram]} {incr fr} {
   set num [lindex $fram $fr]
   molinfo top set frame $num
   set sel [atomselect top all frame $num]
  
   #set sed [measure center $sel]
   set x [expr [molinfo top get a frame $num]+0.0]
   set y [expr [molinfo top get b frame $num]+0.0]
   set z [expr [molinfo top get c frame $num]+0.0]
   
   set fname [format %03d $fr]
   set s_size [open initial/$fname.xsc w]
   
   puts $s_size "# NAMD extended system configuration restart file"
   puts $s_size "#\$LABELS step a_x a_y a_z b_x b_y b_z c_x c_y c_z o_x o_y o_z s_x s_y s_z s_u s_v s_w"
   puts $s_size "$num $x 0 0 0 $y 0 0 0 $z 0 0 0 0 0 0 0 0 0"
   $sel writepdb initial/$fname.pdb
   close $s_size
}

exit