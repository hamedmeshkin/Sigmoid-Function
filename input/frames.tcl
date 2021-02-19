set fd [open sumN.dat r]
  set pair [read $fd]
  close $fd
set a [llength $pair]
mol load psf ../common/pbt_wb.psf dcd pbt_wb_eq.dcd
set nf [molinfo top get numframes]
puts "$a  $nf"
exit