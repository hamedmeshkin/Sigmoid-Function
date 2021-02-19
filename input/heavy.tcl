#mol pdbload 5PTI

mol load psf ../common/pbt_wb.psf pdb ../common/pbt_wb.pdb

set fil [open heavy_atoms.dat w]


##### all the protein's atoms in atomselection
set sel [atomselect top "protein not hydrogen"]         
set xyz [$sel get index]
puts $fil $xyz
close $fil
exit
