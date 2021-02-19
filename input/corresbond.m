clear 
load conditional_pairs.dat;
alpha = conditional_pairs+1;
load usedatoms.dat
used = usedatoms;
clear conditional_pairs usedatoms

for i=1:length(alpha)
  heavy(i,1) =  find(used==(alpha(i,1)));
  heavy(i,2) =  find(used==(alpha(i,2)));
end
 heavy=heavy-1;