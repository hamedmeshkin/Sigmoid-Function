clear
load conditional_pairs.dat
data = conditional_pairs;
clear conditional_pairs;
[c1, c2]=size(data);
data1 = zeros(c1*2,1);
data1(1:c1,:) = data(:,1);
data1(c1+1:end,:) = data(:,2);
data2=0;
l=0;
for i=1:length(data1)
    if any(data1(i)==data2)
    else
        l=l+1; data2(l,1)=data1(i);
    end
end

arng=sort(data2)+1; %+1 is because of 1 and 0 defers in apps I want to use it as pairs atomes in namd TclForce
%save usedatoms.dat arng '-ascii'
