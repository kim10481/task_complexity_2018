function [out]=randnum()

max_num_mat=[4, 70];

out=[];
for i=1:1:length(max_num_mat)
    out=[out ceil(rand(1)*max_num_mat(i))];
end
if length(max_num_mat)==2
    disp(sprintf('- (session#%d,trial#%d) has been selected.',out(1),out(2)));
end
end