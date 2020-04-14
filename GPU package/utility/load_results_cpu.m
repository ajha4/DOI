function [tr] = load_results_cpu(nme)

nme1 = [nme 'Trans.out'];
fid = fopen(nme1,'r');
n = fread(fid,[1 2],'int');
tr = fread(fid,prod(n),'double');
tr = reshape(tr,n);
fclose(fid);
