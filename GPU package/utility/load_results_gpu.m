
function [tr,ref] = load_results_gpu(nme)

nme1 = [nme 'Trans.out'];
fid = fopen(nme1,'r');
n = fread(fid,[1 2],'int');
tr = fread(fid,prod(n),'double');
tr = reshape(tr,n);
fclose(fid);


nme2 = [nme 'Ref.out'];
fid = fopen(nme2,'r');
n = fread(fid,[1 2],'int');
ref = fread(fid,prod(n),'double');
ref = reshape(ref,n);
fclose(fid);


%nme3 = [nme 'Phase.out'];
%fid = fopen(nme3,'r');
%n = fread(fid,[1 2],'int');
%ph = fread(fid,prod(n),'double');
%ph = reshape(ph,n);
%fclose(fid);
