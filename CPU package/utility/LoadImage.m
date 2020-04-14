function [img] = LoadImage(fname);

% Function M-file to read in an output image from UA RTE solver
% ex) img = LoadImage("output file name")
%
% img is the output image in matrix form

fid = fopen(fname,'r');

n = fread(fid, [1,2],'int');    % n is the x and y dimensions of the image
img = fread(fid, prod(n),'double');
img = reshape(img,n);
fclose(fid);