function [map] = create_phantom(nx, ny, nz)

map = ones(nx,ny, nz);

% map( (nz-1) * (nx*ny) + (floor(ny/2) -1) * nx + (floor(nx/2) -1) +1 ) = 3;
map(1:5,:,:) = 1;
map(6:nz-6,:,:) = 2;
map(nz-5:nz,:,:) = 3;

map = reshape(map,nx*ny*nz,1);
fid = fopen('tiss_map_nonuni.bin','w+');
fwrite(fid, map, 'uchar');
fclose(fid);
copyfile tiss_map_nonuni.bin /home/abhinav/RapiDOT/InverseRTE_cuda/multi_gpu/RTE_GPU_Jan20/examples
copyfile tiss_map_nonuni.bin /home/abhinav/RapiDOT/InverseRTE_cuda/multi_gpu/RTE_GPU_Jan20/src

