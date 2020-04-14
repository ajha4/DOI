SDK_DIR=$(HOME)/NVIDIA_GPU_Computing_SDK
INCLUDE=-I$(SDK_DIR)/C/common/inc/
LIB=-L$(SDK_DIR)/C/lib -L$(SDK_DIR)/C/common/lib/linux -l cutil -lcuda 
BIN_DIR=/usr/local/cuda/bin
NVCC=$(BIN_DIR)/nvcc
PTXAS=$(BIN_DIR)/ptxas
#FLAGS=-arch sm_20 -m 64 -g -G -gencode arch=compute_20,code=sm_20 --disable-warnings  
#FLAGS=-arch sm_20 -m 64 -g -use_fast_math -gencode arch=compute_20,code=sm_20 -Xcompiler -fopenmp --disable-warnings
FLAGS=-arch sm_20 -m 64 -use_fast_math -O3 -gencode arch=compute_20,code=sm_20 -Xcompiler -fopenmp --disable-warnings
 

rte=rte_gpu
inv_rte=rte_inv

INPUTS_RTE=rte.cu kernels.cu neumann.cu iniparser.cu dist.cu dictionary.cu inits.cu sph_harm.cu
INPUTS_INV=inverse_rte.cu kernels.cu neumann.cu iniparser.cu dist.cu dictionary.cu inits.cu sph_harm.cu

all: $(rte) $(inv_rte) 

$(rte): $(INPUTS_RTE)
	$(NVCC) $(FLAGS) $(INCLUDE) $(LIB) $(INPUTS_RTE) -o $(rte)

$(inv_rte): $(INPUTS_INV)
	$(NVCC) $(FLAGS) $(INCLUDE) $(LIB) $(INPUTS_INV) -o $(inv_rte)

clean: 
	rm -f *.o $(rte) $(inv_rte)
