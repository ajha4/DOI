#LIB_ARCH         = x86_64
SDK_DIR=$(HOME)/NVIDIA_GPU_Computing_SDK
INCLUDE=-I$(SDK_DIR)/C/common/inc/
LIB=-L$(SDK_DIR)/C/lib -L$(SDK_DIR)/C/common/lib/linux -l cutil -lcuda 
BIN_DIR=/usr/local/cuda/bin
NVCC=$(BIN_DIR)/nvcc
PTXAS=$(BIN_DIR)/ptxas
#FLAGS=-arch sm_20 -g -G -gencode arch=compute_20,code=sm_20 --disable-warnings  
#FLAGS=-arch sm_20 -m 64 -g -use_fast_math -gencode arch=compute_20,code=sm_20 -Xcompiler -fopenmp --disable-warnings -keep --ptxas-options=-v
FLAGS=-arch sm_20 -m 64 -O3 -use_fast_math -gencode arch=compute_20,code=sm_20 -Xcompiler -fopenmp --disable-warnings
 

rte=rte_gpu
inv_rte=rte_inv
rte_test = rte_test
rte_fim = rte_fim
rte_snr = rte_snr

INPUTS_RTE=rte.cu kernels.cu neumann.cu iniparser.cu dist.cu dictionary.cu inits.cu sph_harm.cu
INPUTS_INV=error_map.cu kernels.cu neumann.cu iniparser.cu dist.cu dictionary.cu inits.cu sph_harm.cu
INPUTS_TEST=compute_grad_rte.cu kernels.cu neumann.cu iniparser.cu dist.cu dictionary.cu inits.cu sph_harm.cu
INPUTS_FIM=rte_fim.cu kernels.cu neumann.cu iniparser.cu dist.cu dictionary.cu inits.cu sph_harm.cu
INPUTS_SNR=rte_snr.cu kernels.cu neumann.cu iniparser.cu dist.cu dictionary.cu inits.cu sph_harm.cu


#all: $(inv_rte) $(rte) $(rte_fim) $(rte_snr)
all: $(rte_fim)

$(rte): $(INPUTS_RTE)
	$(NVCC) $(FLAGS) $(INCLUDE) $(LIB) $(INPUTS_RTE) -o $(rte)

$(inv_rte): $(INPUTS_INV)
	$(NVCC) $(FLAGS) $(INCLUDE) $(LIB) $(INPUTS_INV) -o $(inv_rte)

$(rte_fim): $(INPUTS_FIM)
	$(NVCC) $(FLAGS) $(INCLUDE) $(LIB) $(INPUTS_FIM) -o $(rte_fim)

$(rte_snr): $(INPUTS_SNR)
	$(NVCC) $(FLAGS) $(INCLUDE) $(LIB) $(INPUTS_SNR) -o $(rte_snr)

$(rte_test): $(INPUTS_TEST)
	$(NVCC) $(FLAGS) $(INCLUDE) $(LIB) $(INPUTS_TEST) -o $(rte_test)

clean: 
	rm -f *.o $(rte) $(inv_rte) $(rte_test) $(rte_fim) $(rte_snr)
