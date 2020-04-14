CC = gcc
LINKER = gcc
LIBS = -lpthread -lm 

CFLAGS = -O2 -fopenmp -m64 -Wall
#CFLAGS = -g 

RM            = /bin/rm

.c.o:
	$(CC) $(CFLAGS) -c $<

.C.o:
	$(C++) $(CFLAGS) -c $<

EXECS         = rte_cpu 

default: $(EXECS)


#RTE_uniform : RTE_uniform.o Inits.o SphHarm.o Dist.o Kernels_selfsub.o iniparser.o dictionary.o
#	$(LINKER) $(CFLAGS) -o $@ $(LIBS) $^

rte_cpu : RTE_uniform.o Inits.o SphHarm.o Dist.o Kernels_selfsub.o iniparser.o dictionary.o
	$(LINKER) $(CFLAGS) -o $@ $(LIBS) $^

#rte_sph : RTE_uniform.o Inits.o SphHarm.o Dist.o Kernels_sph.o iniparser.o dictionary.o
#	$(LINKER) $(CFLAGS) -o $@ $(LIBS) $^

clean:
	@-$(RM) -f *.o *~ PI* $(EXECS)
