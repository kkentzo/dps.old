CC=gcc
CFLAGS=-Wall -g # -pg
LDFLAGS=-lgsl -lgslcblas -lhdf5 -lhdf5_hl `pkg-config --cflags --libs glib-2.0`

dps: dps.c params.c logger.c cell.c plasmid.c population.c pool.c lib/clib.c lib/sampling.c lib/nvar.c lib/qhdf.c lib/bitvector.c
