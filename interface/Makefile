###########################################################################
###  Femlisp/interface Makefile
###########################################################################
#
# Copyright (C) 2003 Nicolas Neuss, University of Heidelberg.
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
#
#  1. Redistributions of source code must retain the above copyright
# notice, this list of conditions and the following disclaimer.
#
#  2. Redistributions in binary form must reproduce the above copyright
# notice, this list of conditions and the following disclaimer in the
# documentation and/or other materials provided with the distribution.
#
# THIS SOFTWARE IS PROVIDED ``AS IS'' AND ANY EXPRESS OR IMPLIED
# WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
# MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN
# NO EVENT SHALL THE AUTHOR, THE UNIVERSITY OF HEIDELBERG OR OTHER
# CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
# PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
# LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
# NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
###########################################################################

SUPERLU_INCLUDE=/usr/include/superlu
SUPERLU_LIBRARY=/usr/lib
UMFPACK_INCLUDE=/usr/include/suitesparse
UMFPACK_LIBRARY=/usr/lib

SUPERLU_LIB=#-L $(SUPERLU_LIBRARY)
AMD_LIB=#-L $(UMFPACK_LIBRARY)

all: superlu umfpack

superlu: $(SUPERLU_INCLUDE)/slu_ddefs.h $(SUPERLU_LIBRARY)/libsuperlu.so
	gcc -c -fPIC -I$(SUPERLU_INCLUDE) superlu.c;
	ld -shared $(SUPERLU_LIB) -lsuperlu superlu.o -o superlu.so

umfpack: $(UMFPACK_INCLUDE)/umfpack.h
	gcc -I$(UMFPACK_INCLUDE) -Wall -fPIC -c umfpack.c;
	ld -shared $(UMFPACK_LIB) $(AMD_LIB) -lumfpack -lamd -lblas -llapack umfpack.o -o umfpack.so

mkl: mkl.c
	icpc -mkl -fPIC -shared -o mkl.so mkl.c

clean:
	rm -f superlu superlu.o superlu.so superlu.lgs superlu.solution
	rm -f umfpack umfpack.o umfpack.so
	rm -f a.out
