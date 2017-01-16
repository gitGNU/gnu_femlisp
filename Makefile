###########################################################################
###  Femlisp Makefile
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

-include femlisp-config.mk

.PHONY: all clean cleanall configure documentation superlu umfpack femlisp mpi-worker

help:
	echo "Options: all, configure, clean, cleanall, documentation, femlisp,";\
	echo "         help, asdf, superlu, umfpack."

all: configure superlu umfpack femlisp documentation

asdf:
	cd external; make asdf

configure:
	./configure

quickload-libraries:
	$(FEMLISP_CL) --eval "(progn (ql:quickload (list :cl-ppcre :cl-gd :fiveam :bordeaux-threads :lparallel :cffi)) (quit))"

documentation:
	cd doc; $(MAKE) all

superlu:
	cd interface; $(MAKE) superlu

umfpack:
	cd interface; $(MAKE) umfpack

femlisp:
	$(FEMLISP_CL) --eval "(asdf:oos 'asdf:load-op :femlisp-save-core)"

mpi-worker:
	cd ./bin; rm -f mpi-worker; $(FEMLISP_CL) --eval "(asdf:oos 'asdf:load-op :femlisp-mpi-worker)"

clean:
	rm -f *.x86f *.fasl *.ufasl *.fas? *.fas *.o *.amd64f *.lx32fsl;
	cd bin; rm -f *-core mpi-worker mpi-worker-connection-data;
	cd doc; $(MAKE) clean;
	cd src; $(MAKE) clean;
	cd external; $(MAKE) clean;

cleanall: clean
	rm femlisp-config.mk
	cd external; $(MAKE) clean;
	cd interface; $(MAKE) clean;
