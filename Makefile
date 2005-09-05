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

help:
	echo "Options: all, configure, clean, documentation, femlisp, help,";\
	echo "         triangle, superlu, umfpack."

all: configure superlu umfpack triangle femlisp documentation

configure:
	cd bin; ./femlisp-configure

documentation:
	cd doc; make all

download_superlu:
	echo "Downloading SuperLU in femlisp/external."
	cd external; \
	wget http://crd.lbl.gov/~xiaoye/SuperLU/superlu_3.0.tar.gz; \
	tar xzf superlu_3.0.tar.gz; rm -f superlu_3.0.tar.gz;

superlu:
	cd interface; make superlu

triangle:
	echo "Installing Triangle in femlisp/external.  Note that Triangle	\
comes with a separate license which you should read (and accept) before		\
using it."
	cd external; mkdir triangle; cd triangle;\
	wget http://cm.bell-labs.com/netlib/voronoi/triangle.zip;\
	unzip triangle.zip; rm triangle.zip; make

umfpack:
	cd interface; make umfpack

femlisp:
	cd bin; sh ./femlisp --save-core-and-die

slime:
	cd elisp; wget -O - http://common-lisp.net/project/slime/slime-1.0.tar.gz| tar xzvf -

clean:
	rm -f *.x86f *.fasl *.fas?
	cd bin; rm -f *.core
	cd doc; make clean;
	cd external; make clean;
	cd interface; make clean;
	cd private; make clean;
	cd src; make clean;
