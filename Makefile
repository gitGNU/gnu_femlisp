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
	echo "Options: doc, infix, superlu, umfpack, triangle";\
	echo "         femlisp, femlisp-core, slime, all, clean."

all: doc infix superlu umfpack triangle femlisp femlisp-core

doc:
	cd doc; make all

infix:
	cd external;\
	wget http://www.cs.cmu.edu/afs/cs/project/ai-repository/ai/lang/lisp/code/syntax/infix/infix.cl

superlu:
	cd interface; make superlu

umfpack:
	cd interface; make umfpack

triangle:
	echo "Installing Triangle in femlisp/external.  Note that Triangle	\
comes with a separate license which you should read (and accept) before		\
using it."
	cd external; mkdir triangle; cd triangle;\
	wget http://cm.bell-labs.com/netlib/voronoi/triangle.zip;\
	unzip triangle.zip; rm triangle.zip; make

femlisp:
	sed "/^FEMLISP_DIR=.*/c\FEMLISP_DIR=`pwd`" bin/femlisp >bin/femlisp2;\
	mv -f bin/femlisp2 bin/femlisp; chmod +x bin/femlisp;

femlisp-core:
	cd bin; rm -f femlisp.core;\
	sh ./femlisp -eval "(progn (ext:save-lisp \"femlisp.core\" :print-herald nil) (quit))"

slime:
	cd elisp; wget -O - http://common-lisp.net/project/slime/slime-1.0.tar.gz| tar xzvf -

clean:
	cd doc; make clean;
	cd interface; make clean;
	cd src; make clean;
