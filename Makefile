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
	echo "Options: doc, infix, slime, femlisp-core, triangle, clean."

doc:
	cd doc; make all

infix:
	cd external;\
	wget http://www.cs.cmu.edu/afs/cs/project/ai-repository/ai/lang/lisp/code/syntax/infix/infix.cl

slime:
	cd elisp; wget -O - http://common-lisp.net/project/slime/slime-1.0.tar.gz| tar xzvf -

femlisp-core:
	cd bin; rm femlisp.core;\
	sh ./femlisp -eval "(progn (ext:save-lisp \"femlisp.core\" :print-herald nil) (quit))"

triangle:
	echo "Installing Triangle in femlisp/external.  Note that Triangle	\
comes with a separate license which you should read (and accept) before		\
using it."
	cd external; mkdir triangle; cd triangle;\
	wget http://cm.bell-labs.com/netlib/voronoi/triangle.zip;\
	unzip triangle.zip; rm triangle.zip; make

superlu:
	cd interface;\
	gcc -c -fPIC -I/usr/include/superlu superlu.c;\
	ld -lsuperlu -shared superlu.o -o superlu.so

umfpack:
	cd interface;\
	gcc -I/usr/include/umfpack/ -Wall -fPIC -c umfpack.c;\
	ld -lumfpack -lamd -shared umfpack.o -o umfpack.so

clean:
	cd src; make clean;
	cd doc; make clean;

