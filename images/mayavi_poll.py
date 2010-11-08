###########################################################################
### mayavi_poll.py - Polls output.vtk
###########################################################################
###
### Copyright (C) 2010
### Vitaly Polisky, Nicolas Neuss, Karlsruhe Institute of Technology
###
### All rights reserved.
### 
### Redistribution and use in source and binary forms, with or without
### modification, are permitted provided that the following conditions
### are met:
### 
### 1. Redistributions of source code must retain the above copyright
### notice, this list of conditions and the following disclaimer.
### 
### 2. Redistributions in binary form must reproduce the above copyright
### notice, this list of conditions and the following disclaimer in the
### documentation and/or other materials provided with the distribution.
### 
### THIS SOFTWARE IS PROVIDED ``AS IS'' AND ANY EXPRESS OR IMPLIED
### WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
### MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
### IN NO EVENT SHALL THE AUTHOR, THE KARLSRUHE INSTITUTE OF TECHNOLOGY,
### OR OTHER CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
### SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
### LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
### DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
### THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
### (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
### OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
###
###########################################################################

import os
from os.path import join, abspath, dirname
from enthought.mayavi.scripts import mayavi2
from enthought.mayavi.sources.vtk_file_reader import VTKFileReader
#from enthought.mayavi.modules.outline import Outline
from enthought.mayavi.modules.api import Surface
from enthought.mayavi.modules.contour_grid_plane import ContourGridPlane
from enthought.pyface.timer.api import Timer

class Pollster(object):

	def __init__(self, fname, data):
		self.fname = fname
	        self.data = data
	        self.last_stat = os.stat(fname)
	
	def poll_file(self):
	        s = os.stat(self.fname)
	        if s[-2] == self.last_stat[-2]:
	            return
	        else:
	            self.last_stat = s
	            self.update_pipeline()
	
	def update_pipeline(self):
		print "file changed"
	        d = self.data
	        d.reader.modified()
	        d.update()
	        d.data_changed = True
	
		                       
def setup_data(fname):
	mayavi.new_scene()
	d = VTKFileReader()
	d.initialize(fname)
	mayavi.add_source(d)
	return d
	
def view_data():
	s = Surface()
	mayavi.add_module(s)
	s.module_manager.scalar_lut_manager.show_scalar_bar = True
	
@mayavi2.standalone
def main():

	fname = 'output.vtk'
	
	data = setup_data(fname)
	view_data()
	

	p = Pollster(fname, data)
	timer = Timer(1000, p.poll_file)
	mayavi2.savedtimerbug = timer
	
	    # To stop polling the file do:
	    #timer.Stop()
	
if __name__ == '__main__':
	main()
