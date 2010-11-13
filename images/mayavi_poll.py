#########################################################################
# mayavi_poll.py - Polls output.vtk
#########################################################################
#
# Copyright (c) 2006-2007, Enthought Inc.
# License: BSD Style.
#
# 2010 Femlisp adaption by Vitaly Polisky (KIT)
#
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
