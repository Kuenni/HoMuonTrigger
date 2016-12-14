from plotting.PlotStyle import setPlotStyle, drawLabelCmsPrivateData,\
	drawLabelCmsPrivateSimulation,drawWaterMark
from plotting.OutputModule import CommandLineHandler, CliColors
from plotting.RootFileHandler import RootFileHandler

import os
import inspect
import subprocess
import shlex

class Plot:
	def __init__(self,filename = None,data = False, debug = False):
		setPlotStyle()
		self.commandLine = CommandLineHandler('[' + self.__class__.__name__ + ']')
		self.key = 'L1MuonPresent' if data else 'L1MuonTruth'
		self.data = data
		self.DEBUG = debug
		self.filename = filename
		if self.DEBUG:
			self.debug("Creating plot module %s" % self.__class__.__name__)
		if filename != None:
			self.fileHandler = self.createFileHandler(filename)
		pass
	
	def setModuleName(self,moduleName):
		self.fileHandler.setModuleName(moduleName)
	
	def createFileHandler(self,filename):
		fh = RootFileHandler(filename,debug=self.DEBUG)
		fh.printStatus()
		return fh
	
	def getGitCommitHash(self):
		cmd = shlex.split("git show --abbrev-commit")
		process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		gitOutput,err = process.communicate()
		gch = 'None'
		firstLine = gitOutput.split('\n')[0]
		if firstLine.find('commit') != -1:
			gch = firstLine.split(' ')[-1]
		return gch
	
	def createPlotSubdir(self,subdirname):
		if( not os.path.exists('plots')):
			if self.DEBUG:
				self.debug('Creating dir plots')
			os.mkdir('plots')
		gitCommitHash = self.getGitCommitHash()
		if( not os.path.exists('plots/' + gitCommitHash)):
			if self.DEBUG: self.debug('Creating dir plots/' + gitCommitHash)
			os.mkdir('plots/' + gitCommitHash)
		if( not os.path.exists('plots/' + gitCommitHash + '/' + subdirname)):
			if self.DEBUG: self.debug('Creating dir plots/' + gitCommitHash + '/' + subdirname)
			os.mkdir('plots/' + gitCommitHash + '/' + subdirname)
		self.plotSubdir = 'plots/' + gitCommitHash + '/' + subdirname
		return
	
	#Save a canvas as gif file with the source data file name attached
	def storeCanvas(self,canvas,plotname,drawMark = True, drawLabel = True, markPosition = None, labelPosition = None, marginRight=None, marginLeft=0.1):
		# Do the following to overcome the late binding in python which evaluates
		# the default objects only once and persists the information
		if markPosition == None:
			markPosition = {'x1ndc' : marginLeft - 0.007, 'y1ndc' : 0.95, 'x2ndc' : marginLeft + 0.226, 'y2ndc' : 0.99}
			
		#is this way also needed for "elementary data types"?
		if marginRight == None:
			marginRight = 0.1
			
		if labelPosition == None:
			labelPosition = {'x1ndc' : .7 - marginRight, 'y1ndc' : 0.95, 'x2ndc' : 1 - marginRight, 'y2ndc' : 0.98}
		
		if(plotname.find('/') != -1):
			if( not os.path.exists(self.plotSubdir + '/' + plotname[0:plotname.rfind('/')])):
				os.makedirs(self.plotSubdir + '/' + plotname[0:plotname.rfind('/')])
				
		canvas.cd().SetTopMargin(.05)
		canvas.cd().SetRightMargin(marginRight)
		canvas.cd().SetLeftMargin(marginLeft)

		if drawMark:
			mark = drawWaterMark(markPosition = markPosition)
		if drawLabel:		
			label = self.drawLabel(labelPosition = labelPosition)
		
		canvas.Update()
		canvas.SaveAs('%s/%s_%s.gif'%(self.plotSubdir,plotname,self.fileHandler.filename))
		canvas.SaveAs('%s/%s_%s.png'%(self.plotSubdir,plotname,self.fileHandler.filename))
		canvas.SaveAs('%s/%s_%s.pdf'%(self.plotSubdir,plotname,self.fileHandler.filename))
				
		return
	
	def drawLabel(self,labelPosition = {'x1ndc' : 0.6, 'y1ndc' : 0.95, 'x2ndc' : 0.9, 'y2ndc' : 0.98}):
		if self.data:
			label = drawLabelCmsPrivateData(labelPosition['x1ndc'],labelPosition['y1ndc'],labelPosition['x2ndc'],labelPosition['y2ndc'])
		else:
			label = drawLabelCmsPrivateSimulation(labelPosition['x1ndc'],labelPosition['y1ndc'],labelPosition['x2ndc'],labelPosition['y2ndc'])
		return label
	
	def debug(self,string):
		self.commandLine.debug(string)
	
	def warning(self,string):
		self.commandLine.warning(string)
		
	def error(self,string):
		self.commandLine.error(string)
		
	def output(self,string):
		self.commandLine.output(CliColors.BOLD + '<' + inspect.stack()[1][3] + '>  ' + CliColors.ENDC + str(string))
		
	def printProgress(self,done,total,updateEvery = 1000):
		if done == total or not (done % updateEvery):
			self.commandLine.printProgress(done, total)