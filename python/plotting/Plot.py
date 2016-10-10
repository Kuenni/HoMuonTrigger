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
	def storeCanvas(self,canvas,plotname,drawMark = True):
		if(plotname.find('/') != -1):
			if( not os.path.exists(self.plotSubdir + '/' + plotname[0:plotname.rfind('/')])):
				os.makedirs(self.plotSubdir + '/' + plotname[0:plotname.rfind('/')])
				
		canvas.cd()
		if drawMark:
			mark = drawWaterMark()
		canvas.SaveAs('%s/%s_%s.gif'%(self.plotSubdir,plotname,self.fileHandler.filename))
		canvas.SaveAs('%s/%s_%s.png'%(self.plotSubdir,plotname,self.fileHandler.filename))
		return
	
	def drawLabel(self,x1ndc = 0.6, y1ndc = 0.90, x2ndc = 0.9, y2ndc = 0.93):
		label = None
		if self.data:
			label = drawLabelCmsPrivateData(x1ndc,y1ndc,x2ndc,y2ndc)
		else:
			label = drawLabelCmsPrivateSimulation(x1ndc,y1ndc,x2ndc,y2ndc)
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