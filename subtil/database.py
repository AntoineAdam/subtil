try:
	import mysql.connector
except ImportError as ex:
	print ex
	print "You need to install the mysql-connector-python to be able to use a mysql database."

import traceback
import time
import getpass
import sqlite3
import csv
import numpy as np
from os import path

from . import data
from .utils import is_boolean,is_float

_loaded=[]
_cells=np.array([])
_lineages=np.array([])

#if path.exists("dbconfig.py"):
#	execfile("dbconfig.py")
#elif path.exists(path.join(path.dirname(__file__), "dbconfig.py")):
#	execfile(path.join(path.dirname(__file__), "dbconfig.py"))
#else:
_dbtype="sqlite"
_dbfile="subtil.db"

_sqlitedatatypes={"BOOL":"nominal",
							"INTEGER":"numeric",
							"DOUBLE":"numeric",
							"TEXT":"nominal"}

def set_database(dbtype, dblocation=None, dbname=None, dbuser=None, dbpassword=None):

	"""Connect to a different database than the default one ("subtil.db" sqlite file).

	Parameters
	-------------
	dbtype: string
		The type of the sql database to connect to:
		Supported ones are "sqlite" or "mysql"
	
	"""

	global _dbtype
	if dbtype=="sqlite":
		if dblocation is None:
			print("You need to specify a file")
			return

		global _dbfile
		_dbtype = dbtype
		_dbfile = dblocation

	elif dbtype == "mysql":
		if dblocation is None:
			print("You need to specify a host location")
			return
		if dbname is None:
			print("You need to specify a name")
			return
		if dbuser is None:
			print("You need to specify a user")
			return
		if dbpassword is None:
			dbpassword = getpass.getpass("Password for {} on {}:".format(dbuser,dblocation))
	
		global _dbhost,_dbname,_dbuser,_dbpassword
		_dbtype = dbtype
		_dbhost = dblocation
		_dbname = dbname
		_dbuser = dbuser
		_dbpassword = dbpassword

	else:
		print("database type not supported: {}".format(type))


def _connect():

	if _dbtype=="mysql":
		#TODO: check if db exists. Create if not.
	
		return mysql.connector.connect(user=_dbuser, password=_dbpassword, host=_dbhost, database=_dbname)
	if _dbtype=="sqlite":
		#TODO create if doesn't exists

		connection = sqlite3.connect(_dbfile)
		connection.execute("pragma foreign_keys=ON")

		if connection.execute("select name from sqlite_master where type='table' and name = 'Experiment'").fetchone()==None:
			with open(path.join(path.dirname(__file__), "subtildb-structure-sqlite.sql")) as f:
				script = "".join(f.readlines())
			connection.executescript(script)
		

		return connection


def lookup():

	"""Look up for experiments in the database

	Prints and returns the experiments as an array.
	Each line is an experiments.
	Columns are: id, date of the experiment, timestep between 2 time-lapse picture, timestep between 2 fluorescence picture, number of cells in the experiment, comments

	Returns
	-------
	experiments : array
		The experiments previously uploaded to the database.
	"""
	
	cnx = _connect()
	cursor=cnx.cursor()
	
	cursor.execute("select e.id,date,timestep,fluotimestep,count(c.id),e.comments from experiment e left outer join lineage l on l.experimentid= e.id left outer join cell c on c.lineageid=l.id group by e.id")
	rows=cursor.fetchall()
	cursor.close()
	cnx.close()

	print("{:3}|{:11}|{:8}|{:12}|{:8}|{:20}".format("id","date","timestep","fluotimestep","numcells","comments"))
	for exp in rows:
		try:
			if _dbtype == 'sqlite':
				date = exp[1][:10]
			elif _dbtype == 'mysql':
				date = "{:4d}-{:2d}-{:2d}".format(exp[1].year,exp[1].month,exp[1].day)
			print("{:3}|{:11}|{:8}|{:12}|{:8}|{:20}".format(exp[0],date,exp[2],exp[3],exp[4],exp[5]))
		except Exception as ex:
			print traceback.format_exc()
	return rows

def create_experiment(date, comments, scale = 0.11, timestep=1, fluotimestep="NULL"):

	"""Creates a new experiment in the database.

	Parameters
	----------
	date: string
		The date of the experiment
		Format should be 'yyyy-mm-dd'.

	comments: string
		Comments on the experiments: Type of cell, medium, ...

	scale: double
		Spatial scale of the experiments. Default is 0.11 um/px.

	timestep: double
		Time step between 2 frames in minutes. Default is 1.

	fluotimestep: double
		Time step between 2 fluorescence measurment in minutes. Default is 4.

	Returns
	-------
	id: int
		The identifier of the created experiment

	Examples
	--------
	>>> create_experiment("2016-07-23", "e-coli, 7 generations", timestep=1)
	6
	"""
	
	cnx = _connect()
	cursor=cnx.cursor()
	
	cursor.execute("insert into experiment(date,comments, timestep, fluotimestep, scale) values('{}','{}', {}, {}, {})".format(date, comments, timestep, fluotimestep, scale))

	expid = cursor.lastrowid

	cnx.commit()
	cursor.close()
	cnx.close()

	return expid

def delete_experiment(expid):

	"""Permanently removes an experiment from the database.

	Parameters
	----------
	expid: int
		The id of the experiment to delete

	Examples
	--------
	>>> delete_experiment(3)
	"""
	
	cnx = _connect()
	cursor=cnx.cursor()
	
	cursor.execute("delete from experiment where id = {}".format(expid))

	cnx.commit()
	cursor.close()
	cnx.close()


def _upload_old_old(experimentid, cellfile, fluofile=None, last_division=True):

	"""Upload new data to the database.

	Parameters
	----------

	experimentid : int
		The experiment id.
	
	cellfile : string
		The path to the csv file containing the cell measurements.

	fluofile : string
		The path to the csv file containing the fluorescence measurments.

	last_division : boolean
		Whether the last cells were recorded until they divided.
	"""

	cnx = _connect()
	cursor=cnx.cursor()

	try:
		#retrieving timestep and fluotimestep from experiment
		cursor.execute("select timestep, fluotimestep from experiment where experiment.id={}".format(experimentid))
		res = cursor.fetchone()
		if res == None:
			raise Exception("Experiment {} could not be found".format(experimentid))
		else:
			timestep,fluotimestep = res


		cells=np.zeros(5096, dtype='int') #assuming we won't have more than 5000 cells in one import. This might need a fix.
		lineageids=[]

		with open(cellfile) as datafile:
			reader = csv.reader(datafile)
			#lines=datafile.readlines()
			reader.next()		#first line is title of columns
			#for line in lines[1:]:
			for l in reader:
				#l=line.split(',')
				time=(int(l[0])-1)*timestep #converting frame to minutes
#				time="{:02d}:{:02d}:{:02d}".format(int(time/3600),int(time%3600/60),int(time%60)) #converting seconds to SQL TIME format

				if cells[int(l[1])]==0: #new cell detected

					# cell i gives birth to cell 2i and 2i+1
					# the cell inheriting the oldpole has the same parity as it's parent
					parent=cells[int(l[1])/2]

					#other numbering system where cell 'i' gives birth to cell 'i1' and 'i2'
					#parent=cells[int(l[1])/10]
				
					#if no parent, create a new lineage
					if parent==0:
						cursor.execute("insert into lineage (experimentid, comments) values ({}, '{}')".format(experimentid, path.basename(cellfile)))
						lineageid = cursor.lastrowid
						lineageids.append(lineageid)
						#position in classic tree, i.e. tree with old pole on the outside
						side="NULL"
						outtreeposition = 1
						pole="NULL"
						lefttreeposition = 1
					else:
						cursor.execute("update cell set divisiontime='{}' where id = {}".format(time,parent))
						cursor.execute("select lineageid, outtreeposition, side, pole, lefttreeposition from cell where id={}".format(parent))
						lineageid, parentouttreeposition, parentside, parentpole, parentlefttreeposition = cursor.fetchone()
						side = int(int(l[1])%2)
						outtreeposition = parentouttreeposition*2 + side
						if outtreeposition <= 3: 
							pole = "NULL"
							lefttreeposition = outtreeposition
						else:
							pole = 1-int(side%2 == parentside%2)
							lefttreeposition = parentlefttreeposition*2 + pole
				
					#create new cell
					cursor.execute("insert into cell (lineageid, parentid, outtreeposition, side, pole, lefttreeposition, birthtime) values ({},{},{},{},{},{},{})".format(lineageid, parent, outtreeposition, side, pole, lefttreeposition, time if outtreeposition>1 else  "NULL"))
					cid=cursor.lastrowid
					cells[int(l[1])]=cid

				# cell already exists
				else:
					#retrieve cell id
					cid=cells[int(l[1])]

				#add new state:
				cursor.execute("insert into state (cellid, time, length, maxwidth, area, volume) values ({},{},{},{},{},{})".format(cid, time, l[2], l[5], l[3], l[4]))
				sid = cursor.lastrowid

				#if fluorescence was measured:
				if len(l)>6 and l[6] != "" and time%fluotimestep == 0:
					cursor.execute("update state set intensity = {} where id={}".format(l[6], sid))
				#if fluofile != None:
				#	pass
				#if heatshocked:
				if len(l)>8 and l[8] != "":
					cursor.execute("update cell set alive='{}' where id={}".format(l[8]=="1",cid))
					#lagtime of alive cells

					if l[8]=='1' and len(l)>9 and l[9] != "":
						cursor.execute("update cell set lagtime={} where id={}".format(float(l[9]),cid))

		#division of last cells
		if last_division:
			if _dbtype == 'sqlite':
				cursor.execute("select c.id, s.maxtime from cell c join (select cellid, max(time) as maxtime from state group by cellid) as s on c.id=s.cellid where c.lineageid in ({}) and c.divisiontime is null".format(','.join(map(str,lineageids))))
				for cid,maxtime in cursor.fetchall():
					cursor.execute("update cell set divisiontime = {} where id = {}".format(maxtime + timestep, cid))
			elif _dbtype == 'mysql':
				cursor.execute("update cell c join (select cellid, max(time) as maxtime from state group by cellid) as s on c.id=s.cellid set c.divisiontime=s.maxtime+{} where c.lineageid in ({}) and c.divisiontime is null".format(timestep, ','.join(map(str,lineageids))))

		if fluofile is not None:
			#TODO
			with open(fluofile) as datafile:
				reader = csv.reader(datafile)
				#lines=datafile.readlines()
				for i in range(4):
					reader.next()
				#for line in lines[4:]: #start at 4th line. why? might be to change later
				for l in reader:
					#l=line.split(',')
					# print l[4][:-1]
					time=(int(l[0])-1)*timestep # convert frame to time in minute
					#time="%02d:%02d:%02d"%(t/3600,t%3600/60,t%60)
					#cursor.execute("select id from cell c where c.lineageid in {} and positioninclassictree={}"%(','.join(map(str,lineageids)), l[1]))
					#cellid=cursor.fetchall()[0][0]
					cursor.execute("update state set nucleoids={}, relativearea={} where cellid={} and time='{}'".format(int(l[2]), float(l[4].replace('\n','')), cells[int(l[1])], time))

		cnx.commit()

	except Exception as ex:
		print traceback.format_exc()
		cnx.rollback()

	finally:
		cursor.close()
		cnx.close()


def _upload_old(experimentid, cellfile, last_division=True):

	"""Upload new data to the database.
	
	Parameters
	----------

	experimentid : int
		The experiment id.
	
	datafile : string
		The path to the csv file containing the data.
		The file must be a csv file (comma separated).
		The first line contains keywords explaining what each column contains.
		Their order does not matter.
		Possible keywords are the following (those marked with a * are mandatory).
		*frame: integer
			Frame number. Starts from 1.
		*cell: integer
			Cell number. The numbering should be based on the parent-child relation.
			Up to now, only one numbering system is supported: the daughters of cell i are numbered 2*i for the cell inheriting the old pole, 2*i+1 for the cell inheriting the new pole. If a new cell number n is encountered while n/2 does not exist, it is considered as the mother cell of a new lineage.
		length: real
			The length of the cell in pixels.
		area: real
			The area of the cell in squarepixels.
		volume: real
			The volume of the cell in cubicpixels.
		maxwidth: real
			The maximum width of the cell in pixels.
		x: real
			The x coordinate of the center of the cell.
		y: real
			The y coordinate of the center of the cell.
		intensity: real
			The intensity of a fluorescence signal in the whole cell. Unit is arbitrary.
		nucleoids: integer
			The number of fluorescent nucleoid in the cell.
		relativearea: real in [0,1]
			The fraction of the cell area occupied by the nucleoids.
		alive: {0,1}
			In case the cells was exposed to a (heat) shock treatment, whether the cell survived that shock. 0 means it died, 1 means means it's alive.
		
		Not implemented yet: if several fluo signals, allow intensity1, intensity2, nucledoids1, nucleoids2,...

		The rows have to be ordered by frame, but the order within one frame does not matter.

		Example:
			frame,cell,length,maxwidth,alive
			1,1,40.25,8.8,
			2,1,42.53,9.1,
			3,2,22.56,8.6,
			3,3,21.19,9.1,
			3,3,22.86,9.2,1
			3,2,24.03,8.6,0

	last_division : boolean
		Whether the last recorded cells divided at the end.
	"""

	cnx = _connect()
	cursor=cnx.cursor()

	try:
		#retrieving timestep and fluotimestep from experiment
		cursor.execute("select timestep, fluotimestep from experiment where experiment.id={}".format(experimentid))
		res = cursor.fetchone()
		if res == None:
			raise Exception("Experiment {} could not be found".format(experimentid))
		else:
			timestep,fluotimestep = res


		with open(cellfile) as datafile:
			reader = csv.reader(datafile)

			#Checking columns and getting order.
			columns = reader.next()
			validcolumns=["frame", "cell", "length", "area", "volume", "maxwidth", "x", "y", "intensity", "nucleoids", "relativearea", "alive", "lagtime"]
			for c in columns:
				if c not in validcolumns:
					print("Unknown column name: {}".format(c))
					print("Column will be ignored")
			idx = {}
			for i,c in enumerate(columns):
				idx[c]=i

			#########Building dataset############
			cellids = {}
			lineageids=[]

			#Reading the data
			for row in reader:
				time=(int(row[idx["frame"]])-1)*timestep #converting frame to minutes
				
				cell = int(row[idx["cell"]])
				if cell <= 0:
					raise Exception("Line {}: Cell index has to be > 0.".format(reader.line_num))

				#print cellids
				if not cellids.has_key( cell ): #new cell detected

					# cell i gives birth to cell 2i and 2i+1
					# the cell inheriting the oldpole has the same parity as its parent
					# i.e. oldpole on the outside
					parent = cell/2

					#other numbering system where cell 'i' gives birth to cell 'i1' and 'i2'
					#parent=cell/10
				
					#if no parent, create a new lineage
					if not cellids.has_key(parent):
						parentid=0
						cursor.execute("insert into lineage (experimentid, comments) values ({}, '{}')".format(experimentid, cellfile))
						lineageid = cursor.lastrowid
						lineageids.append(lineageid)
						side="NULL"
						outtreeposition = 1
						pole="NULL"
						lefttreeposition = 1
					else:
						parentid=cellids[parent]
						cursor.execute("update cell set divisiontime='{}' where id = {}".format(time,parentid))
						cursor.execute("select lineageid, outtreeposition, side, pole, lefttreeposition from cell where id={}".format(parentid))
						lineageid, parentouttreeposition, parentside, parentpole, parentlefttreeposition = cursor.fetchone()
						
						side = cell%2
						if parentouttreeposition is None or type(parentouttreeposition*2) == "long":
							outtreeposition = "NULL"
							lefttreeposition = "NULL"
							pole = 1-int(side%2 == parentside%2)
						else:
							outtreeposition = parentouttreeposition*2 + side
							if outtreeposition <= 3:
								pole = "NULL"
								lefttreeposition = outtreeposition
							else:
								pole = 1-int(side%2 == parentside%2)
								lefttreeposition = parentlefttreeposition*2 + pole
				
					#create new cell
					cursor.execute("insert into cell (lineageid, parentid, outtreeposition, side, pole, lefttreeposition, birthtime) values ({},{},{},{},{},{},{})".format(lineageid, parentid, outtreeposition, side, pole, lefttreeposition, time if outtreeposition>1 else  "NULL"))
					cellid=cursor.lastrowid
					cellids[cell]=cellid

				# cell already exists
				else:
					#retrieve cell id
					cellid=cellids[cell]

				#add new state:
				cursor.execute("insert into state (cellid, time) values ({},{})".format(cellid, time))
				sid = cursor.lastrowid

				for attr in ["length", "area", "volume", "maxwidth", "x", "y", "intensity"]:
					if idx.has_key(attr):
						cursor.execute("update state set {} = {} where id={}".format(attr,row[idx[attr]], sid))

				#if heatshocked:
				if idx.has_key("alive") and row[idx["alive"]] != "":
					cursor.execute("update cell set alive='{}' where id={}".format(row[idx["alive"]]=="1",cellid))
				
					#lagtime of alive cells
					if idx.has_key("lagtime") and row[idx["lagtime"]] != "":
						cursor.execute("update cell set lagtime={} where id={}".format(float(row[idx["lagtime"]]),cellid))

		#division of last cells
		if last_division:
			if _dbtype == 'sqlite':
				cursor.execute("select c.id, s.maxtime from cell c join (select cellid, max(time) as maxtime from state group by cellid) as s on c.id=s.cellid where c.lineageid in ({}) and c.divisiontime is null".format(','.join(map(str,lineageids))))
				for cid,maxtime in cursor.fetchall():
					cursor.execute("update cell set divisiontime = {} where id = {}".format(maxtime + timestep, cid))
			elif _dbtype == 'mysql':
				cursor.execute("update cell c join (select cellid, max(time) as maxtime from state group by cellid) as s on c.id=s.cellid set c.divisiontime=s.maxtime+{} where c.lineageid in ({}) and c.divisiontime is null".format(timestep, ','.join(map(str,lineageids))))

#		if fluofile is not None:
#			#TODO: to rewrite for new input format!!!!!!
#			with open(fluofile) as datafile:
#				reader = csv.reader(datafile)
#				#lines=datafile.readlines()
#				for i in range(4):
#					reader.next()
#				#for line in lines[4:]: #start at 4th line. why? might be to change later
#				for l in reader:
#					#l=line.split(',')
#					# print l[4][:-1]
#					time=(int(l[0])-1)*timestep # convert frame to time in minute
#					#time="%02d:%02d:%02d"%(t/3600,t%3600/60,t%60)
#					#cursor.execute("select id from cell c where c.lineageid in {} and positioninclassictree={}"%(','.join(map(str,lineageids)), l[1]))
#					#cellid=cursor.fetchall()[0][0]
#					cursor.execute("update state set nucleoids={}, relativearea={} where cellid={} and time='{}'".format(int(l[2]), float(l[4].replace('\n','')), cells[int(l[1])], time))

		cnx.commit()

	except Exception as ex:
		if "row" in locals():
			print row
		print traceback.format_exc()
		cnx.rollback()

	finally:
		cursor.close()
		cnx.close()


def upload(experimentid, cellfile, last_division=True):

	"""Upload new data to the database.
	
	Parameters
	----------

	experimentid : int
		The experiment id.
	
	datafile : string
		The path to the csv file containing the data.
		The file must be a csv file (comma separated).
		The first line contains keywords explaining what each column contains.
		Their order does not matter.
		Possible keywords are the following (those marked with a * are mandatory).
		*frame: integer
			Frame number. Starts from 1.
		*cell: integer
			Cell number. The numbering should be based on the parent-child relation.
			Up to now, only one numbering system is supported: the daughters of cell i are numbered 2*i for the cell inheriting the old pole, 2*i+1 for the cell inheriting the new pole. If a new cell number n is encountered while n/2 does not exist, it is considered as the mother cell of a new lineage.
		length: real
			The length of the cell in pixels.
		area: real
			The area of the cell in squarepixels.
		volume: real
			The volume of the cell in cubicpixels.
		maxwidth: real
			The maximum width of the cell in pixels.
		x: real
			The x coordinate of the center of the cell.
		y: real
			The y coordinate of the center of the cell.
		intensity: real
			The intensity of a fluorescence signal in the whole cell. Unit is arbitrary.
		nucleoids: integer
			The number of fluorescent nucleoid in the cell.
		relativearea: real in [0,1]
			The fraction of the cell area occupied by the nucleoids.
		alive: {0,1}
			In case the cells was exposed to a (heat) shock treatment, whether the cell survived that shock. 0 means it died, 1 means means it's alive.
		
		Not implemented yet: if several fluo signals, allow intensity1, intensity2, nucledoids1, nucleoids2,...

		The rows have to be ordered by frame, but the order within one frame does not matter.

		Example:
			frame,cell,length,maxwidth,alive
			1,1,40.25,8.8,
			2,1,42.53,9.1,
			3,2,22.56,8.6,
			3,3,21.19,9.1,
			3,3,22.86,9.2,1
			3,2,24.03,8.6,0

	last_division : boolean
		Whether the last recorded cells divided at the end.
	"""

	cnx = _connect()
	cursor=cnx.cursor()

	try:
		#retrieving timestep and fluotimestep from experiment
		cursor.execute("select timestep, fluotimestep from experiment where experiment.id={}".format(experimentid))
		res = cursor.fetchone()
		if res == None:
			raise Exception("Experiment {} could not be found".format(experimentid))
		else:
			timestep,fluotimestep = res


		with open(cellfile) as datafile:
			reader = csv.reader(datafile)

			#Checking columns and getting order.
			columns = reader.next()
			knowncolumns=["frame", "cell", "length", "area", "volume", "maxwidth", "x", "y", "intensity", "nucleoids", "relativearea", "alive", "lagtime"]
			newattributes={}
			newattributestypes={}
			newattributesobjects={}
			extraknowncolumns = {}
			for n,t,o in cursor.execute("SELECT name,type,object FROM attributes WHERE experimentid = {}".format(experimentid)).fetchall():
				extraknowncolumns[n]=(t,o)
			for c in columns:
				if c == "":
					continue
				while c not in knowncolumns:
					if c in extraknowncolumns:
						knowncolumns.append(c)
						newattributes[c]={}
						newattributestypes[c]=extraknowncolumns[c][0]
						newattributesobjects[c]="cell" if extraknowncolumns[c][1]== 0 else "state"
					else:
						print("Unknown column name: {}. A new attribute will be created.".format(c))
						answer = raw_input("Proceed? (y/n)")
						if answer == "y":
							knowncolumns.append(c)
							newattributes[c]={}
							newattributestypes[c]="BOOL"
							newattributesobjects[c]="cell"
							#TODO: test if attribute already exists in db
						elif answer == "n":
							return
						else:
							continue
			idx = {}
			for i,c in enumerate(columns):
				idx[c]=i

			#########Building dataset############
			cellids = {}
			lineageids=[]

			#Reading the data
			for row in reader:
				time=(int(row[idx["frame"]])-1)*timestep #converting frame to minutes
				
				cell = int(row[idx["cell"]])
				if cell <= 0:
					raise Exception("Line {}: Cell index has to be > 0.".format(reader.line_num))

				#print cellids
				if not cellids.has_key( cell ): #new cell detected

					# cell i gives birth to cell 2i and 2i+1
					# the cell inheriting the oldpole has the same parity as its parent
					# i.e. oldpole on the outside
					parent = cell/2

					#other numbering system where cell 'i' gives birth to cell 'i1' and 'i2'
					#parent=cell/10
				
					#if no parent, create a new lineage
					if not cellids.has_key(parent):
						parentid=0
						cursor.execute("insert into lineage (experimentid, comments) values ({}, '{}')".format(experimentid, cellfile))
						lineageid = cursor.lastrowid
						lineageids.append(lineageid)
						side="NULL"
						outtreeposition = 1
						pole="NULL"
						lefttreeposition = 1
					else:
						parentid=cellids[parent]
						cursor.execute("update cell set divisiontime='{}' where id = {}".format(time,parentid))
						cursor.execute("select lineageid, outtreeposition, side, pole, lefttreeposition from cell where id={}".format(parentid))
						lineageid, parentouttreeposition, parentside, parentpole, parentlefttreeposition = cursor.fetchone()
						
						side = cell%2
						if parentouttreeposition is None or type(parentouttreeposition*2) == "long":
							outtreeposition = "NULL"
							lefttreeposition = "NULL"
							pole = 1-int(side%2 == parentside%2)
						else:
							outtreeposition = parentouttreeposition*2 + side
							if outtreeposition <= 3:
								pole = "NULL"
								lefttreeposition = outtreeposition
							else:
								pole = 1-int(side%2 == parentside%2)
								lefttreeposition = parentlefttreeposition*2 + pole
				
					#create new cell
					cursor.execute("insert into cell (lineageid, parentid, outtreeposition, side, pole, lefttreeposition, birthtime) values ({},{},{},{},{},{},{})".format(lineageid, parentid, outtreeposition, side, pole, lefttreeposition, time if outtreeposition>1 else  "NULL"))
					cellid=cursor.lastrowid
					cellids[cell]=cellid

				# cell already exists
				else:
					#retrieve cell id
					cellid=cellids[cell]

				#add new state:
				cursor.execute("insert into state (cellid, time) values ({},{})".format(cellid, time))
				sid = cursor.lastrowid

				for attr in ["length", "area", "volume", "maxwidth", "x", "y", "intensity"]:
					if idx.has_key(attr) and row[idx[attr]] != "":
						cursor.execute("update state set {} = {} where id={}".format(attr,row[idx[attr]], sid))
						
				#new cell attributes
				for attr in newattributes:
					if row[idx[attr]] != "":
						if newattributes[attr].has_key(cellid):
							if newattributesobjects[attr] == "cell":
								newattributesobjects[attr] = "state"
							newattributes[attr][cellid][time]=row[idx[attr]]
						else:
								newattributes[attr][cellid]={time:row[idx[attr]]}
							
						#infer type
						if newattributestypes[attr]=="BOOL" and not is_boolean(row[idx[attr]]):
							newattributestypes[attr]="DOUBLE"
						if newattributestypes[attr]=="DOUBLE" and not is_float(row[idx[attr]]):
							newattributestypes[attr]="TEXT"
						
		for aname in newattributes:
			attribute = newattributes[aname]
			
			if newattributesobjects[aname] == "cell":
				if not np.any([r[1]==aname for r in cursor.execute("PRAGMA table_info ('cell')").fetchall()]):
					cursor.execute("ALTER TABLE cell ADD COLUMN {} {} DEFAULT NULL".format(aname,newattributestypes[aname]))
				for c in attribute:
					values = [attribute[c][t] for t in attribute[c]]
					cursor.execute("UPDATE cell SET {} = '{}' WHERE id = {}".format(aname, values[0], c))  #TODO more efficiently for BOOL and TEXT (when similiar values)
				
				cursor.execute("INSERT INTO attributes (experimentid, name, type, object) values ({},'{}','{}',0)".format(experimentid, aname, newattributestypes[aname]))
					
			else: # new state attribute
			
				#make sure all values are same format (as dict with time as key)
				for c in attribute:
					if type(attribute[c]) != dict:
						value,t = attribute[cellid]
						attribute[cellid]={t:value}
					#else: #TO CHECK
						#break
				
				#update data
				if not np.any([r[1]==aname for r in cursor.execute("PRAGMA table_info ('state')").fetchall()]):
					cursor.execute("ALTER TABLE state ADD COLUMN {} {} DEFAULT NULL".format(aname,newattributestypes[aname]))
				for c in attribute:
					for t in attribute[c]:
						cursor.execute("UPDATE state SET {} = '{}' WHERE cellid = {} AND time = {}".format(aname, attribute[c][t], c, t))    #TODO more efficiently for BOOL and TEXT (when similiar values)
					
				cursor.execute("INSERT INTO attributes (experimentid, name, type, object) values ({},'{}','{}',1)".format(experimentid, aname, newattributestypes[aname]))
					
		#division of last cells
		if last_division:
			if _dbtype == 'sqlite':
				cursor.execute("select c.id, s.maxtime from cell c join (select cellid, max(time) as maxtime from state group by cellid) as s on c.id=s.cellid where c.lineageid in ({}) and c.divisiontime is null".format(','.join(map(str,lineageids))))
				for cid,maxtime in cursor.fetchall():
					cursor.execute("update cell set divisiontime = {} where id = {}".format(maxtime + timestep, cid))
			elif _dbtype == 'mysql':
				cursor.execute("update cell c join (select cellid, max(time) as maxtime from state group by cellid) as s on c.id=s.cellid set c.divisiontime=s.maxtime+{} where c.lineageid in ({}) and c.divisiontime is null".format(timestep, ','.join(map(str,lineageids))))

		cnx.commit()

	except Exception as ex:
		if "row" in locals():
			print row
		print traceback.format_exc()
		cnx.rollback()

	finally:
		cursor.close()
		cnx.close()


		
def load(experimentids=[]):

	"""Loads data in the memory.

	Loads the given experiments as array of Cell and Lineage objects.
	
	Parameters
	----------
	experimentids : list
		The list of id of the experiments to load.

	Returns
	-------
	cells : array
		numpy array of Cell objects.
	lineages : array
		numpy array of Lineage objects.

	Examples
	--------
	>>> cells, lineages = load([2,4])
	"""

	assert len(experimentids)>0, "No experiment specified."

	#To not reload what has just been loaded
	global _loaded
	global _cells
	global _lineages
	if _loaded==experimentids:
		return _cells,_lineages
		
	
	print("loading...")
	start = time.time()

	cnx = _connect()
	cursor=cnx.cursor()
	
	
	######Retrieving the custom attributes for the required experiments
	#print "getting attributes...",
	cursor.execute("SELECT name,type,object FROM attributes WHERE experimentid in ({})".format(",".join(map(str,experimentids))))
	extracellattributes=[]
	extrastateattributes=[]
	for name,type,object in cursor.fetchall():
		if object == 0: #cell attribute
			extracellattributes.append(name)
			data.Cell.attribute_types[name]=_sqlitedatatypes[type]
		else: #object==1, state attribute
			extrastateattributes.append(name)
			data.Cell.attribute_types[name]="list.{}".format(_sqlitedatatypes[type])
			
	#print "{}".format(time.time()-start)
			
	#Retrieving lineages id
	#print "getting lineages...",
	start = time.time()
	
	cursor.execute("select l.id from lineage l where experimentid in ({})".format(",".join(map(str,experimentids))))
	lineageids=[l[0] for l in cursor.fetchall()]

	#print "{}".format(time.time()-start)

	######Retrieving the cells
	#print "getting cells...",
	#start = time.time()
	
	#TODO: load by lineages to go faster?
	
	cursor.execute("select c.id, c.lineageid, c.pole, c.lefttreeposition, c.birthtime, c.divisiontime, c.alive, c.lagtime{extra} from cell c where c.lineageid in ({ids}) order by c.id".format(extra="".join([", c.{}".format(attr) for attr in extracellattributes]), ids=','.join([str(s) for s in lineageids])))
	rows=cursor.fetchall()
	cells=[]
	cellids = []
	
	#print "{}".format(time.time()-start)
	#print "getting states and building cell objects...",
	start = time.time()

	#SLOW!!!!
	for r in rows:
		cursor.execute("select s.time, s.length, s.maxwidth, s.area, s.intensity, s.nucleoids, s.relativearea {extra} from state s where s.cellid={cid} order by s.time".format(cid=r[0], extra="".join([", s.{}".format(attr) for attr in extrastateattributes])))
		states = np.array(cursor.fetchall())
		stime=states[:,0]
		length=states[:,1]
		maxwidth=states[:,2]
		area=states[:,3]
		intensity=states[:,4]
		nucleoids=states[:,5]
		relativearea=states[:,6]

		cellids.append(r[0])
		
		cells.append(data.Cell(r[0], r[1], None, r[2], r[3], r[4] if r[4] is not None else None, r[5] if r[5] is not None else None, time=stime, length=length, maxwidth=maxwidth, area=area, fluorescence=intensity, nucleoids=nucleoids, relativearea=relativearea))
		if r[6] != None:
			cells[-1].alive=1 if r[6]=='True' else 0
		if r[7] != None:
			cells[-1].lagtime=r[7]
		a=8
		for attr in extracellattributes:
			if r[a] != None:
				setattr(cells[-1], attr, r[a])
			a=a+1
		a=7
		for attr in extrastateattributes:
			setattr(cells[-1], attr, states[:,a])
			a=a+1
		
	
	#print "{}".format(time.time()-start)
	
	#print "building parent-child relations...",
	start = time.time()

	_cells=np.array(cells)
	cursor.execute("select c.id, c.parentid from cell c where c.lineageid in ({}) order by c.id".format(','.join([str(s) for s in lineageids])))
	rows=cursor.fetchall()
	for i,(cell,r) in enumerate(zip(_cells,rows)):
		if r[1] is None or r[1]==0:
			continue
		parent=_cells[cellids.index(r[1])]
		cell.parent=parent
		parent.children.append(cell)
		if cell.lefttreeposition is None and cell.pole is not None:
			cell.lefttreeposition = parent.lefttreeposition*2+cell.pole
			cell.generation=parent.generation+1
			cell.age = parent.age+1 if cell.pole == 0 else 0
			cell.youth = parent.youth+1 if cell.youth == 0 else 1
	
	#print "{}".format(time.time()-start)
	
	#print "building lineage objects...",
	start = time.time()

	lineages = []
	for lid in lineageids:
		cursor.execute("select experiment.id,timestep, fluotimestep from experiment join lineage on lineage.experimentid=experiment.id where lineage.id={}".format(lid))
		eid,timestep,fluotimestep = cursor.fetchone()
		lineages.append(data.Lineage(lid, 0, timestep, fluotimestep, data.select(_cells,"lineageid == {}".format(lid))))
	_lineages = np.array(lineages)

	#print "{}".format(time.time()-start)
	
	print("done.")

	cnx.close()
	
	_loaded=experimentids
	
	return _cells, _lineages
	