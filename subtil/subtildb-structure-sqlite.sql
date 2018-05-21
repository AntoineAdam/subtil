
--
-- Database : `subtildb`
--

-- --------------------------------------------------------

--
-- Table `experiment`
--

CREATE  TABLE IF NOT EXISTS "Experiment" (
  "Id" INTEGER PRIMARY KEY  AUTOINCREMENT  NOT NULL ,
  "Date" DATETIME DEFAULT CURRENT_TIMESTAMP,
  "Comments" TEXT, 
  "Timestep" FLOAT NOT NULL ,
  "Fluotimestep" FLOAT DEFAULT NULL,
  "Scale" FLOAT NOT NULL
) ;
-- --------------------------------------------------------

--
-- Table `attributes`
--

CREATE  TABLE IF NOT EXISTS "Attributes" (
  "Id" INTEGER PRIMARY KEY  AUTOINCREMENT  NOT NULL ,
  "ExperimentId" INTEGER NOT NULL ,
  "Name" TEXT NOT NULL , 
  "Object" INTEGER NOT NULL ,
  "Type" TEXT NOT NULL, 
  FOREIGN KEY("ExperimentId") REFERENCES Experiment("Id") ON DELETE CASCADE ON UPDATE CASCADE
) ;
-- --------------------------------------------------------

--
-- Table `lineage`
--

CREATE  TABLE IF NOT EXISTS "Lineage" (
  "Id" INTEGER PRIMARY KEY  AUTOINCREMENT  NOT NULL ,
  "ExperimentId" INTEGER NOT NULL ,
  "Comments" TEXT,
  FOREIGN KEY("ExperimentId") REFERENCES Experiment("Id") ON DELETE CASCADE ON UPDATE CASCADE
) ;

-- --------------------------------------------------------

--
-- Table `cell`
--

CREATE TABLE IF NOT EXISTS "Cell" (
  "Id" INTEGER PRIMARY KEY  AUTOINCREMENT  NOT NULL,
  "LineageId" INTEGER NOT NULL,
  "ParentId" INTEGER DEFAULT NULL,
  "Pole" BOOL DEFAULT NULL,
--  "PositionInAgeTree" INTEGER NOT NULL DEFAULT 1,
  "LeftTreePosition" INTEGER NOT NULL DEFAULT 1,
  "Side" BOOL DEFAULT NULL,
--   "PositionInClassicTree" INTEGER NOT NULL DEFAULT 1,
  "OutTreePosition" INTEGER NOT NULL DEFAULT 1,
  "BirthTime" FLOAT DEFAULT NULL,
  "DivisionTime" FLOAT DEFAULT NULL,
  "LagTime" FLOAT DEFAULT NULL,
  "Alive" BOOL DEFAULT NULL,
  FOREIGN KEY("LineageId") REFERENCES Lineage("Id") ON DELETE CASCADE ON UPDATE CASCADE
--  FOREIGN KEY("ParentId") REFERENCES Cell("Id")
) ;

CREATE UNIQUE INDEX "cellindex" ON "Cell" ("Id" ASC) ;

-- --------------------------------------------------------

--
-- Table `state`
--

CREATE TABLE IF NOT EXISTS "State" (
  "Id" INTEGER PRIMARY KEY  AUTOINCREMENT  NOT NULL,
  "CellId" INTEGER NOT NULL,
  "Time" FLOAT DEFAULT NULL,
  "MaxWidth" DOUBLE DEFAULT NULL,
  "Length" DOUBLE DEFAULT NULL,
--  "Curvature" DOUBLE DEFAULT NULL,
--  "Perimeter" DOUBLE DEFAULT NULL,
  "Area" DOUBLE DEFAULT NULL,
  "Volume" DOUBLE DEFAULT NULL,
  "X" DOUBLE DEFAULT NULL,
  "Y" DOUBLE DEFAULT NULL,
  "Intensity" DOUBLE DEFAULT NULL,
  "Nucleoids" INTEGER DEFAULT NULL,
  "RelativeArea" DOUBLE DEFAULT NULL,
  FOREIGN KEY("CellId") REFERENCES Cell("Id") ON DELETE CASCADE ON UPDATE CASCADE
) ;

CREATE UNIQUE INDEX "stateindex" ON "State" ("Id" ASC) ;
CREATE UNIQUE INDEX "statecellindex" ON "State" ("CellId" ASC, "Time" ASC) ;

