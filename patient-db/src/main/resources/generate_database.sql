CREATE TABLE patientInfo
 ( id int NOT NULL AUTO_INCREMENT, cpctId varchar(255) DEFAULT NULL, drupId varchar(255) DEFAULT NULL, sex varchar(10) DEFAULT NULL, birthYear int DEFAULT NULL, hospital varchar(255) DEFAULT NULL, ethnicity varchar(255) DEFAULT NULL, deathDate DATE, PRIMARY KEY (id) );

CREATE TABLE tumorData
 ( id int NOT NULL AUTO_INCREMENT, location varchar(255) DEFAULT NULL, entryStage varchar(255) DEFAULT NULL, PRIMARY KEY (id), patientId int NOT NULL, FOREIGN KEY (patientId) REFERENCES patientInfo(id) );

CREATE TABLE biopsyLocations
 ( id int NOT NULL AUTO_INCREMENT, location varchar(255) DEFAULT NULL, tumorId int NOT NULL, PRIMARY KEY (id), FOREIGN KEY (tumorId) REFERENCES tumorData(id) );

CREATE TABLE systemicTherapyData
 ( id int NOT NULL AUTO_INCREMENT, startDate DATE, endDate DATE, type varchar(255), treatment varchar(255), bestResponse varchar(255), PRIMARY KEY (id), patientId int NOT NULL, FOREIGN KEY (patientId) REFERENCES patientInfo(id) );

CREATE TABLE radioTherapyData
 ( id int NOT NULL AUTO_INCREMENT, endDate DATE, site varchar(255), patientId int NOT NULL, PRIMARY KEY (id), FOREIGN KEY (patientId) references patientInfo(id) );

CREATE TABLE treatmentData
 ( id int NOT NULL AUTO_INCREMENT, startDate DATE, endDate DATE, name varchar(255), earlyResponse varchar(255), radiotherapyStartDate DATE, radiotherapyEndDate DATE, patientId int NOT NULL, PRIMARY KEY (id), FOREIGN KEY (patientId) references patientInfo(id) );