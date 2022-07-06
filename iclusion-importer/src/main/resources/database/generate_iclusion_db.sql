SET FOREIGN_KEY_CHECKS = 0;

DROP TABLE IF EXISTS tumorLocations;
DROP TABLE IF EXISTS blacklistedTumorLocations;
DROP TABLE IF EXISTS mutationConditions;

DROP TABLE IF EXISTS study;
CREATE TABLE study
(   id int NOT NULL AUTO_INCREMENT,
    idDB varchar(50) NOT NULL,
    acronym varchar(100) NOT NULL,
    title varchar(500) NOT NULL,
    eudra varchar(50),
    nct varchar(50),
    ipn varchar(50),
    ccmo varchar(50),
    PRIMARY KEY (id)
);

DROP TABLE IF EXISTS tumorLocation;
CREATE TABLE tumorLocation
(   id int NOT NULL AUTO_INCREMENT,
    tumorLocationId int NOT NULL,
    tumorLocation varchar(50) NOT NULL,
    FOREIGN KEY (tumorLocationId) REFERENCES study(id),
    PRIMARY KEY (id)
);

DROP TABLE IF EXISTS blacklistedTumorLocation;
CREATE TABLE blacklistedTumorLocation
(   id int NOT NULL AUTO_INCREMENT,
    blacklistedTumorLocationId int NOT NULL,
    blacklistedTumorLocation varchar(50) NOT NULL,
    FOREIGN KEY (blacklistedTumorLocationId) REFERENCES study(id),
    PRIMARY KEY (id)
);

DROP TABLE IF EXISTS mutationCondition;
CREATE TABLE mutationCondition
(   id int NOT NULL AUTO_INCREMENT,
    mutationConditionId int NOT NULL,
    gene varchar(50) NOT NULL,
    mutation varchar(50) NOT NULL,
    FOREIGN KEY (mutationConditionId) REFERENCES study(id),
    PRIMARY KEY (id)
);

SET FOREIGN_KEY_CHECKS = 1;