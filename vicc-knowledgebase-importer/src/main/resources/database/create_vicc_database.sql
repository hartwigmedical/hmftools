SET FOREIGN_KEY_CHECKS = 0;

DROP TABLE IF EXISTS viccEntry;
CREATE TABLE viccEntry
(   id int NOT NULL AUTO_INCREMENT,
    source varchar(255) NOT NULL,
    PRIMARY KEY (id)
);

DROP TABLE IF EXISTS gene;
CREATE TABLE gene
(   id int NOT NULL AUTO_INCREMENT,
    viccEntryId int NOT NULL,
    geneName varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (viccEntryId) REFERENCES viccEntry(id)
);

DROP TABLE IF EXISTS geneIdentifier;
CREATE TABLE geneIdentifier
(   id int NOT NULL AUTO_INCREMENT,
    viccEntryId int NOT NULL,
    symbol varchar(255) NOT NULL,
    entrezId varchar(255) NOT NULL,
    ensemblGeneId varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (viccEntryId) REFERENCES viccEntry(id)
);

DROP TABLE IF EXISTS featureName;
CREATE TABLE featureName
(   id int NOT NULL AUTO_INCREMENT,
    viccEntryId int NOT NULL,
    featureName varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (viccEntryId) REFERENCES viccEntry(id)
);

DROP TABLE IF EXISTS tag;
CREATE TABLE tag
(   id int NOT NULL AUTO_INCREMENT,
    viccEntryId int NOT NULL,
    tagName varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (viccEntryId) REFERENCES viccEntry(id)
);

DROP TABLE IF EXISTS devTag;
CREATE TABLE devTag
(   id int NOT NULL AUTO_INCREMENT,
    viccEntryId int NOT NULL,
    devTagName varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (viccEntryId) REFERENCES viccEntry(id)
);

SET FOREIGN_KEY_CHECKS = 1;