ALTER TABLE copyNumber CHANGE structuralVariantSupport segmentStartSupport varchar(255) NOT NULL;
ALTER TABLE copyNumber ADD COLUMN segmentEndSupport varchar(255) NOT NULL AFTER segmentStartSupport;
ALTER TABLE copyNumberRegion CHANGE structuralVariantSupport segmentStartSupport varchar(255) NOT NULL;
ALTER TABLE copyNumberRegion ADD COLUMN actualTumorBaf DOUBLE PRECISION not null AFTER modelTumorRatio;

ALTER TABLE copyNumber DROP COLUMN ratioSupport;
ALTER TABLE copyNumber ADD COLUMN inferred BOOLEAN NOT NULL AFTER end;

DROP TABLE IF EXISTS structuralVariantCluster;
DROP TABLE IF EXISTS copyNumberCluster;
CREATE TABLE copyNumberCluster
(   id int NOT NULL AUTO_INCREMENT,
    modified DATETIME NOT NULL,
    sampleId varchar(255) NOT NULL,
    chromosome varchar(255) NOT NULL,
    start int not null,
    end int not null,
    firstVariant int,
    finalVariant int,
    variantCount int not null,
    firstRatio int,
    finalRatio int,
    ratioCount int not null,
    firstBaf int,
    finalBaf int,
    bafCount int not null,
    type varchar(255) NOT NULL,
    PRIMARY KEY (id),
    INDEX(sampleId)
);

ALTER TABLE copyNumberRegion ADD COLUMN gcContent DOUBLE PRECISION not null AFTER observedTumorRatioCount;


ALTER TABLE copyNumber ADD COLUMN copyNumberMethod varchar(255) NOT NULL AFTER copyNumber;
ALTER TABLE copyNumber DROP COLUMN inferred;