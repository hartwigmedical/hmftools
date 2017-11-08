ALTER TABLE copyNumber CHANGE structuralVariantSupport segmentStartSupport varchar(255) NOT NULL;
ALTER TABLE copyNumber ADD COLUMN segmentEndSupport varchar(255) NOT NULL AFTER segmentStartSupport;
ALTER TABLE copyNumberRegion CHANGE structuralVariantSupport segmentStartSupport varchar(255) NOT NULL;


ALTER TABLE copyNumberRegion ADD COLUMN actualTumorBaf DOUBLE PRECISION not null AFTER modelTumorRatio;