ALTER TABLE somaticVariant
    ADD COLUMN minorAllelePloidy DOUBLE PRECISION NOT NULL AFTER germlineStatus;


ALTER TABLE copyNumberRegion
    ADD COLUMN ploidyPenalty DOUBLE PRECISION not null AFTER totalDeviation;