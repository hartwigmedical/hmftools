ALTER TABLE svLink
    CHANGE ploidy junctionCopyNumber DOUBLE PRECISION,
    CHANGE ploidyUncertainty junctionCopyNumberUncertainty DOUBLE PRECISION,
    ADD ecDna BOOLEAN;

ALTER TABLE svAnnotation
    CHANGE ploidyMin junctionCopyNumberMin DOUBLE PRECISION,
    CHANGE ploidyMax junctionCopyNumberMax DOUBLE PRECISION;

