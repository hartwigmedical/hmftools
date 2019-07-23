ALTER TABLE svFusion
    CHANGE skippedExons skippedExonsUp INT,
    ADD skippedExonsDown INT;

ALTER TABLE svLink
    ADD ploidyUncertainty DOUBLE PRECISION,
    CHANGE chainIndex chainIndex VARCHAR(50);
