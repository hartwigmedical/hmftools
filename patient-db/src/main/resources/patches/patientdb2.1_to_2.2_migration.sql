ALTER TABLE sample
    ADD COLUMN dnaNanograms int AFTER samplingDate;

CREATE TABLE metric
(   sampleId varchar(255) NOT NULL,
    refMeanCoverage DOUBLE PRECISION NOT NULL,
    refCoverage10xPercentage DOUBLE PRECISION NOT NULL,
    refCoverage20xPercentage DOUBLE PRECISION NOT NULL,
    tumorMeanCoverage DOUBLE PRECISION NOT NULL,
    tumorCoverage30xPercentage DOUBLE PRECISION NOT NULL,
    tumorCoverage60xPercentage DOUBLE PRECISION NOT NULL,
    PRIMARY KEY (sampleId)
);