ALTER TABLE somaticVariant
    ADD kataegis varchar(20) NOT NULL AFTER recovered,
    DROP clonality,
    ADD subclonalLikelihood DOUBLE PRECISION NOT NULL AFTER trinucleotideContext;

ALTER TABLE driverCatalog
    ADD chromosome varchar(255) NOT NULL AFTER sampleId,
    ADD chromosomeBand varchar(255) NOT NULL AFTER chromosome,
    ADD minCopyNumber DOUBLE PRECISION NOT NULL AFTER biallelic,
    ADD maxCopyNumber DOUBLE PRECISION NOT NULL AFTER minCopyNumber;