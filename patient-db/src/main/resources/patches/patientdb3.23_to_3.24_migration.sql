ALTER TABLE somaticVariant
    ADD kataegis varchar(20) NOT NULL AFTER recovered,
    DROP clonality,
    ADD subclonalLikelihood DOUBLE PRECISION NOT NULL AFTER trinucleotideContext;

ALTER TABLE driverCatalog
    ADD minCopyNumber DOUBLE PRECISION NOT NULL AFTER biallelic;