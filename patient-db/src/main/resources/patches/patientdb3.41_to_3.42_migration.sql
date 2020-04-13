ALTER TABLE somaticVariant
    ADD referenceAlleleReadCount int AFTER kataegis,
    ADD referenceTotalReadCount int AFTER referenceAlleleReadCount,
    ADD rnaAlleleReadCount int AFTER referenceTotalReadCount,
    ADD rnaTotalReadCount int AFTER rnaAlleleReadCount,
    ADD qual double precision not null AFTER rnaTotalReadCount;

ALTER TABLE treatmentResponse
    MODIFY COLUMN response varchar(500);