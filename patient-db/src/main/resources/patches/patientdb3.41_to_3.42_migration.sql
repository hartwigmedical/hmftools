ALTER TABLE somaticVariant
    ADD referenceAlleleReadCount int AFTER kataegis,
    ADD referenceTotalReadCount int AFTER referenceAlleleReadCount,
    ADD rnaAlleleReadCount int AFTER referenceTotalReadCount,
    ADD rnaTotalReadCount int AFTER rnaAlleleReadCount;