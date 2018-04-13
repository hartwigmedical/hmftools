ALTER TABLE baseline
    ADD COLUMN preTreatmentsType varchar(255) AFTER preTreatments;

ALTER TABLE clinicalFindings
    DROP COLUMN ecrfItem;