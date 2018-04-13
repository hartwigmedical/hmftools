ALTER TABLE baseline
    ADD COLUMN preTreatmentsType varchar(510) AFTER preTreatments;

ALTER TABLE clinicalFindings
    DROP COLUMN ecrfItem;