ALTER TABLE biopsy
    ADD COLUMN biopsyType varchar(255) AFTER biopsyEvaluable;

ALTER TABLE baseline
    DROP COLUMN cancerType;

ALTER TABLE baseline
    ADD COLUMN primaryTumorLocation varchar(255) AFTER birthYear;
