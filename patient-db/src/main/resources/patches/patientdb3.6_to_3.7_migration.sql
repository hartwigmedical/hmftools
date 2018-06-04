ALTER TABLE drug
    ADD COLUMN mechanism varchar(255) AFTER type;

ALTER TABLE preTreatmentDrug
    ADD COLUMN mechanism varchar(255) AFTER type;

ALTER TABLE sample
    ADD COLUMN limsPrimaryTumor varchar(255) after dnaNanograms;

ALTER TABLE baseline
    ADD COLUMN preTreatmentsMechanism varchar(510) after preTreatmentsType;

ALTER TABLE treatment
    ADD COLUMN mechanism varchar(255) after type;

ALTER TABLE baseline DROP COLUMN preTreatments;
ALTER TABLE baseline
    ADD COLUMN preTreatments varchar(800) after hasRadiotherapyPreTreatment;

ALTER TABLE preTreatmentDrug DROP COLUMN name;
ALTER TABLE preTreatmentDrug
    ADD COLUMN name varchar(800) after endDate;

ALTER TABLE treatment DROP COLUMN name;
ALTER TABLE treatment
    ADD COLUMN name varchar(800) after endDate;

ALTER TABLE drug DROP COLUMN name;
ALTER TABLE drug
    ADD COLUMN name varchar(800)after endDate;