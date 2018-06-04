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