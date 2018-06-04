ALTER TABLE drug
    ADD COLUMN mechanism varchar(255) AFTER type;

ALTER TABLE preTreatmentDrug
    ADD COLUMN mechanism varchar(255) AFTER type;

ALTER TABLE sample
    ADD COLUMN limsPrimaryTumor varchar(255) after dnaNanograms;