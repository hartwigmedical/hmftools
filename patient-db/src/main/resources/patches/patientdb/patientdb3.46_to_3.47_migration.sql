ALTER TABLE chord
    DROP COLUMN nothing,
    DROP COLUMN predictedResponse,
    ADD COLUMN hrStatus varchar(255) NOT NULL,
    ADD COLUMN hrdType varchar(255) NOT NULL,
    ADD COLUMN remarksHrStatus varchar(255) NOT NULL,
    ADD COLUMN remarksHrdType varchar(255) NOT NULL;

UPDATE chord SET hrStatus = "N/A";
UPDATE chord SET hrdType = "N/A";
UPDATE chord SET remarksHrStatus = "N/A";
UPDATE chord SET remarksHrdType = "N/A";

