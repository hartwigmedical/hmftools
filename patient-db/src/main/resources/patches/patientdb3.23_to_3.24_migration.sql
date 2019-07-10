ALTER TABLE somaticVariant
    ADD kataegis varchar(20) NOT NULL AFTER recovered;

UPDATE somaticVariant SET kataegis = "NONE";