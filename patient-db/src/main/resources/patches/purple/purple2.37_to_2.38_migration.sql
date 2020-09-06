ALTER TABLE purity
    ADD COLUMN tmbPerMb DOUBLE PRECISION not null AFTER msStatus,
    ADD COLUMN tmbStatus varchar(10) not null AFTER tmbPerMb,
    ADD COLUMN tml DOUBLE PRECISION not null AFTER tmbStatus,
    ADD COLUMN tmlStatus varchar(10) not null AFTER tml;


UPDATE purity SET tmbPerMb = 0, tml = 0, tmbStatus = "UNKNOWN", tmlStatus = "UNKNOWN";
