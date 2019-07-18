ALTER TABLE purity
  ADD COLUMN msIndelsPerMb DOUBLE PRECISION not null AFTER maxDiploidProportion,
  ADD COLUMN msStatus varchar(10) not null AFTER msIndelsPerMb;

UPDATE purity set msStatus = "UNKNOWN";