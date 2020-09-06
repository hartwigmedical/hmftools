ALTER TABLE geneCopyNumber
    DROP COLUMN nonsenseBiallelicVariants,
    DROP COLUMN nonsenseNonBiallelicVariants,
    DROP COLUMN nonsenseNonBiallelicPloidy,
    DROP COLUMN spliceBiallelicVariants,
    DROP COLUMN spliceNonBiallelicVariants,
    DROP COLUMN spliceNonBiallelicPloidy,
    DROP COLUMN missenseBiallelicVariants,
    DROP COLUMN missenseNonBiallelicVariants,
    DROP COLUMN missenseNonBiallelicPloidy;

ALTER TABLE copyNumber
    ADD COLUMN minorAllelePloidy DOUBLE PRECISION not null AFTER copyNumber,
    ADD COLUMN majorAllelePloidy DOUBLE PRECISION not null AFTER minorAllelePloidy,
    CHANGE actualBAF baf DOUBLE PRECISION not null;

UPDATE copyNumber
SET minorAllelePloidy = round(greatest(0, (1- baf) * copyNumber), 3), majorAllelePloidy = round(copyNumber - minorAllelePloidy, 3);


ALTER TABLE copyNumberGermline
    ADD COLUMN minorAllelePloidy DOUBLE PRECISION not null AFTER copyNumber,
    ADD COLUMN majorAllelePloidy DOUBLE PRECISION not null AFTER minorAllelePloidy,
    CHANGE actualBAF baf DOUBLE PRECISION not null;

UPDATE copyNumberGermline
SET minorAllelePloidy = round(greatest(0, (1- baf) * copyNumber), 3), majorAllelePloidy = round(copyNumber - minorAllelePloidy, 3);


ALTER TABLE purityRange
    CHANGE somaticDeviation somaticPenalty DOUBLE PRECISION not null;

ALTER TABLE purity
    CHANGE somaticDeviation somaticPenalty DOUBLE PRECISION not null;

ALTER TABLE copyNumberRegion
    CHANGE totalDeviation deviationPenalty DOUBLE PRECISION not null,
    CHANGE ploidyPenalty eventPenalty DOUBLE PRECISION not null;