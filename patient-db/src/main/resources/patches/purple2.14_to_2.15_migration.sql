ALTER TABLE purity
    ADD somaticDeviation DOUBLE PRECISION not null AFTER score;


ALTER TABLE purityRange
    ADD somaticDeviation DOUBLE PRECISION not null AFTER score;


ALTER TABLE copyNumberRegion
    DROP COLUMN modelPloidy,
    DROP COLUMN modelBaf,
    DROP COLUMN modelTumorRatio,
    DROP COLUMN cnvDeviation,
    DROP COLUMN bafDeviation,
    ADD COLUMN minorAllelePloidyDeviation DOUBLE PRECISION not null AFTER refNormalisedTumorCopyNumber,
    ADD COLUMN majorAllelePloidyDeviation DOUBLE PRECISION not null AFTER minorAllelePloidyDeviation,
    ADD COLUMN minorAllelePloidy DOUBLE PRECISION not null AFTER gcContent,
    ADD COLUMN majorAllelePloidy DOUBLE PRECISION not null AFTER minorAllelePloidy;