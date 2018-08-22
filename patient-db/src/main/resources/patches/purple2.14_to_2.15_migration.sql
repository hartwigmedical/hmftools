ALTER TABLE purity
    ADD somaticDeviation DOUBLE PRECISION not null AFTER score;


ALTER TABLE purityRange
    ADD somaticDeviation DOUBLE PRECISION not null AFTER score;