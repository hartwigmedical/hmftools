ALTER TABLE geneCopyNumber
    ADD COLUMN minRegions int not null AFTER germlineHetRegions,
    ADD COLUMN minRegionStart int not null AFTER minRegions,
    ADD COLUMN minRegionEnd int not null AFTER minRegionStart,
    ADD COLUMN minRegionStartSupport varchar(255) NOT NULL AFTER minRegionEnd,
    ADD COLUMN minRegionEndSupport varchar(255) NOT NULL AFTER minRegionStartSupport,
    ADD COLUMN minRegionMethod varchar(255) NOT NULL AFTER minRegionEndSupport;