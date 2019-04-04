ALTER TABLE structuralVariant
    ADD COLUMN startAnchoringSupportDistance int after insertSequenceRepeatCoverage,
    ADD COLUMN endAnchoringSupportDistance int after startAnchoringSupportDistance,
    CHANGE COLUMN startTumourVariantFragmentCount startTumorVariantFragmentCount int,
    CHANGE COLUMN startTumourReferenceFragmentCount startTumorReferenceFragmentCount int,
    CHANGE COLUMN endTumourVariantFragmentCount endTumorVariantFragmentCount int,
    CHANGE COLUMN endTumourReferenceFragmentCount endTumorReferenceFragmentCount int;