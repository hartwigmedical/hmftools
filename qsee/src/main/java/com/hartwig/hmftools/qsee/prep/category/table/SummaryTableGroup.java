package com.hartwig.hmftools.qsee.prep.category.table;

public enum SummaryTableGroup
{
    TMB("Mutational burden"),
    CONTAMINATION("Contamination"),
    COPY_NUMBER("Copy number"),
    MAPPING("Mapping");

    private final String mDisplayName;

    SummaryTableGroup(String displayName)
    {
        mDisplayName = displayName;
    }

    public String displayName() { return mDisplayName; }
}
