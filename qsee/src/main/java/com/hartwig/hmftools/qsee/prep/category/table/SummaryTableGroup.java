package com.hartwig.hmftools.qsee.prep.category.table;

enum SummaryTableGroup
{
    TMB("Mutational burden"),
    CONTAMINATION("Contamination"),
    COPY_NUMBER("Copy number"),
    MAPPING("Mapping");

    private final String mHumanReadableName;

    SummaryTableGroup(String humanReadableName)
    {
        mHumanReadableName = humanReadableName;
    }

    public String humanReadableName() { return mHumanReadableName; }
}
