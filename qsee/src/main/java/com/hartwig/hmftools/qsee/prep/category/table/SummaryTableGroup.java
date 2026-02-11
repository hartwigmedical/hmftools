package com.hartwig.hmftools.qsee.prep.category.table;

enum SummaryTableGroup
{
    GENERAL("General"),
    TMB("Mutation burden"),
    COVERAGE("Coverage stats"),
    READ("Read stats");

    private final String mHumanReadableName;

    SummaryTableGroup(String humanReadableName)
    {
        mHumanReadableName = humanReadableName;
    }

    public String humanReadableName() { return mHumanReadableName; }
}
