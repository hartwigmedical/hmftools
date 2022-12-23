package com.hartwig.hmftools.bamtools.markdups;

import static java.lang.Math.min;

public class DuplicateStats
{
    public int MarkedCount;
    public int IdentifiedCount;
    public int MisMatches;
    public int PartiallyMarked;

    public int[] FrequencyCounts;

    public DuplicateStats()
    {
        MarkedCount = 0;
        IdentifiedCount = 0;
        MisMatches = 0;
        PartiallyMarked = 0;
        FrequencyCounts = new int[100];
    }

    public void addFrequency(int frequency)
    {
        ++FrequencyCounts[min(frequency, FrequencyCounts.length - 1)];
    }
}
