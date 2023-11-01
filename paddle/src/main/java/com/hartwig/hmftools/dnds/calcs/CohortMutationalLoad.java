package com.hartwig.hmftools.dnds.calcs;

import java.util.Map;

import com.hartwig.hmftools.dnds.SampleMutationalLoad;

public class CohortMutationalLoad
{
    public final int SampleCount;
    public final int SnvCount;
    public final int IndelCount;

    public CohortMutationalLoad(final int sampleCount, final int snvCount, final int indelCount)
    {
        SampleCount = sampleCount;
        SnvCount = snvCount;
        IndelCount = indelCount;
    }

    public static CohortMutationalLoad fromSampleMap(final Map<String,SampleMutationalLoad> sampleMutationalLoadMap)
    {
        int snvTotal = 0;
        int indelTotal = 0;

        for(SampleMutationalLoad sampleMutationalLoad : sampleMutationalLoadMap.values())
        {
            snvTotal += sampleMutationalLoad.snvTotal();
            indelTotal += sampleMutationalLoad.indelTotal();
        }

        return new CohortMutationalLoad(sampleMutationalLoadMap.size(), snvTotal, indelTotal);
    }
}
