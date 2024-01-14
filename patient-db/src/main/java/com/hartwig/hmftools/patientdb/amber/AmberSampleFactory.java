package com.hartwig.hmftools.patientdb.amber;

import java.util.List;

import com.hartwig.hmftools.common.amber.NormalHeterozygousCheck;

public class AmberSampleFactory
{
    private final int mMinReadDepth;
    private final NormalHeterozygousCheck mHeterozygousFilter;

    public AmberSampleFactory(final int minReadDepth, final double minHetAFPercentage, final double maxHetAFPercentage)
    {
        mMinReadDepth = minReadDepth;
        mHeterozygousFilter = new NormalHeterozygousCheck(minHetAFPercentage, maxHetAFPercentage);
    }

    public AmberSample createSampleData(final String sample, final List<SiteEvidence> baseDepths)
    {
        byte[] entries = new byte[baseDepths.size()];
        for(int i = 0; i < baseDepths.size(); i++)
        {
            entries[i] = asByte(baseDepths.get(i));
        }

        return ImmutableAmberSample.builder().sampleId(sample).entries(entries).build();
    }

    public byte asByte(final SiteEvidence siteEvidence)
    {
        if(siteEvidence.ReadDepth < mMinReadDepth)
        {
            return AmberSample.DO_NOT_MATCH;
        }

        if(siteEvidence.RefSupport == siteEvidence.ReadDepth)
        {
            return (byte) 1;
        }

        if(mHeterozygousFilter.test(siteEvidence.ReadDepth, siteEvidence.RefSupport, siteEvidence.AltSupport, 0))
        {
            return (byte) 2;
        }

        if(siteEvidence.AltSupport == siteEvidence.ReadDepth)
        {
            return (byte) 3;
        }

        return AmberSample.DO_NOT_MATCH;
    }
}
