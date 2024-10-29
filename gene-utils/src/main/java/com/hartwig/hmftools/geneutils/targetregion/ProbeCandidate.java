package com.hartwig.hmftools.geneutils.targetregion;

import com.hartwig.hmftools.common.region.ChrBaseRegion;

public class ProbeCandidate
{
    private final ChrBaseRegion chrBaseRegion;
    private final String sequence;
    private final double gcContent;
    private double sumBlastnBitScore = Double.NaN;

    private String filterReason = "";

    public ProbeCandidate(final ChrBaseRegion chrBaseRegion, final String sequence, final double gcContent)
    {
        this.chrBaseRegion = chrBaseRegion;
        this.sequence = sequence;
        this.gcContent = gcContent;
    }

    public ChrBaseRegion getChrBaseRegion()
    {
        return chrBaseRegion;
    }

    public int getStart()
    {
        return chrBaseRegion.start();
    }

    public int getEnd()
    {
        return chrBaseRegion.end();
    }

    public String getSequence()
    {
        return sequence;
    }

    public double getGcContent()
    {
        return gcContent;
    }

    public double getSumBlastnBitScore()
    {
        return sumBlastnBitScore;
    }

    public String getFilterReason()
    {
        return filterReason;
    }

    public void setSumBlastnBitScore(final double sumBlastnBitScore)
    {
        this.sumBlastnBitScore = sumBlastnBitScore;
    }

    public void setFilterReason(final String filterReason)
    {
        this.filterReason = filterReason;
    }

    public boolean passFilter()
    {
        return filterReason.isEmpty();
    }
}
