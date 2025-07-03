package com.hartwig.hmftools.geneutils.paneldesign;

import static java.lang.String.format;

import static com.hartwig.hmftools.geneutils.paneldesign.BlastnResult.INVALID_SCORE;

import com.hartwig.hmftools.common.genome.gc.GcCalcs;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.region.ChrBaseRegion;

public class ProbeCandidate
{
    private final ChrBaseRegion mChrBaseRegion;
    private final String mSequence;
    private final double mGcContent;

    private double mSumBlastnBitScore = Double.NaN;
    private String mFilterReason = "";

    public ProbeCandidate(final ChrBaseRegion chrBaseRegion, final String sequence, final double gcContent)
    {
        mChrBaseRegion = chrBaseRegion;
        mSequence = sequence;
        mGcContent = gcContent;
        mSumBlastnBitScore = INVALID_SCORE;
    }

    public ChrBaseRegion region()
    {
        return mChrBaseRegion;
    }

    public int getStart()
    {
        return mChrBaseRegion.start();
    }
    public int getEnd()
    {
        return mChrBaseRegion.end();
    }

    public String getSequence()
    {
        return mSequence;
    }

    public double getGcContent()
    {
        return mGcContent;
    }

    public double getSumBlastnBitScore()
    {
        return mSumBlastnBitScore;
    }

    public void setSumBlastnBitScore(final double sumBlastnBitScore)
    {
        mSumBlastnBitScore = sumBlastnBitScore;
    }

    public void setFilterReason(final String filterReason)
    {
        this.mFilterReason = filterReason;
    }

    public boolean passFilter()
    {
        return mFilterReason.isEmpty();
    }

    public String toString()
    {
        return format("region(%s) gc(%.3f) score(%.3f) sequence(%s)", mChrBaseRegion, mGcContent, mSumBlastnBitScore, mSequence);
    }

    public static ProbeCandidate createProbeCandidate(final ChrBaseRegion chrBaseRegion, final RefGenomeInterface refGenome)
    {
        // get the sequence from the ref genome
        String sequence = refGenome.getBaseString(chrBaseRegion.chromosome(), chrBaseRegion.start(), chrBaseRegion.end());
        double gcContent = GcCalcs.calcGcPercent(sequence);

        return new ProbeCandidate(chrBaseRegion, sequence, gcContent);
    }

}
