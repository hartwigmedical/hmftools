package com.hartwig.hmftools.geneutils.paneldesign;

import static java.lang.Double.NaN;
import static java.lang.String.format;

import java.util.Optional;

import com.hartwig.hmftools.common.genome.gc.GcCalcs;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.region.ChrBaseRegion;

public class ProbeCandidate
{
    private final ChrBaseRegion mChrBaseRegion;
    private final String mSequence;
    private final double mGcContent;

    private Optional<Double> mQualityScore = Optional.empty();
    private String mFilterReason = "";

    public ProbeCandidate(final ChrBaseRegion chrBaseRegion, final String sequence, final double gcContent)
    {
        mChrBaseRegion = chrBaseRegion;
        mSequence = sequence;
        mGcContent = gcContent;
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

    public Optional<Double> getQualityScore()
    {
        return mQualityScore;
    }

    public void setQualityScore(double qualityScore)
    {
        mQualityScore = Optional.of(qualityScore);
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
        return format("region(%s) gc(%.3f) score(%.3f) sequence(%s)", mChrBaseRegion, mGcContent, mQualityScore.orElse(NaN), mSequence);
    }

    public static ProbeCandidate createProbeCandidate(final ChrBaseRegion chrBaseRegion, final RefGenomeInterface refGenome)
    {
        // get the sequence from the ref genome
        String sequence = refGenome.getBaseString(chrBaseRegion.chromosome(), chrBaseRegion.start(), chrBaseRegion.end());
        double gcContent = GcCalcs.calcGcPercent(sequence);

        return new ProbeCandidate(chrBaseRegion, sequence, gcContent);
    }

}
