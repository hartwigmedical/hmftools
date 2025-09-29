package com.hartwig.hmftools.cobalt.calculations;

import com.hartwig.hmftools.cobalt.count.DepthReading;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;

public class BamRatio
{
    public final Chromosome mChromosome;
    public final int Position;
    private final double mReadDepth;
    private double Ratio;
    private double DiploidAdjustedRatio = -1.0;
    private final double GcContent;
    private boolean Included;

    public BamRatio(Chromosome chromosome, DepthReading readDepth, boolean inTargetRegion)
    {
        mChromosome = chromosome;
        Position = readDepth.StartPosition;
        mReadDepth = Double.isFinite(readDepth.ReadDepth) ? readDepth.ReadDepth : -1.0;
        Ratio = mReadDepth;
        GcContent = Double.isFinite(readDepth.ReadGcContent) ? readDepth.ReadGcContent : -1.0;
        Included = inTargetRegion;
        if (!Included)
        {
            Ratio = -1.0;
        }
    }

    public void normaliseForGc(double medianReadDepthForGcBucket)
    {
        normalise(medianReadDepthForGcBucket);
    }

    public void applyEnrichment(double enrichment)
    {
        normalise(enrichment);
    }

    public void normaliseByMean(double mean)
    {
        normalise(mean);
    }

    public void setDiploidAdjustedRatio(double ratio)
    {
        DiploidAdjustedRatio = ratio;
    }

    public double getDiploidAdjustedRatio()
    {
        return DiploidAdjustedRatio;
    }

    private void normalise(final double factor)
    {
        if (!Included)
        {
            return;
        }
        if(factor <= 0 || Double.isNaN(factor))
        {
            Included = false;
            Ratio = -1.0;
        }
        else
        {
            Ratio = Ratio / factor;
        }
    }

    public double readDepth()
    {
        return Included ? mReadDepth : -1.0;
    }

    public double ratio()
    {
        return Included ? Ratio : -1.0;
    }

    public double gcContent()
    {
        return Included ? GcContent : -1.0;
    }

    @Override
    public String toString()
    {
        return "BamRatio{" +
                "mChromosome=" + mChromosome +
                ", Position=" + Position +
                ", mReadDepth=" + mReadDepth +
                ", Ratio=" + Ratio +
                ", GcContent=" + GcContent +
                ", Included=" + Included +
                '}';
    }
}
