package com.hartwig.hmftools.cobalt.calculations;

import com.hartwig.hmftools.cobalt.count.ReadDepth;
import com.hartwig.hmftools.common.cobalt.CobaltRatio;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;

public class BamRatio
{
    public final Chromosome mChromosome;
    public final int Position;
    private double mReadDepth;
    private double Ratio;
    private final double GcContent;
    private boolean Included;

    public BamRatio(Chromosome chromosome, ReadDepth readDepth, boolean inTargetRegion)
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

    public CobaltRatio toTumorRatio(RefGenomeVersion version)
    {
        double r = Ratio;
        double g = GcContent;
        double rd = mReadDepth;
        return new CobaltRatio(version.versionedChromosome(mChromosome), Position, -1.0, -1.0, -1.0, -1.0, rd, r, g);
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
