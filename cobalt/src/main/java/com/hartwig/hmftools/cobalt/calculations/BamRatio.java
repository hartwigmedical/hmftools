package com.hartwig.hmftools.cobalt.calculations;

import java.util.Objects;

import com.hartwig.hmftools.cobalt.count.DepthReading;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.utils.Doubles;

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
        this(chromosome, readDepth.StartPosition, readDepth.ReadDepth, readDepth.ReadGcContent, inTargetRegion);
    }

    public BamRatio(Chromosome chromosome, int position, double readDepth, double gcContent)
    {
        this(chromosome, position, readDepth, readDepth, gcContent);
    }

    public BamRatio(Chromosome chromosome, int position, double readDepth, double ratio, double gcContent)
    {
        mChromosome = chromosome;
        Position = position;
        mReadDepth = readDepth;
        Ratio = ratio;
        GcContent = gcContent;
        Included = true;
    }

    public BamRatio(Chromosome chromosome, int position, double readDepth, double gcContent, boolean included)
    {
        mChromosome = chromosome;
        Position = position;
        mReadDepth = Double.isFinite(readDepth) ? readDepth : -1.0;
        Ratio = mReadDepth;
        GcContent = Double.isFinite(gcContent) ? gcContent : -1.0;
        Included = included;
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

    public void normaliseDiploidAdjustedRatio(double factor)
    {
        // ratio 0 => {-1 if not relevant, 0 if relevant}
        if (DiploidAdjustedRatio == 0.0)
        {
            return;
        }
        if (factor <= 0 || Double.isNaN(factor) || DiploidAdjustedRatio < 0)
        {
            DiploidAdjustedRatio = -1.0;
        }
        else
        {
            DiploidAdjustedRatio = DiploidAdjustedRatio / factor;
        }
    }

    public void setDiploidAdjustedRatio(double ratio)
    {
        if (Ratio == 0.0)
        {
            DiploidAdjustedRatio = 0.0;
        }
        else
        {
            DiploidAdjustedRatio = ratio;
        }
    }

    public double getDiploidAdjustedRatio()
    {
        return DiploidAdjustedRatio;
    }

    public void overrideRatio(double ratio)
    {
        Ratio = ratio;
        if (Ratio > 0)
        {
            Included = true;
        }
    }

    private void normalise(final double factor)
    {
        if(Double.isNaN(factor))
        {
            Included = false;
            Ratio = -1.0;
        }
        if (Doubles.isZero(Ratio))
        {
            return;
        }
        if (!Included | Ratio < 0)
        {
            Ratio = -1.0;
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

    public int position()
    {
        return Position;
    }

    public double readDepth()
    {
        return mReadDepth;
    }

    public double ratio()
    {
        return Included ? Ratio : -1.0;
    }

    public double gcContent()
    {
        return GcContent;
    }

    @Override
    public boolean equals(final Object o)
    {
        if(o == null || getClass() != o.getClass())
        {
            return false;
        }
        final BamRatio bamRatio = (BamRatio) o;
        return Position == bamRatio.Position && Objects.equals(mChromosome, bamRatio.mChromosome);
    }

    @Override
    public int hashCode()
    {
        return Objects.hash(mChromosome, Position);
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
