package com.hartwig.hmftools.cobalt.calculations;

import java.util.Objects;

import com.hartwig.hmftools.cobalt.count.DepthReading;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.utils.Doubles;

public class BamRatio
{
    public final HumanChromosome Chromosome;
    public final int Position;
    public final double ReadDepth;
    public final double GcContent;

    private double mRatio;
    private double mDiploidAdjustedRatio;
    private boolean mIncluded;

    public BamRatio(HumanChromosome chromosome, DepthReading readDepth, boolean inTargetRegion)
    {
        this(chromosome, readDepth.StartPosition, readDepth.ReadDepth, readDepth.ReadGcContent, inTargetRegion);
    }

    public BamRatio(HumanChromosome chromosome, int position, double readDepth, double gcContent)
    {
        this(chromosome, position, readDepth, readDepth, gcContent);
    }

    public BamRatio(HumanChromosome chromosome, int position, double readDepth, double ratio, double gcContent)
    {
        Chromosome = chromosome;
        Position = position;
        ReadDepth = readDepth;
        GcContent = gcContent;

        mRatio = ratio;
        mIncluded = true;
        mDiploidAdjustedRatio = -1.0;
    }

    public BamRatio(HumanChromosome chromosome, int position, double readDepth, double gcContent, boolean included)
    {
        Chromosome = chromosome;
        Position = position;
        ReadDepth = Double.isFinite(readDepth) ? readDepth : -1.0;
        GcContent = Double.isFinite(gcContent) ? gcContent : -1.0;

        mRatio = ReadDepth;
        mDiploidAdjustedRatio = -1.0;
        mIncluded = included;

        if(!mIncluded)
        {
            mRatio = -1.0;
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
        if(mDiploidAdjustedRatio == 0.0)
        {
            return;
        }
        if(factor <= 0 || Double.isNaN(factor) || mDiploidAdjustedRatio < 0)
        {
            mDiploidAdjustedRatio = -1.0;
        }
        else
        {
            mDiploidAdjustedRatio = mDiploidAdjustedRatio / factor;
        }
    }

    public void setDiploidAdjustedRatio(double ratio)
    {
        if(mRatio == 0.0)
        {
            mDiploidAdjustedRatio = 0.0;
        }
        else
        {
            mDiploidAdjustedRatio = ratio;
        }
    }

    public double getDiploidAdjustedRatio()
    {
        return mDiploidAdjustedRatio;
    }

    public void overrideRatio(double ratio)
    {
        mRatio = ratio;
        if(mRatio > 0)
        {
            mIncluded = true;
        }
    }

    private void normalise(final double factor)
    {
        if(Double.isNaN(factor))
        {
            mIncluded = false;
            mRatio = -1.0;
            return;
        }
        if(Doubles.isZero(mRatio))
        {
            return;
        }
        if(factor <= 0 || !mIncluded || mRatio < 0)
        {
            mIncluded = false;
            mRatio = -1.0;
        }
        else
        {
            mRatio /= factor;
        }
    }

    public int position()
    {
        return Position;
    }

    public double readDepth()
    {
        return ReadDepth;
    }

    public double ratio()
    {
        return mIncluded ? mRatio : -1.0;
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
        return Position == bamRatio.Position && Objects.equals(Chromosome, bamRatio.Chromosome);
    }

    @Override
    public int hashCode()
    {
        return Objects.hash(Chromosome, Position);
    }

    @Override
    public String toString()
    {
        return "BamRatio{" +
                "mChromosome=" + Chromosome +
                ", Position=" + Position +
                ", mReadDepth=" + ReadDepth +
                ", Ratio=" + mRatio +
                ", GcContent=" + GcContent +
                ", Included=" + mIncluded +
                '}';
    }
}
