package com.hartwig.hmftools.cobalt.calculations;

import java.util.Objects;

import com.google.common.base.Preconditions;
import com.hartwig.hmftools.cobalt.count.DepthReading;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;

class CobaltWindow
{
    public final Chromosome mChromosome;
    public final int Position;
    private final Double ReadDepth;
    public final Double GcContent;
    public final GCPail GcBucket;
    public final boolean IsInExcludedRegion;
    public final boolean IsInTargetRegion;

    CobaltWindow(final Chromosome chromosome, final DepthReading depth, boolean isInExcludedRegion, boolean isInTargetRegion)
    {
        mChromosome = chromosome;
        Position = depth.StartPosition;
        ReadDepth = depth.ReadDepth;
        GcContent = depth.ReadGcContent;
        this.GcBucket = null;
        IsInExcludedRegion = isInExcludedRegion;
        IsInTargetRegion = isInTargetRegion;
    }

    CobaltWindow(Chromosome chromosome, int position, double readDepth, double gcContent, GCPail GcBucket, boolean isInTargetRegion)
    {
        mChromosome = chromosome;
        Position = position;
        ReadDepth = readDepth;
        GcContent = gcContent;
        this.GcBucket = GcBucket;
        IsInExcludedRegion = false;
        IsInTargetRegion = isInTargetRegion;
    }

    boolean include()
    {
        return !IsInExcludedRegion && IsInTargetRegion;
    }

    BamRatio toBamRatio()
    {
        return new BamRatio(mChromosome, Position, ReadDepth, GcContent, include());
    }

    CobaltWindow correctedByReferenceValue(WindowStatuses genomeData)
    {
        Preconditions.checkArgument(this.GcBucket == null);
        if(ReadDepth >= 1.0 || ReadDepth < 0)
        {
            return this;
        }
        if(!include())
        {
            return this;
        }
        double referenceValue;
        try
        {
            referenceValue = genomeData.referenceGcValueForWindow(mChromosome, Position);
        }
        catch(Exception e)
        {
            return this;
        }
        return new CobaltWindow(mChromosome, Position, ReadDepth, referenceValue, null, IsInTargetRegion);
    }

    CobaltWindow bucketed(final GCPail bucket)
    {
        Preconditions.checkArgument(this.GcBucket == null);
        Preconditions.checkArgument(GCPail.bucketIndex(GcContent) == bucket.mGC);
        if(IsInExcludedRegion)
        {
            return this;
        }
        else
        {
            if(this.mChromosome.isAutosome() && IsInTargetRegion && ReadDepth > 0)
            {
                bucket.addReading(ReadDepth);
            }
            return new CobaltWindow(mChromosome, Position, ReadDepth, GcContent, bucket, IsInTargetRegion);
        }
    }

    @Override
    public boolean equals(final Object o)
    {
        if(o == null || getClass() != o.getClass())
        {
            return false;
        }
        final CobaltWindow that = (CobaltWindow) o;
        return Position == that.Position && mChromosome == that.mChromosome;
    }

    @Override
    public int hashCode()
    {
        return Objects.hash(mChromosome, Position);
    }

    @Override
    public String toString()
    {
        return "CobaltWindow{" +
                "Chromosome=" + mChromosome +
                ", Position=" + Position +
                ", ReadDepth=" + ReadDepth +
                ", GcContent=" + GcContent +
                ", IsInExcludedRegion=" + IsInExcludedRegion +
                ", IsInTargetRegion=" + IsInTargetRegion +
                '}';
    }
}
