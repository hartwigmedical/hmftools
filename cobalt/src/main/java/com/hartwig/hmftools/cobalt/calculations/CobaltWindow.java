package com.hartwig.hmftools.cobalt.calculations;

import java.util.Objects;

import com.google.common.base.Preconditions;
import com.hartwig.hmftools.cobalt.count.DepthReading;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;

public class CobaltWindow
{
    public final Chromosome Chromosome;
    public final int Position;
    public final DepthReading mDepthReading;
    public final GCPail GcBucket;
    public final boolean IsInExcludedRegion;
    public final boolean IsInTargetRegion;

    public CobaltWindow(final Chromosome chromosome, final DepthReading depth, boolean isInExcludedRegion, boolean isInTargetRegion)
    {
        Chromosome = chromosome;
        Position = depth.StartPosition;
        this.mDepthReading = depth;
        this.GcBucket = null;
        IsInExcludedRegion = isInExcludedRegion;
        IsInTargetRegion = isInTargetRegion;
    }

    public CobaltWindow(final Chromosome chromosome, DepthReading depthReading, final GCPail GcBucket, boolean isInTargetRegion)
    {
        Chromosome = chromosome;
        Position = depthReading.StartPosition;
        this.mDepthReading = depthReading;
        this.GcBucket = GcBucket;
        IsInExcludedRegion = false;
        IsInTargetRegion = isInTargetRegion;
    }

    public CobaltWindow bucketed(final GCPail bucket)
    {
        Preconditions.checkArgument(this.GcBucket == null);
        Preconditions.checkArgument(GCPail.bucketIndex(mDepthReading.ReadGcContent) == bucket.mGC);
        if(IsInExcludedRegion)
        {
            return this;
        }
        else
        {
            if(this.Chromosome.isAutosome() && IsInTargetRegion && mDepthReading.ReadDepth > 0)
            {
                bucket.addReading(mDepthReading.ReadDepth);
            }
            return new CobaltWindow(Chromosome, mDepthReading, bucket, IsInTargetRegion);
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
        return Position == that.Position && Chromosome == that.Chromosome;
    }

    @Override
    public int hashCode()
    {
        return Objects.hash(Chromosome, Position);
    }

    @Override
    public String toString()
    {
        return "CobaltWindow{" +
                "Chromosome=" + Chromosome +
                ", Position=" + Position +
                ", ReadDepth=" + mDepthReading +
                '}';
    }
}
