package com.hartwig.hmftools.svannotation.analysis;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.structural.EnrichedStructuralVariant;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantData;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantType;

import java.util.List;

public class SvClusterData
{
    private final String mId; // sourced from either VCF or DB

    private final String mStartChromosome;
    private final String mEndChromosome;
    private final long mStartPosition;
    private final long mEndPosition;

    // other fields as required
    private final byte mStartOrientation;
    private final byte mEndOrientation;
    private StructuralVariantType mType;

    // clustering info
    private List<SvClusterData> mSubSVs; // other SVs wholy contained within this SV
    private boolean mIsSubSV;
    private SvCluster mStartCluster;
    private SvCluster mEndCluster;

    public SvClusterData(
            final String id,
            final String startChromosome,
            final String endChromosome,
            final long startPosition,
            final long endPosition,
            final byte startOrientation,
            final byte endOrientation,
            final StructuralVariantType type)
    {
        mId = id;
        mStartChromosome = startChromosome;
        mEndChromosome = endChromosome;
        mStartPosition = startPosition;
        mEndPosition = endPosition;

        mStartOrientation = startOrientation;
        mEndOrientation = endOrientation;
        mType = type;

        mSubSVs = Lists.newArrayList();
        mIsSubSV = false;
        mStartCluster = null;
        mEndCluster = null;
    }

    public static SvClusterData from(final EnrichedStructuralVariant enrichedSV)
    {
        return new SvClusterData(
                enrichedSV.id(),
                enrichedSV.chromosome(true),
                enrichedSV.chromosome(false),
                enrichedSV.position(true),
                enrichedSV.position(false),
                enrichedSV.orientation(true),
                enrichedSV.orientation(false),
                enrichedSV.type());
    }

    public static SvClusterData from(final StructuralVariantData svRecord)
    {
        return new SvClusterData(
                svRecord.id(),
                svRecord.startChromosome(),
                svRecord.endChromosome(),
                svRecord.startPosition(),
                svRecord.endPosition(),
                svRecord.startOrientation(),
                svRecord.endOrientation(),
                svRecord.type());
    }

    public final String id() { return mId; }

    public final String posId()
    {
        return String.format("id(%s) position(%s:%d -> %s:%d)",
                id(), chromosome(true), position(true), chromosome(false), position(false));
    }

    public boolean equals(final SvClusterData other) { return id().equals(other.id()); }

    public final String chromosome(boolean isStart) { return isStart ? mStartChromosome : mEndChromosome; }
    public final long position(boolean isStart) { return isStart ? mStartPosition : mEndPosition; }
    public final byte orientation(boolean isStart) { return isStart ? mStartOrientation: mEndOrientation; }
    public  StructuralVariantType type() { return mType; }

    public final boolean isLocal()
    {
        // means that both ends are within the same chromosomal arm
        return chromosome(true).equals(chromosome(false));
    }

    public final long getSpan()
    {
        if(!chromosome(true).equals(chromosome(false)))
            return -1;

        return position(false) - position(true);
    }

    public boolean addSubSV(SvClusterData subVariant)
    {
        for(final SvClusterData existing : mSubSVs)
        {
            if(existing.id() == subVariant.id())
                return false;
        }

        mSubSVs.add(subVariant);
        subVariant.setIsSubSV(true);
        return true;
    }

    public List<SvClusterData> getSubSVs() { return mSubSVs; }
    public boolean hasSubSVs() { return !mSubSVs.isEmpty(); }

    public void setIsSubSV(final boolean isSubSV) { mIsSubSV = isSubSV; }
    public final boolean isSubSV() { return mIsSubSV; }

    public void setStartCluster(SvCluster cluster) { mStartCluster = cluster; }
    public void setEndCluster(SvCluster cluster) { mEndCluster = cluster; }
    public final SvCluster getStartCluster() { return mStartCluster; }
    public final SvCluster getEndCluster() { return mEndCluster; }

    public final boolean areClustersSet() { return mStartCluster != null && mEndCluster != null; }
}
