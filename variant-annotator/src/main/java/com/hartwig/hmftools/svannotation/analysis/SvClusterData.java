package com.hartwig.hmftools.svannotation.analysis;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.structural.EnrichedStructuralVariant;
import com.hartwig.hmftools.common.variant.structural.ImmutableStructuralVariantData;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantData;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantType;

import java.util.List;

import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.STRUCTURALVARIANT;

public class SvClusterData
{
    private final String mId; // sourced from either VCF or DB

    // full set of DB fields
    private final StructuralVariantData mSVData;
    private String mStartArm;
    private String mEndArm;

    // clustering info
    private List<SvClusterData> mSubSVs; // other SVs wholy contained within this SV
    private boolean mIsSubSV;
    private SvCluster mStartCluster;
    private SvCluster mEndCluster;

    public SvClusterData(final StructuralVariantData svData)
    {
        mId = svData.id();
        mSVData = svData;
        mStartArm = "";
        mEndArm = "";

        mSubSVs = Lists.newArrayList();
        mIsSubSV = false;
        mStartCluster = null;
        mEndCluster = null;
    }

    public static SvClusterData from(final EnrichedStructuralVariant enrichedSV)
    {
        StructuralVariantData svData =
            ImmutableStructuralVariantData.builder()
                .id(enrichedSV.id())
                .startChromosome(enrichedSV.chromosome(true))
                .endChromosome(enrichedSV.chromosome(false))
                .startPosition(enrichedSV.position(true))
                .endPosition(enrichedSV.position(false))
                .startOrientation(enrichedSV.orientation(true))
                .endOrientation(enrichedSV.orientation(false))
                .startAF(enrichedSV.start().alleleFrequency())
                .adjustedStartAF(enrichedSV.start().adjustedAlleleFrequency())
                .adjustedStartCopyNumber(enrichedSV.start().adjustedCopyNumber())
                .adjustedStartCopyNumberChange(enrichedSV.start().adjustedCopyNumberChange())
                .endAF(enrichedSV.end().alleleFrequency())
                .adjustedEndAF(enrichedSV.end().adjustedAlleleFrequency())
                .adjustedEndCopyNumber(enrichedSV.end().adjustedCopyNumber())
                .adjustedEndCopyNumberChange(enrichedSV.end().adjustedCopyNumberChange())
                .ploidy(enrichedSV.ploidy())
                .type(enrichedSV.type())
                .build();

        return new SvClusterData(svData);
    }

    public final String id() { return mId; }
    public final StructuralVariantData getSvData() { return mSVData; }

    // for convenience
    public final byte orientation(boolean isStart){ return isStart ? mSVData.startOrientation() : mSVData.endOrientation(); }
    public final StructuralVariantType type() { return mSVData.type(); }

    public final String getStartArm() { return mStartArm; }
    public final String getEndArm() { return mEndArm; }
    public void setChromosomalArms(final String start, final String end)
    {
        mStartArm = start;
        mEndArm = end;
    }


    public final String posId()
    {
        return String.format("id(%s) position(%s:%d -> %s:%d)",
                id(), chromosome(true), position(true), chromosome(false), position(false));
    }

    public boolean equals(final SvClusterData other) { return id().equals(other.id()); }

    public final String chromosome(boolean isStart) { return isStart ? mSVData.startChromosome() : mSVData.endChromosome(); }
    public final long position(boolean isStart) { return isStart ? mSVData.startPosition() : mSVData.endPosition(); }

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
