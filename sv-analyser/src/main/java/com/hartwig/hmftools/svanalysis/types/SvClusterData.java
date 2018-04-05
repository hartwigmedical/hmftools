package com.hartwig.hmftools.svanalysis.types;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.structural.EnrichedStructuralVariant;
import com.hartwig.hmftools.common.variant.structural.ImmutableStructuralVariantData;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantData;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantLeg;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantType;

import java.util.List;


public class SvClusterData
{
    private final String mId; // sourced from either VCF or DB

    // full set of DB fields
    private final StructuralVariantData mSVData;
    private String mStartArm;
    private String mEndArm;
    private int mPonCount;
    private int mPonRegionCount; // allowing for a small buffer either side of a PON
    private String mStartFragileSite;
    private String mEndFragileSite;
    private String mStartLineElement;
    private String mEndLineElement;

    private List<StructuralVariantLeg> mUniqueBreakends;

    // clustering info
    private List<SvClusterData> mSubSVs; // other SVs wholy contained within this SV
    private boolean mIsSubSV;
//    private SvCluster mStartCluster;
//    private SvCluster mEndCluster;

    private int mNearestSVLength;
    private String mNearestSVLinkType;
    private int mNearestTILength; // templated insertion if exists
    private int mNearestDBLength; // deletion-bridge link if exists

    public SvClusterData(final StructuralVariantData svData)
    {
        mId = svData.id();
        mSVData = svData;
        mStartArm = "";
        mEndArm = "";
        mPonCount = 0;
        mPonRegionCount = 0;
        mStartFragileSite = "";
        mEndFragileSite = "";
        mStartLineElement = "";
        mEndLineElement = "";

        mNearestSVLength = -1;
        mNearestSVLinkType = "NONE";
        mNearestTILength = -1;
        mNearestDBLength = -1;

        mSubSVs = Lists.newArrayList();
        mIsSubSV = false;
//        mStartCluster = null;
//        mEndCluster = null;
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
    public boolean equals(final SvClusterData other) { return id().equals(other.id()); }

    public final String chromosome(boolean isStart) { return isStart ? mSVData.startChromosome() : mSVData.endChromosome(); }
    public final long position(boolean isStart) { return isStart ? mSVData.startPosition() : mSVData.endPosition(); }
    public final byte orientation(boolean isStart){ return isStart ? mSVData.startOrientation() : mSVData.endOrientation(); }
    public final StructuralVariantType type() { return mSVData.type(); }

    public final String posId()
    {
        return String.format("id(%s) position(%s:%d:%d -> %s:%d:%d)",
                id(), chromosome(true), orientation(true), position(true),
                chromosome(false), orientation(false), position(false));
    }

    public final String arm(boolean isStart) { return isStart ? mStartArm : mEndArm; }
    public final String getStartArm() { return mStartArm; }
    public final String getEndArm() { return mEndArm; }
    public void setChromosomalArms(final String start, final String end)
    {
        mStartArm = start;
        mEndArm = end;
    }

    public void setPonCount(int count) { mPonCount = count; }
    public int getPonCount() { return mPonCount; }

    public void setPonRegionCount(int count) { mPonRegionCount = count; }
    public int getPonRegionCount() { return mPonRegionCount; }

    public void setFragileSites(String typeStart, String typeEnd) { mStartFragileSite = typeStart; mEndFragileSite = typeEnd; }
    public String isStartFragileSite() { return mStartFragileSite; }
    public String isEndFragileSite() { return mEndFragileSite; }

    public void setLineElements(String typeStart, String typeEnd) { mStartLineElement = typeStart; mEndLineElement = typeEnd; }
    public String isStartLineElement() { return mStartLineElement; }
    public String isEndLineElement() { return mEndLineElement; }

    public int getNearestSVLength() { return mNearestSVLength; }
    public void setNearestSVLength(int length) { mNearestSVLength = length; }

    public final String getNearestSVLinkType() { return mNearestSVLinkType; }
    public void setNearestSVLinkType(String type) { mNearestSVLinkType = type; }

    public int getNearestTILength() { return mNearestTILength; }
    public void setNearestTILength(int length) { mNearestTILength = length; }

    public int getNearestDBLength() { return mNearestDBLength; }
    public void setNearestDBLength(int length) { mNearestDBLength = length; }

    public final String typeStr()
    {
        if(type() != StructuralVariantType.BND && mStartArm != mEndArm)
        {
            return "CRS";
        }
        else
        {
            return type().toString();
        }
    }

    public final boolean isLocal()
    {
        // means that both ends are within the same chromosomal arm
        return chromosome(true).equals(chromosome(false));
    }

    public final boolean isCrossArm() { return mStartArm != mEndArm; }

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

//    public void setStartCluster(SvCluster cluster) { mStartCluster = cluster; }
//    public void setEndCluster(SvCluster cluster) { mEndCluster = cluster; }
//    public final SvCluster getStartCluster() { return mStartCluster; }
//    public final SvCluster getEndCluster() { return mEndCluster; }

    // public final boolean areClustersSet() { return mStartCluster != null && mEndCluster != null; }
}
