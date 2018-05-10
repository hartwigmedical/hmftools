package com.hartwig.hmftools.svanalysis.types;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.structural.EnrichedStructuralVariant;
import com.hartwig.hmftools.common.variant.structural.ImmutableStructuralVariantData;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantData;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantLeg;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantType;

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

    private int mNearestSVLength;
    private String mNearestSVLinkType;
    private int mNearestTILength; // templated insertion if exists
    private int mNearestDBLength; // deletion-bridge link if exists

    private boolean mDupBEStart;
    private boolean mDupBEEnd;

    private String mTransType;
    private int mTransLength;
    private String mTransSvLinks;

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

        mDupBEStart = false;
        mDupBEEnd = false;

        mTransType = "";
        mTransLength = 0;
        mTransSvLinks = "";

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

    public final String posId(boolean useStart)
    {
        return String.format("%s: %s %s:%d:%d)",
                id(), useStart ? "start" :"end", chromosome(useStart), orientation(useStart), position(useStart));
    }

    public final String toCsv()
    {
        return String.format("%s,%s,%d,%d,%s,%d,%d,%s)",
                id(), chromosome(true), orientation(true), position(true),
                chromosome(false), orientation(false), position(false), type());
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

    public boolean isDupBEStart() { return mDupBEStart; }
    public boolean isDupBEEnd() { return mDupBEEnd; }
    public void setIsDupBEStart(boolean toggle) { mDupBEStart = toggle; }
    public void setIsDupBEEnd(boolean toggle) { mDupBEEnd = toggle; }

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

    public String getTransType() { return mTransType; }
    public int getTransLength() { return mTransLength; }
    public String getTransSvLinks() { return mTransSvLinks; }

    public void setTransData(String transType, int transLength, final String transSvLinks)
    {
        mTransType = transType;
        mTransLength = transLength;
        mTransSvLinks = transSvLinks;
    }

}
