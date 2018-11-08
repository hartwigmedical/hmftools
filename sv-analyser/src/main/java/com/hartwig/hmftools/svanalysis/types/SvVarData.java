package com.hartwig.hmftools.svanalysis.types;

import static java.lang.Math.abs;
import static java.lang.Math.round;

import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.BND;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.DEL;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.DUP;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.INS;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.SGL;
import static com.hartwig.hmftools.svanalysis.annotators.FragileSiteAnnotator.NO_FS;
import static com.hartwig.hmftools.svanalysis.annotators.LineElementAnnotator.NO_LINE_ELEMENT;
import static com.hartwig.hmftools.svanalysis.types.SvLinkedPair.ASSEMBLY_MATCH_MATCHED;
import static com.hartwig.hmftools.svanalysis.types.SvLinkedPair.ASSEMBLY_MATCH_NONE;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.structural.EnrichedStructuralVariant;
import com.hartwig.hmftools.common.variant.structural.ImmutableStructuralVariantData;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantData;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantType;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class SvVarData
{
    private final String mId; // sourced from either VCF or DB

    // full set of DB fields
    private final StructuralVariantData mSVData;
    private boolean mNoneSegment; // created from a NONE copy number segment
    private String mStartArm;
    private String mEndArm;
    private int mPonCount;
    private int mPonRegionCount; // allowing for a small buffer either side of a PON
    private String mStartFragileSite;
    private String mEndFragileSite;
    private String mStartLineElement;
    private String mEndLineElement;

    private String mAssemblyStartData;
    private String mAssemblyEndData;

    private boolean mDupBEStart;
    private boolean mDupBEEnd;

    private String mTransType;
    private int mTransLength;
    private String mTransSvLinks;

    private SvCluster mCluster;
    private String mClusterReason;

    private String mFoldbackLinkStart;
    private String mFoldbackLinkEnd;
    private int mFoldbackLenStart;
    private int mFoldbackLenEnd;

    private long mNearestSvDistance;
    private String mNearestSvRelation;

    private SvLinkedPair mStartLink;
    private SvLinkedPair mEndLink;
    private List<String> mStartTempInsertionAssemblies;
    private List<String> mStartDsbAssemblies;
    private List<String> mStartOtherAssemblies;
    private List<String> mEndTempInsertionAssemblies;
    private List<String> mEndDsbAssemblies;
    private List<String> mEndOtherAssemblies;
    private String mStartAssemblyMatchType;
    private String mEndAssemblyMatchType;
    private boolean mIsReplicatedSv;
    private SvVarData mReplicatedSv;
    private int mReplicatedCount;

    public static String ASSEMBLY_TYPE_DSB = "dsb";
    public static String ASSEMBLY_TYPE_TI = "asm";
    public static String ASSEMBLY_TYPE_OTHER = "bpb";
    public static String ASSEMBLY_TYPE_EQV = "eqv";

    public static String RELATION_TYPE_NEIGHBOUR = "NHBR";
    public static String RELATION_TYPE_OVERLAP = "OVRL";

    // iterators for start and end data
    public static int SVI_START = 0;
    public static int SVI_END = 1;

    private static final Logger LOGGER = LogManager.getLogger(SvVarData.class);

    public SvVarData(final StructuralVariantData svData)
    {
        mId = svData.id();
        mSVData = svData;

        init();

        setAssemblyData(mSVData.startLinkedBy(), mSVData.endLinkedBy());
    }

    private void init()
    {
        mNoneSegment = false;
        mStartArm = "";
        mEndArm = "";
        mPonCount = 0;
        mPonRegionCount = 0;
        mStartFragileSite = NO_FS;
        mEndFragileSite = NO_FS;
        mStartLineElement = NO_LINE_ELEMENT;
        mEndLineElement = NO_LINE_ELEMENT;

        mNearestSvDistance = -1;
        mNearestSvRelation = "";

        mIsReplicatedSv = false;
        mReplicatedSv = null;
        mReplicatedCount = 0;

        mDupBEStart = false;
        mDupBEEnd = false;

        mTransType = "";
        mTransLength = 0;
        mTransSvLinks = "";

        mClusterReason = "";
        mCluster = null;

        mStartLink = null;
        mEndLink = null;

        mFoldbackLinkStart = "";
        mFoldbackLinkEnd = "";
        mFoldbackLenStart = -1;
        mFoldbackLenEnd = -1;
    }

    public static SvVarData from(final EnrichedStructuralVariant enrichedSV)
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

        return new SvVarData(svData);
    }

    public SvVarData(final SvVarData other)
    {
        init();

        mId = other.getSvData().id() + "r";
        mSVData = other.getSvData();
        mNoneSegment = other.isNoneSegment();
        mStartArm = other.arm(true);
        mEndArm = other.arm(false);
        mPonCount = other.getPonCount();
        mPonRegionCount = other.getPonRegionCount();
        mStartFragileSite = other.isFragileSite(true);
        mEndFragileSite = other.isFragileSite(false);
        mStartLineElement = other.getLineElement(true);
        mEndLineElement = other.getLineElement(false);
        mNearestSvDistance = other.getNearestSvDistance();
        mNearestSvRelation = other.getNearestSvRelation();
        setAssemblyData(mSVData.startLinkedBy(), mSVData.endLinkedBy());
        mStartAssemblyMatchType = other.getAssemblyMatchType(true);
        mEndAssemblyMatchType = other.getAssemblyMatchType(false);
        mDupBEStart = other.isDupBreakend(true);
        mDupBEEnd = other.isDupBreakend(false);
        mIsReplicatedSv = true;
        mReplicatedSv = other;
        mClusterReason = other.getClusterReason();
        mCluster = other.getCluster();
    }

    public final String id() { return mId; }
    public final StructuralVariantData getSvData() { return mSVData; }
    public void setNoneSegment(boolean toggle) { mNoneSegment = toggle; }
    public boolean isNoneSegment() { return mNoneSegment; }

    // for convenience
    public boolean equals(final SvVarData other) { return id().equals(other.id()); }

    public final String chromosome(boolean isStart) { return isStart ? mSVData.startChromosome() : mSVData.endChromosome(); }
    public final long position(boolean isStart) { return isStart ? mSVData.startPosition() : mSVData.endPosition(); }
    public final byte orientation(boolean isStart){ return isStart ? mSVData.startOrientation() : mSVData.endOrientation(); }
    public final double copyNumber(boolean isStart){ return isStart ? mSVData.adjustedStartCopyNumber() : mSVData.adjustedEndCopyNumber(); }
    public final double copyNumberChange(boolean isStart){ return isStart ? mSVData.adjustedStartCopyNumberChange() : mSVData.adjustedEndCopyNumberChange(); }
    public final StructuralVariantType type() { return mSVData.type(); }

    public boolean isNullBreakend() { return type() == SGL; }

    public final String posId()
    {
        if(isNullBreakend())
        {
            return String.format("id(%s) pos(%s:%d:%d)",
                    id(), chromosome(true), orientation(true), position(true));
        }
        else
        {
            return String.format("id(%s) pos(%s:%d:%d -> %s:%d:%d)",
                    id(), chromosome(true), orientation(true), position(true),
                    chromosome(false), orientation(false), position(false));
        }
    }

    public final String posId(boolean useStart)
    {
        return String.format("%s: %s %s:%d:%d)",
                id(), useStart ? "start" :"end", chromosome(useStart), orientation(useStart), position(useStart));
    }

    public final String arm(boolean isStart) { return isStart ? mStartArm : mEndArm; }
    public void setChromosomalArms(final String start, final String end)
    {
        mStartArm = start;
        mEndArm = end;
    }

    public final SvCluster getCluster() { return mCluster; }
    public void setCluster(final SvCluster cluster) { mCluster = cluster; }

    public final long length()
    {
        if(type() == BND || isNullBreakend())
            return 0;

        return abs(position(false) - position(true));
    }

    public void addClusterReason(final String reason, final String otherId)
    {
        if(mClusterReason.contains(reason))
            return;

        if(!mClusterReason.isEmpty())
            mClusterReason += ";";

        mClusterReason += reason;

        if(!otherId.isEmpty())
            mClusterReason += "_" + otherId;

        if(otherId.equals(mId))
        {
            LOGGER.warn("SV({}) reason({}) setting to own ID", mId, reason);
        }
    }

    public final String getClusterReason() { return mClusterReason; }

    public boolean isReplicatedSv() { return mIsReplicatedSv; }
    public final SvVarData getReplicatedSv() { return mReplicatedSv; }
    public int getReplicatedCount() { return mReplicatedCount; }
    public void setReplicatedCount(int count) { mReplicatedCount = count; }
    public final String origId() { return mReplicatedSv != null ? mReplicatedSv.id() : mId; }
    public boolean equals(final SvVarData other, boolean allowReplicated)
    {
        if(this == other)
            return true;

        if(allowReplicated)
        {
            if(this == other.getReplicatedSv() || mReplicatedSv == other)
                return true;

            if(mReplicatedSv != null && mReplicatedSv == other.getReplicatedSv())
                return true;
        }

        return false;
    }

    public int impliedCopyNumber(boolean useStart)
    {
        return (int)round(copyNumberChange(useStart));
    }

    public boolean hasInconsistentCopyNumberChange(boolean useStart)
    {
        return mSVData.ploidy() - copyNumberChange(useStart) > 0.8;
        // return round(copyNumberChange(true)) != round(copyNumberChange(false));
    }

    public long getNearestSvDistance() { return mNearestSvDistance; }
    public void setNearestSvDistance(long distance) { mNearestSvDistance = distance; }
    public String getNearestSvRelation() { return mNearestSvRelation; }
    public void setNearestSvRelation(final String rel) { mNearestSvRelation = rel; }

    public void setPonCount(int count) { mPonCount = count; }
    public int getPonCount() { return mPonCount; }

    public void setPonRegionCount(int count) { mPonRegionCount = count; }
    public int getPonRegionCount() { return mPonRegionCount; }

    public void setFragileSites(String typeStart, String typeEnd) { mStartFragileSite = typeStart; mEndFragileSite = typeEnd; }
    public String isFragileSite(boolean useStart) { return useStart ? mStartFragileSite : mEndFragileSite; }

    public void setLineElement(String type, boolean isStart)
    {
        if(isStart)
            mStartLineElement = type;
        else
            mEndLineElement = type;
    }

    public boolean isLineElement(boolean useStart) { return useStart ? !mStartLineElement.equals(NO_LINE_ELEMENT) : !mEndLineElement.equals(NO_LINE_ELEMENT); }
    public boolean inLineElement() { return isLineElement(true) || isLineElement(false); }
    public final String getLineElement(boolean useStart) { return useStart ? mStartLineElement : mEndLineElement; }

    public boolean isDupBreakend(boolean useStart) { return useStart ? mDupBEStart : mDupBEEnd; }

    public void setIsDupBreakend(boolean toggle, boolean isStart)
    {
        if(isStart)
            mDupBEStart = toggle;
        else
            mDupBEEnd = toggle;
    }

    public final SvLinkedPair getLinkedPair(boolean isStart) { return isStart ? mStartLink : mEndLink; }

    public void setLinkedPair(final SvLinkedPair link, boolean isStart)
    {
        if(isStart)
            mStartLink = link;
        else
            mEndLink = link;
    }

    public String getFoldbackLink(boolean useStart) { return useStart ? mFoldbackLinkStart : mFoldbackLinkEnd; }
    public int getFoldbackLen(boolean useStart) { return useStart ? mFoldbackLenStart : mFoldbackLenEnd; }
    public void setFoldbackLink(boolean isStart, String link, int len)
    {
        if(isStart)
        {
            mFoldbackLinkStart = link;
            mFoldbackLenStart = len;
        }
        else
        {
            mFoldbackLinkEnd = link;
            mFoldbackLenEnd = len;
        }

        if(mCluster != null)
        {
            if (mFoldbackLinkStart.isEmpty() && mFoldbackLinkEnd.isEmpty())
            {
                mCluster.deregisterFoldback(this);
            }
            else
            {
                mCluster.registerFoldback(this);
            }
        }
    }

    public final String typeStr()
    {
        //if(chromosome(true).equals(chromosome(false)) && mStartArm != mEndArm )
        //    return "CRS";
        if(isNoneSegment())
            return "NONE";
        else
            return type().toString();
    }

    public final boolean isLocal()
    {
        // means that both ends are within the same chromosomal arm
        return chromosome(true).equals(chromosome(false)) && mStartArm.equals(mEndArm);
    }

    public final boolean isSimpleType()
    {
        return (type() == DEL || type() == DUP || type() == INS);
    }

    public static boolean isTranslocation(StructuralVariantType type) { return type == BND; }
    public boolean isTranslocation()
    {
        return isTranslocation(type());
    }

    public static boolean isStart(int svIter) { return svIter == SVI_START; }

    public String getTransType() { return mTransType; }
    public int getTransLength() { return mTransLength; }
    public String getTransSvLinks() { return mTransSvLinks; }

    public void setTransData(String transType, int transLength, final String transSvLinks)
    {
        mTransType = transType;
        mTransLength = transLength;
        mTransSvLinks = transSvLinks;
    }

    public String getAssemblyData(boolean useStart) { return useStart ? mAssemblyStartData : mAssemblyEndData; }

    public List<String> getTempInsertionAssemblies(boolean useStart) { return useStart ? mStartTempInsertionAssemblies : mEndTempInsertionAssemblies; }
    public List<String> getDsbAssemblies(boolean useStart) { return useStart ? mStartDsbAssemblies : mEndDsbAssemblies; }
    public List<String> getOtherAssemblies(boolean useStart) { return useStart ? mStartOtherAssemblies : mEndOtherAssemblies; }

    public final String getAssemblyMatchType(boolean useStart) { return useStart ? mStartAssemblyMatchType : mEndAssemblyMatchType; }
    public boolean isAssemblyMatched(boolean useStart) { return getAssemblyMatchType(useStart).equals(ASSEMBLY_MATCH_MATCHED); }

    public void setAssemblyMatchType(String type, boolean useStart)
    {
        if(useStart)
            mStartAssemblyMatchType = type;
        else
            mEndAssemblyMatchType = type;
    }

    private void setAssemblyData(final String assemblyStart, final String assemblyEnd)
    {
        mStartTempInsertionAssemblies = Lists.newArrayList();
        mStartDsbAssemblies = Lists.newArrayList();
        mStartOtherAssemblies = Lists.newArrayList();
        mEndTempInsertionAssemblies = Lists.newArrayList();
        mEndDsbAssemblies = Lists.newArrayList();
        mEndOtherAssemblies = Lists.newArrayList();
        mAssemblyStartData = "";
        mAssemblyEndData = "";
        mStartAssemblyMatchType = ASSEMBLY_MATCH_NONE;
        mEndAssemblyMatchType = ASSEMBLY_MATCH_NONE;

        if(!mSVData.startLinkedBy().isEmpty() && !mSVData.startLinkedBy().equals("."))
        {
            mAssemblyStartData = mSVData.startLinkedBy().replaceAll(",", ";");

            String[] assemblyList = mAssemblyStartData.split(";");

            for(int i = 0; i < assemblyList.length; ++i)
            {
                if(assemblyList[i].contains(ASSEMBLY_TYPE_TI))
                    mStartTempInsertionAssemblies.add(assemblyList[i]);
                else if(assemblyList[i].contains(ASSEMBLY_TYPE_DSB))
                    mStartDsbAssemblies.add(assemblyList[i]);
                else
                    mStartOtherAssemblies.add(assemblyList[i]);
            }
        }

        if(!mSVData.endLinkedBy().isEmpty() && !mSVData.endLinkedBy().equals("."))
        {
            mAssemblyEndData = mSVData.endLinkedBy().replaceAll(",", ";");

            String[] assemblyList = mAssemblyEndData.split(";");
            for(int i = 0; i < assemblyList.length; ++i)
            {
                if(assemblyList[i].contains(ASSEMBLY_TYPE_TI))
                    mEndTempInsertionAssemblies.add(assemblyList[i]);
                else if(assemblyList[i].contains(ASSEMBLY_TYPE_DSB))
                    mEndDsbAssemblies.add(assemblyList[i]);
                else
                    mEndOtherAssemblies.add(assemblyList[i]);
            }
        }
    }

    public static boolean haveLinkedAssemblies(final SvVarData var1, final SvVarData var2, boolean v1Start, boolean v2Start)
    {
        for(String assemb1 : var1.getTempInsertionAssemblies(v1Start))
        {
            if(var2.getTempInsertionAssemblies(v2Start).contains(assemb1))
                return true;
        }

        return false;
    }

    public static SvVarData findVariantById(final String id, List<SvVarData> svList)
    {
        for(SvVarData var : svList)
        {
            if(var.id().equals(id))
                return var;
        }

        return null;
    }

    public static boolean haveSameChrArms(final SvVarData var1, final SvVarData var2)
    {
        // tests if 2 variants (including BNDs) link the same 2 chromosomal arms
        if(var1.chromosome(true).equals(var2.chromosome(true)) && var1.chromosome(false).equals(var2.chromosome(false))
        && var1.arm(true).equals(var2.arm(true)) && var1.arm(false).equals(var2.arm(false)))
        {
            return true;
        }
        else if(var1.chromosome(true).equals(var2.chromosome(false)) && var1.chromosome(false).equals(var2.chromosome(true))
        && var1.arm(true).equals(var2.arm(false)) && var1.arm(false).equals(var2.arm(true)))
        {
            return true;
        }

        return false;
    }

}
