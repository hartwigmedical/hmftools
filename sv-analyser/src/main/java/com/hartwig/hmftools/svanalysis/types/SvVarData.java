package com.hartwig.hmftools.svanalysis.types;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.round;

import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.BND;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.DEL;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.DUP;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.INS;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.SGL;
import static com.hartwig.hmftools.svanalysis.analysis.SvUtilities.appendStr;
import static com.hartwig.hmftools.svanalysis.annotators.FragileSiteAnnotator.NO_FS;
import static com.hartwig.hmftools.svanalysis.annotators.LineElementAnnotator.NO_LINE_ELEMENT;
import static com.hartwig.hmftools.svanalysis.types.SvLinkedPair.ASSEMBLY_MATCH_MATCHED;
import static com.hartwig.hmftools.svanalysis.types.SvLinkedPair.ASSEMBLY_MATCH_NONE;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantData;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantType;
import com.hartwig.hmftools.common.variant.structural.annotation.GeneAnnotation;
import com.hartwig.hmftools.common.variant.structural.annotation.Transcript;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class SvVarData
{
    private final String mIdStr; // sourced from DB so could be converted to int

    // full set of DB fields
    private final StructuralVariantData mSVData;
    private String mStartArm;
    private String mEndArm;
    private SvBreakend mBreakendStart;
    private SvBreakend mBreakendEnd;
    private String mStartFragileSite;
    private String mEndFragileSite;
    private String mStartLineElement;
    private String mEndLineElement;

    private String mAssemblyStartData;
    private String mAssemblyEndData;

    private boolean mDupBEStart;
    private boolean mDupBEEnd;

    private SvCluster mCluster;
    private String mClusterReason;

    private SvBreakend mFoldbackBeStart; // either the 2 breakends for this SV, or another SV's breaked
    private SvBreakend mFoldbackBeEnd;
    private int mFoldbackLenStart;
    private int mFoldbackLenEnd;
    private String mFoldbackInfoStart;
    private String mFoldbackInfoEnd;

    private long mNearestSvDistance;
    private String mNearestSvRelation;

    private SvLinkedPair mLinkStart; // templated insertion formed from this breakend to another
    private SvLinkedPair mLinkEnd;

    private SvLinkedPair mDBStart; // deletion bridge formed from this breakend to another
    private SvLinkedPair mDBEnd;
    private List<String> mStartTempInsertionAssemblies;
    private List<String> mEndTempInsertionAssemblies;
    private String mStartAssemblyMatchType;
    private String mEndAssemblyMatchType;
    private boolean mIsReplicatedSv;
    private SvVarData mReplicatedSv;
    private int mReplicatedCount;

    private String mConsecBEStart;
    private String mConsecBEEnd;

    private List<GeneAnnotation> mGenesStart;
    private List<GeneAnnotation> mGenesEnd;

    private String mDriverGeneStart;
    private String mDriverGeneEnd;

    private double mReplicationOriginStart;
    private double mReplicationOriginEnd;

    // copy number related data
    private double mCopyNumberStart; // cached from SV data but modifiable
    private double mCopyNumberEnd;
    private double mCopyNumberChangeStart;
    private double mCopyNumberChangeEnd;
    private double mPloidy;

    private boolean mHasCalcPloidy;
    private double mPloidyMin;
    private double mPloidyMax;
    private SvCNData mCnDataPrevStart; // segment leading to the start position
    private SvCNData mCnDataPostStart; // segment starting with the position position
    private SvCNData mCnDataPrevEnd;
    private SvCNData mCnDataPostEnd;

    public static final String NONE_SEGMENT_INFERRED = "INFERRED";

    public static String ASSEMBLY_TYPE_DSB = "dsb";
    public static String ASSEMBLY_TYPE_TI = "asm";
    public static String ASSEMBLY_TYPE_EQV = "eqv";

    public static String RELATION_TYPE_NEIGHBOUR = "NHBR";
    public static String RELATION_TYPE_OVERLAP = "OVRL";

    // iterators for start and end data
    public static int SVI_START = 0;
    public static int SVI_END = 1;

    private static final Logger LOGGER = LogManager.getLogger(SvVarData.class);

    public SvVarData(final StructuralVariantData svData)
    {
        mIdStr = svData.id();

        mSVData = svData;

        init();

        setAssemblyData(false);
    }

    private void init()
    {
        mStartArm = "";
        mEndArm = "";
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

        mClusterReason = "";
        mCluster = null;

        mLinkStart = null;
        mLinkEnd = null;
        mDBStart = null;
        mDBEnd = null;

        mFoldbackBeStart = null;
        mFoldbackBeEnd = null;
        mFoldbackLenStart = -1;
        mFoldbackLenEnd = -1;
        mFoldbackInfoStart = "";
        mFoldbackInfoEnd = "";

        mConsecBEStart = "";
        mConsecBEEnd = "";

        mGenesStart = Lists.newArrayList();
        mGenesEnd = Lists.newArrayList();

        mDriverGeneStart = "";
        mDriverGeneEnd = "";

        mReplicationOriginStart = 0;
        mReplicationOriginEnd = 0;

        mCopyNumberStart = mSVData.adjustedStartCopyNumber();
        mCopyNumberEnd = mSVData.adjustedEndCopyNumber();
        mCopyNumberChangeStart = mSVData.adjustedStartCopyNumberChange();
        mCopyNumberChangeEnd = mSVData.adjustedEndCopyNumberChange();
        mPloidy = mSVData.ploidy();

        mHasCalcPloidy = false;
        mPloidyMin = 0;
        mPloidyMax = 0;
        mCnDataPrevStart = null;
        mCnDataPostStart = null;
        mCnDataPrevEnd = null;
        mCnDataPostEnd = null;
    }

    public SvVarData(final SvVarData other)
    {
        mSVData = other.getSvData();

        init();

        mIdStr = other.getSvData().id() + "r";
        mStartArm = other.arm(true);
        mEndArm = other.arm(false);

        mBreakendStart = new SvBreakend(this, true);
        mBreakendStart.setChrPosIndex(other.getBreakend(true).getChrPosIndex());

        if(!isNullBreakend())
        {
            mBreakendEnd = new SvBreakend(this, false);
            mBreakendEnd.setChrPosIndex(other.getBreakend(false).getChrPosIndex());
        }

        mStartFragileSite = other.isFragileSite(true);
        mEndFragileSite = other.isFragileSite(false);
        mStartLineElement = other.getLineElement(true);
        mEndLineElement = other.getLineElement(false);
        mNearestSvDistance = other.getNearestSvDistance();
        mNearestSvRelation = other.getNearestSvRelation();

        mAssemblyStartData = other.getAssemblyData(true);
        mAssemblyEndData = other.getAssemblyData(false);
        setAssemblyData(true);

        mStartAssemblyMatchType = other.getAssemblyMatchType(true);
        mEndAssemblyMatchType = other.getAssemblyMatchType(false);
        mDupBEStart = other.isDupBreakend(true);
        mDupBEEnd = other.isDupBreakend(false);
        mIsReplicatedSv = true;
        mReplicatedSv = other;
        mClusterReason = other.getClusterReason();
        mCluster = other.getCluster();
    }

    public final String id() { return mIdStr; }

    public final int dbId()
    {
        return Integer.parseInt(mSVData.id());
    }

    public final StructuralVariantData getSvData() { return mSVData; }
    public boolean isNoneSegment() { return mSVData.filter().equals(NONE_SEGMENT_INFERRED); }

    // for convenience
    public final String chromosome(boolean isStart) { return isStart ? mSVData.startChromosome() : mSVData.endChromosome(); }
    public final long position(boolean isStart) { return isStart ? mSVData.startPosition() : mSVData.endPosition(); }
    public final byte orientation(boolean isStart){ return isStart ? mSVData.startOrientation() : mSVData.endOrientation(); }
    public final double copyNumber(boolean isStart){ return isStart ? mCopyNumberStart : mCopyNumberEnd; }
    public final StructuralVariantType type() { return mSVData.type(); }

    public SvBreakend getBreakend(boolean isStart) { return isStart ? mBreakendStart : mBreakendEnd; }

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
        return String.format("%s: %s %s:%d:%d",
                id(), useStart ? "start" :"end", chromosome(useStart), orientation(useStart), position(useStart));
    }

    public final String arm(boolean isStart) { return isStart ? mStartArm : mEndArm; }
    public void setChromosomalArms(final String start, final String end)
    {
        mStartArm = start;
        mEndArm = end;

        mBreakendStart = new SvBreakend(this, true);

        if(!isNullBreakend())
            mBreakendEnd = new SvBreakend(this, false);
    }

    public final SvCluster getCluster() { return mCluster; }
    public void setCluster(final SvCluster cluster)
    {
        // isSpecificSV(this);
        mCluster = cluster;
    }

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

        mClusterReason = appendStr(mClusterReason, reason, ';');

        if(!otherId.isEmpty())
            mClusterReason += "_" + otherId;

        if(otherId.equals(mIdStr))
        {
            LOGGER.warn("SV({}) reason({}) setting to own ID", mIdStr, reason);
        }
    }

    public final String getClusterReason() { return mClusterReason; }

    public boolean isReplicatedSv() { return mIsReplicatedSv; }
    public final SvVarData getReplicatedSv() { return mReplicatedSv; }
    public final SvVarData getOrigSV() { return mIsReplicatedSv ? mReplicatedSv : this; }
    public int getReplicatedCount() { return mReplicatedCount; }
    public void setReplicatedCount(int count) { mReplicatedCount = count; }
    public final String origId() { return getOrigSV().id(); }
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

    public double ploidy()
    {
        // use the estimated ploidy when present
        return mHasCalcPloidy ? (mPloidyMax + mPloidyMin) * 0.5 : mPloidy;
    }

    public double copyNumberChange(boolean isStart)
    {
        // TEMP: precise DBs cause incorrect copy number change, so in this case use ploidy
        if(isStart)
        {
            if(mDBStart != null && mDBStart.length() == 0)
                return mPloidy;
            else
                return mCopyNumberChangeStart;
        }
        else
        {
            if(mDBEnd != null && mDBEnd.length() == 0)
                return mPloidy;
            else
                return mCopyNumberChangeEnd;
        }
    }

    public void setCopyNumberData(boolean isStart, double copyNumber, double copyNumberChange)
    {
        if(isStart)
        {
            mCopyNumberStart = copyNumber;
            mCopyNumberChangeStart = copyNumberChange;
        }
        else
        {
            mCopyNumberEnd = copyNumber;
            mCopyNumberChangeEnd = copyNumberChange;
        }
    }

    public static double SUSPECT_CN_CHANGE = 0.2;
    private static double CN_ROUND_SIZE = 0.5;

    public double getRoundedCNChange()
    {
        double cnChgStart = copyNumberChange(true);
        double cnChange;

        if(!isNullBreakend())
        {
            // TEMP: before more rebust CN Change correction logic introduced
            double cnChgEnd = copyNumberChange(false);

            if(cnChgStart > SUSPECT_CN_CHANGE * 2 && cnChgEnd < SUSPECT_CN_CHANGE)
                cnChange = cnChgStart;
            else if(cnChgEnd > SUSPECT_CN_CHANGE * 2 && cnChgStart < SUSPECT_CN_CHANGE)
                cnChange = cnChgEnd;
            else
                cnChange = (cnChgStart + cnChgEnd) * 0.5;
        }
        else
        {
            cnChange = cnChgStart;
        }

        if(cnChange >= 0.6)
            return round(cnChange);
        else
            return round(cnChange/CN_ROUND_SIZE) * CN_ROUND_SIZE;
    }

    public boolean hasInconsistentCopyNumberChange(boolean useStart)
    {
        return ploidy() - copyNumberChange(useStart) > 0.8;
    }

    public long getNearestSvDistance() { return mNearestSvDistance; }
    public void setNearestSvDistance(long distance) { mNearestSvDistance = distance; }
    public String getNearestSvRelation() { return mNearestSvRelation; }
    public void setNearestSvRelation(final String rel) { mNearestSvRelation = rel; }

    public void setFragileSites(String typeStart, String typeEnd) { mStartFragileSite = typeStart; mEndFragileSite = typeEnd; }
    public String isFragileSite(boolean useStart) { return useStart ? mStartFragileSite : mEndFragileSite; }

    public void setLineElement(String type, boolean isStart)
    {
        if(isStart)
        {
            if(mStartLineElement.contains(type))
                return;
            else if(mStartLineElement.equals(NO_LINE_ELEMENT))
                mStartLineElement = type;
            else
                mStartLineElement = mStartLineElement + ";" + type;
        }
        else
        {
            if(mEndLineElement.contains(type))
                return;
            else if(mEndLineElement.equals(NO_LINE_ELEMENT))
                mEndLineElement = type;
            else
                mEndLineElement = mEndLineElement + ";" + type;
        }
    }

    public boolean isLineElement(boolean useStart)
    {
        if(useStart)
            return !mStartLineElement.equals(NO_LINE_ELEMENT);
        else
            return !mEndLineElement.equals(NO_LINE_ELEMENT);
    }

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

    public final SvLinkedPair getLinkedPair(boolean isStart) { return isStart ? mLinkStart : mLinkEnd; }

    public void setLinkedPair(final SvLinkedPair link, boolean isStart)
    {
        if(isStart)
            mLinkStart = link;
        else
            mLinkEnd = link;
    }

    public final SvLinkedPair getDBLink(boolean isStart) { return isStart ? mDBStart : mDBEnd; }

    public void setDBLink(final SvLinkedPair link, boolean isStart)
    {
        if(isStart)
            mDBStart = link;
        else
            mDBEnd = link;
    }

    public final String getFoldbackLink(boolean useStart)
    {
        if(useStart && mFoldbackBeStart != null)
            return mFoldbackBeStart.getSV().id();
        else if(!useStart && mFoldbackBeEnd != null)
            return mFoldbackBeEnd.getSV().id();
        else
            return "";
    }

    public final SvBreakend getFoldbackBreakend(boolean useStart) { return useStart ? mFoldbackBeStart : mFoldbackBeEnd; }
    public int getFoldbackLength(boolean useStart) { return useStart ? mFoldbackLenStart : mFoldbackLenEnd; }
    public final String getFoldbackInfo(boolean useStart) { return useStart ? mFoldbackInfoStart : mFoldbackInfoEnd; }
    public boolean isFoldback() { return mFoldbackBeEnd != null || mFoldbackBeStart != null; }
    public boolean isChainedFoldback()
    {
        if(mFoldbackBeEnd != null && mFoldbackBeEnd != mBreakendStart)
            return true;
        else if(mFoldbackBeStart != null && mFoldbackBeStart != mBreakendEnd)
            return true;
        else
            return false;
    }

    public final SvBreakend getChainedFoldbackBreakend()
    {
        return mFoldbackBeStart != null ? mFoldbackBeStart : mFoldbackBeEnd;
    }

    public void setFoldbackLink(boolean isStart, final SvBreakend link, int length, String linkInfo)
    {
        if(isStart)
        {
            mFoldbackBeStart = link;
            mFoldbackLenStart = length;
            mFoldbackInfoStart = linkInfo;
        }
        else
        {
            mFoldbackBeEnd = link;
            mFoldbackLenEnd = length;
            mFoldbackInfoEnd = linkInfo;
        }

        if(mCluster != null)
        {
            if (mFoldbackBeStart == null && mFoldbackBeEnd == null)
            {
                mCluster.deregisterFoldback(this);
            }
            else
            {
                mCluster.registerFoldback(this);
            }
        }
    }

    public void setFoldbackInfo(boolean isStart, final String info)
    {
        if(isStart)
            mFoldbackInfoStart = info;
        else
            mFoldbackInfoEnd = info;
    }

    public final SvVarData getChainedFoldbackSv()
    {
        if(!isChainedFoldback())
            return this;

        return mFoldbackBeStart != null ? mFoldbackBeStart.getSV() : mFoldbackBeEnd.getSV();
    }

    public String getConsecBEStart(boolean useStart) { return useStart ? mConsecBEStart : mConsecBEEnd; }
    public void setConsecBEStart(final String info, boolean useStart)
    {
        if(useStart)
        {
            mConsecBEStart = mConsecBEStart.isEmpty() ? info : mConsecBEStart + ";" + info;
        }
        else
        {
            mConsecBEEnd = mConsecBEEnd.isEmpty() ? info : mConsecBEEnd + ";" + info;
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

    public final boolean isCrossArm()
    {
        return type() != SGL && !isLocal();
    }

    public final boolean isSimpleType()
    {
        return (type() == DEL || type() == DUP || type() == INS);
    }

    public static boolean isStart(int svIter) { return svIter == SVI_START; }

    public String getAssemblyData(boolean useStart) { return useStart ? mAssemblyStartData : mAssemblyEndData; }

    // unit testing only
    public void setAssemblyData(boolean useStart, final String data)
    {
        if (useStart)
            mAssemblyStartData = data;
        else
            mAssemblyEndData = data;

        setAssemblyData(true);
    }

    public final List<String> getTempInsertionAssemblies(boolean useStart) { return useStart ? mStartTempInsertionAssemblies : mEndTempInsertionAssemblies; }

    public final List<GeneAnnotation> getGenesList(boolean useStart) { return useStart ? mGenesStart : mGenesEnd; }
    public void setGenesList(final List<GeneAnnotation> genesList, boolean isStart)
    {
        if(isStart)
            mGenesStart.addAll(genesList);
        else
            mGenesEnd.addAll(genesList);
    }

    public final String getDriverGene(boolean useStart)
    {
        return useStart ? mDriverGeneStart : mDriverGeneEnd;
    }

    public void setDriveGene(final String geneInfo, boolean isStart)
    {
        if(isStart)
            mDriverGeneStart = geneInfo;
        else
            mDriverGeneEnd = geneInfo;
    }

    public final String getGeneInBreakend(boolean useStart)
    {
        final List<GeneAnnotation> genesList = getGenesList(useStart);

        if(genesList.isEmpty())
            return "";

        Transcript longestTrans = null;
        for(final GeneAnnotation gene : genesList)
        {
            for(final Transcript trans : gene.transcripts())
            {
                if(!trans.isIntronic() && !trans.isExonic())
                    continue;

                if(longestTrans == null || trans.length() > longestTrans.length())
                    longestTrans = trans;
            }
        }

        return longestTrans != null ? longestTrans.parent().GeneName : "";
    }

    public final String getAssemblyMatchType(boolean useStart) { return useStart ? mStartAssemblyMatchType : mEndAssemblyMatchType; }
    public boolean isAssemblyMatched(boolean useStart) { return getAssemblyMatchType(useStart).equals(ASSEMBLY_MATCH_MATCHED); }

    public void setAssemblyMatchType(String type, boolean useStart)
    {
        if(useStart)
            mStartAssemblyMatchType = type;
        else
            mEndAssemblyMatchType = type;
    }

    private void setAssemblyData(boolean useExisting)
    {
        mStartTempInsertionAssemblies = Lists.newArrayList();
        mEndTempInsertionAssemblies = Lists.newArrayList();

        mStartAssemblyMatchType = ASSEMBLY_MATCH_NONE;
        mEndAssemblyMatchType = ASSEMBLY_MATCH_NONE;

        if(!useExisting)
        {
            mAssemblyStartData = "";
            mAssemblyEndData = "";

            if(!mSVData.startLinkedBy().isEmpty() && !mSVData.startLinkedBy().equals("."))
            {
                mAssemblyStartData = mSVData.startLinkedBy().replaceAll(",", ";");
            }

            if(!mSVData.endLinkedBy().isEmpty() && !mSVData.endLinkedBy().equals("."))
            {
                mAssemblyEndData = mSVData.endLinkedBy().replaceAll(",", ";");
            }
        }

        if(!mAssemblyStartData.isEmpty())
        {
            String[] assemblyList = mAssemblyStartData.split(";");

            for(int i = 0; i < assemblyList.length; ++i)
            {
                if(assemblyList[i].contains(ASSEMBLY_TYPE_TI))
                    mStartTempInsertionAssemblies.add(assemblyList[i]);
            }
        }

        if(!mAssemblyEndData.isEmpty())
        {
            String[] assemblyList = mAssemblyEndData.split(";");
            for(int i = 0; i < assemblyList.length; ++i)
            {
                if(assemblyList[i].contains(ASSEMBLY_TYPE_TI))
                    mEndTempInsertionAssemblies.add(assemblyList[i]);
            }
        }
    }

    public double getReplicationOrigin(boolean isStart) { return isStart ? mReplicationOriginStart : mReplicationOriginEnd; }
    public void setReplicationOrigin(boolean isStart, double value)
    {
        if(isStart)
            mReplicationOriginStart = value;
        else
            mReplicationOriginEnd = value;
    }

    public void setPloidyRecalcData(double minPloidy, double maxPloidy)
    {
        mPloidyMin = minPloidy;
        mPloidyMax = maxPloidy;
        mHasCalcPloidy = true;
    }

    public boolean hasCalculatedPloidy() { return mHasCalcPloidy; }
    public double ploidyMax() { return mHasCalcPloidy ? mPloidyMax : mPloidy; }
    public double ploidyMin() { return mHasCalcPloidy ? mPloidyMin : mPloidy; }

    public final SvCNData getCopyNumberData(boolean isStart, boolean isPrevious)
    {
        return isStart ? (isPrevious ? mCnDataPrevStart : mCnDataPostStart) : (isPrevious ? mCnDataPrevEnd : mCnDataPostEnd);
    }

    public void setCopyNumberData(boolean isStart, final SvCNData prevData, final SvCNData postData)
    {
        if(isStart)
        {
            mCnDataPrevStart = prevData;
            mCnDataPostStart = postData;
        }
        else
        {
            mCnDataPrevEnd = prevData;
            mCnDataPostEnd = postData;
        }
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

    // private static String SPECIFIC_VAR_ID = "9559903";
    private static String SPECIFIC_VAR_ID = "";

    public static boolean isSpecificSV(final SvVarData var)
    {
        if(var.id().equals(SPECIFIC_VAR_ID))
            return true;

        return false;
    }

}
