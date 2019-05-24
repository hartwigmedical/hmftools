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
import static com.hartwig.hmftools.svanalysis.annotators.LineElementAnnotator.NO_LINE_ELEMENT;
import static com.hartwig.hmftools.svanalysis.types.SvLinkedPair.ASSEMBLY_MATCH_MATCHED;
import static com.hartwig.hmftools.svanalysis.types.SvLinkedPair.ASSEMBLY_MATCH_NONE;
import static com.hartwig.hmftools.svanalysis.types.SvaConfig.SPECIFIC_SV_ID;
import static com.hartwig.hmftools.svanalysis.types.SvaConstants.MIN_TEMPLATED_INSERTION_LENGTH;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantData;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantType;
import com.hartwig.hmftools.common.variant.structural.annotation.GeneAnnotation;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class SvVarData
{
    private final String mIdStr; // sourced from DB so could be converted to int

    // full set of DB fields
    private final StructuralVariantData mSVData;
    private String[] mArm;
    private SvBreakend[] mBreakend;
    private boolean[] mFragileSite;
    private String[] mLineElement;

    private String[] mAssemblyData;

    private SvCluster mCluster;
    private String mClusterReason;

    private SvBreakend[] mFoldbackBreakends; // either the 2 breakends for this SV, or another SV's breaked
    private int[] mFoldbackLength;
    private String[] mFoldbackInfo;

    private long mNearestSvDistance;
    private String mNearestSvRelation;

    private SvLinkedPair mLinkStart; // templated insertion formed from this breakend to another
    private SvLinkedPair mLinkEnd;

    private SvLinkedPair[] mDbLink; // deletion bridge formed from this breakend to another
    private List<String> mStartTIAssemblies;
    private List<String> mEndTIAssemblies;
    private String[] mAssemblyMatchType;
    private boolean mIsReplicatedSv;
    private SvVarData mReplicatedSv;
    private int mReplicatedCount;

    private List<GeneAnnotation> mGenesStart;
    private List<GeneAnnotation> mGenesEnd;

    private String[] mDriverGene;

    private double[] mReplicationOrigin;

    // copy number related data
    private double[] mCopyNumber; // cached from SV data but modifiable
    private double[] mCopyNumberChange;
    private double mPloidy;

    private boolean mHasCalcPloidy;
    private double mPloidyMin;
    private double mPloidyMax;
    private SvCNData mCnDataPrevStart; // segment leading to the start position
    private SvCNData mCnDataPostStart; // segment starting with the position position
    private SvCNData mCnDataPrevEnd;
    private SvCNData mCnDataPostEnd;

    public static final String NONE_SEGMENT_INFERRED = "INFERRED";
    public static final String INF_SV_TYPE = "INF";

    public static String ASSEMBLY_TYPE_TI = "asm";
    public static String ASSEMBLY_TYPE_EQV = "eqv";

    public static String RELATION_TYPE_NEIGHBOUR = "NHBR";
    public static String RELATION_TYPE_OVERLAP = "OVRL";

    // iterators for start and end data
    public static int SE_START = 0;
    public static int SE_END = 1;
    public static int SE_PAIR = 2;

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
        mArm = new String[SE_PAIR];
        mFragileSite = new boolean[SE_PAIR];
        mLineElement = new String[] {NO_LINE_ELEMENT, NO_LINE_ELEMENT};

        mNearestSvDistance = -1;
        mNearestSvRelation = "";

        mIsReplicatedSv = false;
        mReplicatedSv = null;
        mReplicatedCount = 0;

        mClusterReason = "";
        mCluster = null;

        mLinkStart = null;
        mLinkEnd = null;
        mDbLink = new SvLinkedPair[SE_PAIR];

        mFoldbackBreakends = new SvBreakend[SE_PAIR];
        mFoldbackLength = new int[] {-1, -1};
        mFoldbackInfo = new String[] {"", ""};

        mGenesStart = Lists.newArrayList();
        mGenesEnd = Lists.newArrayList();

        mDriverGene = new String[] {"", ""};

        mReplicationOrigin = new double[SE_PAIR];

        mCopyNumber[SE_START] = mSVData.adjustedStartCopyNumber();
        mCopyNumber[SE_END] = mSVData.adjustedEndCopyNumber();
        mCopyNumberChange[SE_START] = mSVData.adjustedStartCopyNumberChange();
        mCopyNumberChange[SE_END] = mSVData.adjustedEndCopyNumberChange();
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
        mArm[SE_START] = other.arm(true);
        mArm[SE_END] = other.arm(false);

        mBreakend[SE_START] = new SvBreakend(this, true);
        mBreakend[SE_START].setChrPosIndex(other.getBreakend(true).getChrPosIndex());

        if(!isNullBreakend())
        {
            mBreakend[SE_END] = new SvBreakend(this, false);
            mBreakend[SE_END].setChrPosIndex(other.getBreakend(false).getChrPosIndex());
        }

        mFragileSite[SE_START] = other.isFragileSite(true);
        mFragileSite[SE_END] = other.isFragileSite(false);
        mLineElement[SE_START] = other.getLineElement(true);
        mLineElement[SE_END] = other.getLineElement(false);
        mNearestSvDistance = other.getNearestSvDistance();
        mNearestSvRelation = other.getNearestSvRelation();

        mAssemblyData[SE_START] = other.getAssemblyData(true);
        mAssemblyData[SE_END] = other.getAssemblyData(false);
        setAssemblyData(true);

        mAssemblyMatchType[SE_START] = other.getAssemblyMatchType(true);
        mAssemblyMatchType[SE_END] = other.getAssemblyMatchType(false);
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
    public final double copyNumber(boolean isStart){ return mCopyNumber[seIndex(isStart)]; }
    public final StructuralVariantType type() { return mSVData.type(); }

    public SvBreakend getBreakend(boolean isStart) { return mBreakend[seIndex(isStart)]; }

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

    public final String arm(boolean isStart) { return mArm[seIndex(isStart)]; }

    public void setChromosomalArms(final String start, final String end)
    {
        mArm[SE_START] = start;
        mArm[SE_END] = end;

        mBreakend[SE_START] = new SvBreakend(this, true);

        if(!isNullBreakend())
            mBreakend[SE_END] = new SvBreakend(this, false);
    }

    public final SvCluster getCluster() { return mCluster; }
    public void setCluster(final SvCluster cluster)
    {
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

    public void clearClusterReason() { mClusterReason = ""; }
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
        if(mDbLink[seIndex(isStart)] != null && mDbLink[seIndex(isStart)].length() == 0)
            return mPloidy;
        else
            return mCopyNumberChange[seIndex(isStart)];
    }

    public void setCopyNumberData(boolean isStart, double copyNumber, double copyNumberChange)
    {
        mCopyNumber[seIndex(isStart)] = copyNumber;
        mCopyNumberChange[seIndex(isStart)] = copyNumberChange;
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

    public void setFragileSites(boolean isFragleStart, boolean isFragileEnd)
    {
        mFragileSite[SE_START] = isFragleStart;
        mFragileSite[SE_END] = isFragileEnd;
    }

    public boolean isFragileSite(boolean useStart) { return mFragileSite[seIndex(useStart)]; }

    public void setLineElement(String type, boolean isStart)
    {
        int seIndex = seIndex(isStart);

        if(mLineElement[seIndex].contains(type))
            return;
        else if(mLineElement[seIndex].equals(NO_LINE_ELEMENT))
            mLineElement[seIndex] = type;
        else
            mLineElement[seIndex] = mLineElement[seIndex] + ";" + type;
    }

    public boolean isLineElement(boolean useStart)
    {
        return !mLineElement[seIndex(useStart)].equals(NO_LINE_ELEMENT);
    }

    public boolean inLineElement() { return isLineElement(true) || isLineElement(false); }
    public final String getLineElement(boolean useStart) { return mLineElement[seIndex(useStart)]; }

    public final SvLinkedPair getLinkedPair(boolean isStart) { return isStart ? mLinkStart : mLinkEnd; }

    public void setLinkedPair(final SvLinkedPair link, boolean isStart)
    {
        if(isStart)
            mLinkStart = link;
        else
            mLinkEnd = link;
    }

    public final SvLinkedPair getDBLink(boolean isStart) { return mDbLink[seIndex(isStart)]; }

    public void setDBLink(final SvLinkedPair link, boolean isStart)
    {
        mDbLink[seIndex(isStart)] = link;
    }

    public final String getFoldbackLink(boolean useStart)
    {
        if(mFoldbackBreakends[seIndex(useStart)] != null)
            return mFoldbackBreakends[seIndex(useStart)].getSV().id();
        else
            return "";
    }

    public final SvBreakend getFoldbackBreakend(boolean isStart) { return mFoldbackBreakends[seIndex(isStart)]; }
    public int getFoldbackLength(boolean isStart) { return mFoldbackLength[seIndex(isStart)]; }
    public final String getFoldbackInfo(boolean isStart) { return mFoldbackInfo[seIndex(isStart)]; }
    public boolean isFoldback() { return mFoldbackBreakends[SE_START] != null || mFoldbackBreakends[SE_END] != null; }
    public boolean isChainedFoldback()
    {
        if(mFoldbackBreakends[SE_END] != null && mFoldbackBreakends[SE_END] != mBreakend[SE_START])
            return true;
        else if(mFoldbackBreakends[SE_START] != null && mFoldbackBreakends[SE_START] != mBreakend[SE_END])
            return true;
        else
            return false;
    }

    public final SvBreakend getChainedFoldbackBreakend()
    {
        return mFoldbackBreakends[SE_START] != null ? mFoldbackBreakends[SE_START] : mFoldbackBreakends[SE_END];
    }

    public void setFoldbackLink(boolean isStart, final SvBreakend link, int length, String linkInfo)
    {
        mFoldbackBreakends[seIndex(isStart)] = link;
        mFoldbackLength[seIndex(isStart)] = length;
        mFoldbackInfo[seIndex(isStart)] = linkInfo;

        if(mCluster != null)
        {
            if (mFoldbackBreakends[SE_START] == null && mFoldbackBreakends[SE_END] == null)
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
        mFoldbackInfo[seIndex(isStart)] = info;
    }

    public final SvVarData getChainedFoldbackSv()
    {
        if(!isChainedFoldback())
            return this;

        return mFoldbackBreakends[SE_START] != null ? mFoldbackBreakends[SE_START].getSV() : mFoldbackBreakends[SE_END].getSV();
    }

    public final String typeStr()
    {
        if(isNoneSegment())
            return INF_SV_TYPE;
        else
            return type().toString();
    }

    public final boolean isLocal()
    {
        // means that both ends are within the same chromosomal arm
        return chromosome(true).equals(chromosome(false)) && mArm[SE_START].equals(mArm[SE_END]);
    }

    public final boolean isCrossArm()
    {
        return type() != SGL && !isLocal();
    }

    public final boolean isSimpleType()
    {
        return (type() == DEL || type() == DUP || type() == INS);
    }

    public static boolean isStart(int svIter) { return svIter == SE_START; }
    public static int seIndex(boolean isStart) { return isStart ? SE_START : SE_END; }

    public String getAssemblyData(boolean isStart) { return mAssemblyData[seIndex(isStart)]; }

    // unit testing only
    public void setAssemblyData(boolean isStart, final String data)
    {
        mAssemblyData[seIndex(isStart)] = data;
        setAssemblyData(true);
    }

    public final List<String> getTIAssemblies(boolean isStart)
    {
        return isStart ? mStartTIAssemblies : mEndTIAssemblies;
    }

    public boolean isEquivBreakend(boolean isStart)
    {
        return getAssemblyData(true).contains(ASSEMBLY_TYPE_EQV);
    }

    public boolean hasEquivBreakend()
    {
        return isEquivBreakend(true) || isEquivBreakend(false);
    }

    public final List<GeneAnnotation> getGenesList(boolean isStart) { return isStart ? mGenesStart : mGenesEnd; }
    public void setGenesList(final List<GeneAnnotation> genesList, boolean isStart)
    {
        if(isStart)
            mGenesStart.addAll(genesList);
        else
            mGenesEnd.addAll(genesList);
    }

    public final String getDriverGene(boolean isStart) { return mDriverGene[seIndex(isStart)]; }

    public void setDriveGene(final String geneInfo, boolean isStart)
    {
        mDriverGene[seIndex(isStart)] = geneInfo;
    }

    public final String getGeneInBreakend(boolean isStart)
    {
        final List<GeneAnnotation> genesList = getGenesList(isStart);

        String genesStr = "";
        for(final GeneAnnotation gene : genesList)
        {
            genesStr = appendStr(genesStr, gene.GeneName, ';');
        }

        return genesStr;
    }

    public final String getAssemblyMatchType(boolean isStart) { return mAssemblyMatchType[seIndex(isStart)]; }
    public boolean isAssemblyMatched(boolean isStart) { return getAssemblyMatchType(isStart).equals(ASSEMBLY_MATCH_MATCHED); }

    public void setAssemblyMatchType(String type, boolean isStart)
    {
        mAssemblyMatchType[seIndex(isStart)] = type;
    }

    private void setAssemblyData(boolean useExisting)
    {
        mStartTIAssemblies = Lists.newArrayList();
        mEndTIAssemblies = Lists.newArrayList();

        mAssemblyMatchType[SE_START] = ASSEMBLY_MATCH_NONE;
        mAssemblyMatchType[SE_END] = ASSEMBLY_MATCH_NONE;

        if(!useExisting)
        {
            mAssemblyData[SE_START] = "";
            mAssemblyData[SE_END] = "";

            if(!mSVData.startLinkedBy().isEmpty() && !mSVData.startLinkedBy().equals("."))
            {
                mAssemblyData[SE_START] = mSVData.startLinkedBy().replaceAll(",", ";");
            }

            if(!mSVData.endLinkedBy().isEmpty() && !mSVData.endLinkedBy().equals("."))
            {
                mAssemblyData[SE_END] = mSVData.endLinkedBy().replaceAll(",", ";");
            }
        }

        if(!mAssemblyData[SE_START].isEmpty())
        {
            String[] assemblyList = mAssemblyData[SE_START].split(";");

            for(int i = 0; i < assemblyList.length; ++i)
            {
                if(assemblyList[i].contains(ASSEMBLY_TYPE_TI))
                    mStartTIAssemblies.add(assemblyList[i]);
            }
        }

        if(!mAssemblyData[SE_END].isEmpty())
        {
            String[] assemblyList = mAssemblyData[SE_END].split(";");
            for(int i = 0; i < assemblyList.length; ++i)
            {
                if(assemblyList[i].contains(ASSEMBLY_TYPE_TI))
                    mEndTIAssemblies.add(assemblyList[i]);
            }
        }
    }

    public int getMinTemplatedLength(boolean isStart)
    {
        int assembleLength = isStart ? mSVData.startAnchoringSupportDistance() : mSVData.endAnchoringSupportDistance();
        return max(assembleLength - MIN_TEMPLATED_INSERTION_LENGTH, MIN_TEMPLATED_INSERTION_LENGTH);
    }

    public double getReplicationOrigin(boolean isStart) { return mReplicationOrigin[seIndex(isStart)]; }
    public void setReplicationOrigin(boolean isStart, double value)
    {
        mReplicationOrigin[seIndex(isStart)] = value;
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

    public boolean sglToCentromereOrTelomere()
    {
        if(!isNullBreakend() || isNoneSegment())
            return false;

        if(mSVData.insertSequenceRepeatClass().equals("Satellite/centr") || mSVData.insertSequenceRepeatClass().equals("Satellite/telo"))
            return true;

        if(mSVData.insertSequenceRepeatClass().equals("Simple_repeat"))
        {
            if (mSVData.insertSequenceRepeatType().equals("(CCCTAA)n") || mSVData.insertSequenceRepeatType().equals("(TTAGGG)n"))
                return true;
        }

        return false;
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

    public static boolean isSpecificSV(final SvVarData var)
    {
        if(var.id().equals(SPECIFIC_SV_ID))
            return true;

        return false;
    }

}
