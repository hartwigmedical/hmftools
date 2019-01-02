package com.hartwig.hmftools.svanalysis.types;

import static java.lang.Math.abs;
import static java.lang.Math.round;

import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.BND;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.DEL;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.DUP;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.INS;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.SGL;
import static com.hartwig.hmftools.svanalysis.annotators.FragileSiteAnnotator.NO_FS;
import static com.hartwig.hmftools.svanalysis.annotators.LineElementAnnotator.KNOWN_LINE_ELEMENT;
import static com.hartwig.hmftools.svanalysis.annotators.LineElementAnnotator.NO_LINE_ELEMENT;
import static com.hartwig.hmftools.svanalysis.types.SvLinkedPair.ASSEMBLY_MATCH_MATCHED;
import static com.hartwig.hmftools.svanalysis.types.SvLinkedPair.ASSEMBLY_MATCH_NONE;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalog;
import com.hartwig.hmftools.common.variant.structural.EnrichedStructuralVariant;
import com.hartwig.hmftools.common.variant.structural.ImmutableStructuralVariantData;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantData;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantType;
import com.hartwig.hmftools.common.variant.structural.annotation.GeneAnnotation;
import com.hartwig.hmftools.common.variant.structural.annotation.Transcript;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.ensembl.database.homo_sapiens_core.tables.Gene;

public class SvVarData
{
    private final String mIdStr; // sourced from DB so could be converted to int

    // full set of DB fields
    private final StructuralVariantData mSVData;
    private boolean mNoneSegment; // created from a NONE copy number segment
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

    private String mFoldbackLinkStart;
    private String mFoldbackLinkEnd;
    private int mFoldbackLenStart;
    private int mFoldbackLenEnd;
    private String mFoldbackLinkInfoStart;
    private String mFoldbackLinkInfoEnd;

    private long mNearestSvDistance;
    private String mNearestSvRelation;

    private SvLinkedPair mStartLink; // templated insertion formed from this breakend to another
    private SvLinkedPair mEndLink;

    private SvLinkedPair mStartDB; // deletion bridge formed from this breakend to another
    private SvLinkedPair mEndDB;
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

        setAssemblyData();
    }

    private void init()
    {
        mNoneSegment = false;
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

        mStartLink = null;
        mEndLink = null;
        mStartDB = null;
        mEndDB = null;

        mFoldbackLinkStart = "";
        mFoldbackLinkEnd = "";
        mFoldbackLenStart = -1;
        mFoldbackLenEnd = -1;
        mFoldbackLinkInfoStart = "";
        mFoldbackLinkInfoEnd = "";

        mConsecBEStart = "";
        mConsecBEEnd = "";

        mGenesStart = Lists.newArrayList();
        mGenesEnd = Lists.newArrayList();

        mDriverGeneStart = "";
        mDriverGeneEnd = "";
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

        mIdStr = other.getSvData().id() + "r";
        mSVData = other.getSvData();
        mBreakendStart = other.getBreakend(true);
        mBreakendEnd = other.getBreakend(false);
        mNoneSegment = other.isNoneSegment();
        mStartArm = other.arm(true);
        mEndArm = other.arm(false);
        mStartFragileSite = other.isFragileSite(true);
        mEndFragileSite = other.isFragileSite(false);
        mStartLineElement = other.getLineElement(true);
        mEndLineElement = other.getLineElement(false);
        mNearestSvDistance = other.getNearestSvDistance();
        mNearestSvRelation = other.getNearestSvRelation();
        setAssemblyData();
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
    public void setNoneSegment(boolean toggle) { mNoneSegment = toggle; }
    public boolean isNoneSegment() { return mNoneSegment; }

    // for convenience
    public final String chromosome(boolean isStart) { return isStart ? mSVData.startChromosome() : mSVData.endChromosome(); }
    public final long position(boolean isStart) { return isStart ? mSVData.startPosition() : mSVData.endPosition(); }
    public final byte orientation(boolean isStart){ return isStart ? mSVData.startOrientation() : mSVData.endOrientation(); }
    public final double copyNumber(boolean isStart){ return isStart ? mSVData.adjustedStartCopyNumber() : mSVData.adjustedEndCopyNumber(); }
    public final double copyNumberChange(boolean isStart){ return isStart ? mSVData.adjustedStartCopyNumberChange() : mSVData.adjustedEndCopyNumberChange(); }
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

        if(otherId.equals(mIdStr))
        {
            LOGGER.warn("SV({}) reason({}) setting to own ID", mIdStr, reason);
        }
    }

    public final String getClusterReason() { return mClusterReason; }

    public boolean isReplicatedSv() { return mIsReplicatedSv; }
    public final SvVarData getReplicatedSv() { return mReplicatedSv; }
    public int getReplicatedCount() { return mReplicatedCount; }
    public void setReplicatedCount(int count) { mReplicatedCount = count; }
    public final String origId() { return mReplicatedSv != null ? mReplicatedSv.id() : mIdStr; }
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

    public int getCopyNumberChange(boolean useStart)
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

    public final SvLinkedPair getLinkedPair(boolean isStart) { return isStart ? mStartLink : mEndLink; }

    public void setLinkedPair(final SvLinkedPair link, boolean isStart)
    {
        if(isStart)
            mStartLink = link;
        else
            mEndLink = link;
    }

    public final SvLinkedPair getDBLink(boolean isStart) { return isStart ? mStartDB : mEndDB; }

    public void setDBLink(final SvLinkedPair link, boolean isStart)
    {
        if(isStart)
            mStartDB = link;
        else
            mEndDB = link;
    }

    public final String getFoldbackLink(boolean useStart) { return useStart ? mFoldbackLinkStart : mFoldbackLinkEnd; }
    public int getFoldbackLen(boolean useStart) { return useStart ? mFoldbackLenStart : mFoldbackLenEnd; }
    public final String getFoldbackLinkInfo(boolean useStart) { return useStart ? mFoldbackLinkInfoStart : mFoldbackLinkInfoEnd; }
    public void setFoldbackLink(boolean isStart, String link, int length, String linkInfo)
    {
        if(isStart)
        {
            mFoldbackLinkStart = link;
            mFoldbackLenStart = length;
            mFoldbackLinkInfoStart = linkInfo;
        }
        else
        {
            mFoldbackLinkEnd = link;
            mFoldbackLenEnd = length;
            mFoldbackLinkInfoEnd = linkInfo;
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
    public void setAssemblyData(boolean useStart, final String data)
    {
        if (useStart)
            mAssemblyStartData = data;
        else
            mAssemblyEndData = data;
    }

    public List<String> getTempInsertionAssemblies(boolean useStart) { return useStart ? mStartTempInsertionAssemblies : mEndTempInsertionAssemblies; }

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

        return longestTrans != null ? longestTrans.parent().geneName() : "";
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

    private void setAssemblyData()
    {
        mStartTempInsertionAssemblies = Lists.newArrayList();
        mEndTempInsertionAssemblies = Lists.newArrayList();
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

    // private static String SPECIFIC_VAR_ID = "4504083";
    private static String SPECIFIC_VAR_ID = "";

    public static boolean isSpecificSV(final String id)
    {
        if(id.equals(SPECIFIC_VAR_ID))
            return true;

        return false;
    }

}
