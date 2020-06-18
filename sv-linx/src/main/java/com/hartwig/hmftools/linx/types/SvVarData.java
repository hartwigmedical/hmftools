package com.hartwig.hmftools.linx.types;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.round;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.stripChromosome;
import static com.hartwig.hmftools.common.utils.Strings.appendStr;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_PAIR;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.seIndex;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.BND;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.DEL;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.DUP;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.INF;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.INS;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.SGL;
import static com.hartwig.hmftools.linx.LinxConfig.LNX_LOGGER;
import static com.hartwig.hmftools.linx.LinxOutput.SUBSET_DELIM;
import static com.hartwig.hmftools.linx.LinxOutput.SUBSET_SPLIT;
import static com.hartwig.hmftools.linx.annotators.LineElementAnnotator.NO_LINE_ELEMENT;
import static com.hartwig.hmftools.linx.types.SglMapping.convertFromInsertSequenceAlignments;
import static com.hartwig.hmftools.linx.types.SvConstants.MIN_TEMPLATED_INSERTION_LENGTH;

import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantData;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantType;
import com.hartwig.hmftools.common.fusion.GeneAnnotation;
import com.hartwig.hmftools.linx.analysis.ClusteringReason;
import com.hartwig.hmftools.linx.cn.SvCNData;

public class SvVarData
{
    // full set of DB fields
    private final StructuralVariantData mSVData;
    private final String[] mChr; // stripped of 'chr' for logging
    private final ChromosomeArm[] mArm;
    private final SvBreakend[] mBreakend;
    private final boolean[] mFragileSite;
    private final String[] mLineElement;

    private final String[] mAssemblyData;

    private SvCluster mCluster;
    private String mClusterReason;

    private final SvBreakend[] mFoldbackBreakends; // either the 2 breakends for this SV, or another SV's breaked
    private final int[] mFoldbackLength;
    private final String[] mFoldbackInfo;

    private int mNearestSvDistance;
    private String mNearestSvRelation;

    private final List<SvLinkedPair>[] mTiLinks; // start and end lists of inferred or assembled TIs

    private final SvLinkedPair[] mDbLink; // deletion bridge formed from this breakend to another
    private final List<String>[] mTIAssemblies;

    private final List<GeneAnnotation>[] mGenes;

    private final double[] mReplicationOrigin;

    // copy number related data
    private final double[] mCopyNumber; // cached from SV data but modifiable
    private final double[] mCopyNumberChange;
    private double mJcn;

    private boolean mHasCalcJcn;
    private double mJcnMin;
    private double mJcnMax;
    private SvCNData mCnDataPrevStart; // segment leading to the start position
    private SvCNData mCnDataPostStart; // segment starting with the position position
    private SvCNData mCnDataPrevEnd;
    private SvCNData mCnDataPostEnd;

    private final List<SglMapping> mSglMappings;

    public static final String NONE_SEGMENT_INFERRED = "INFERRED";
    public static final String INF_SV_TYPE = "INF";

    public static String ASSEMBLY_TYPE_TI = "asm";
    public static String ASSEMBLY_TYPE_EQV = "eqv";

    public static String RELATION_TYPE_NEIGHBOUR = "NHBR";
    public static String RELATION_TYPE_OVERLAP = "OVRL";

    public SvVarData(final StructuralVariantData svData)
    {
        mSVData = svData;

        mArm = new ChromosomeArm[SE_PAIR];
        mChr = new String[] { stripChromosome(chromosome(true)), stripChromosome(chromosome(false)) };

        mFragileSite = new boolean[SE_PAIR];
        mLineElement = new String[] {NO_LINE_ELEMENT, NO_LINE_ELEMENT};
        mBreakend = new SvBreakend[SE_PAIR];

        mNearestSvDistance = -1;
        mNearestSvRelation = "";

        mClusterReason = "";
        mCluster = null;

        mDbLink = new SvLinkedPair[SE_PAIR];
        mTiLinks = new List[] { Lists.newArrayList(), Lists.newArrayList() };

        mFoldbackBreakends = new SvBreakend[SE_PAIR];
        mFoldbackLength = new int[] {-1, -1};
        mFoldbackInfo = new String[] {"", ""};

        mGenes = new List[] { Lists.newArrayList(), Lists.newArrayList() };

        mReplicationOrigin = new double[SE_PAIR];

        mAssemblyData = new String[SE_PAIR];
        mTIAssemblies = new List[] { Lists.newArrayList(), Lists.newArrayList() };

        mCopyNumber = new double[] { mSVData.adjustedStartCopyNumber(),  mSVData.adjustedEndCopyNumber() };
        mCopyNumberChange = new double[] {mSVData.adjustedStartCopyNumberChange(), mSVData.adjustedEndCopyNumberChange() };
        mJcn = mSVData.junctionCopyNumber();

        mHasCalcJcn = false;
        mJcnMin = 0;
        mJcnMax = 0;
        mCnDataPrevStart = null;
        mCnDataPostStart = null;
        mCnDataPrevEnd = null;
        mCnDataPostEnd = null;

        setAssemblyData(false);

        if(isSglBreakend())
        {
            mSglMappings = Lists.newArrayList();
            convertFromInsertSequenceAlignments(mSglMappings, svData.insertSequenceAlignments());
        }
        else
        {
            mSglMappings = null;
        }
    }

    public final int id() { return mSVData.id(); }
    public final String idStr() { return String.valueOf(mSVData.id()); }

    public String toString() { return posId() + " " + typeStr(); }

    public final StructuralVariantData getSvData() { return mSVData; }

    // for convenience
    public final String chromosome(boolean isStart) { return isStart ? mSVData.startChromosome() : mSVData.endChromosome(); }
    public final int position(boolean isStart) { return isStart ? mSVData.startPosition() : mSVData.endPosition(); }
    public final byte orientation(boolean isStart){ return isStart ? mSVData.startOrientation() : mSVData.endOrientation(); }
    public final double copyNumber(boolean isStart){ return mCopyNumber[seIndex(isStart)]; }
    public final StructuralVariantType type() { return mSVData.type(); }

    public SvBreakend getBreakend(int seIndex) { return mBreakend[seIndex]; }
    public SvBreakend getBreakend(boolean isStart) { return mBreakend[seIndex(isStart)]; }

    public boolean isSglBreakend() { return type() == SGL || type() == INF; }
    public boolean isInferredSgl() { return type() == INF; }

    public final String posId()
    {
        if(isSglBreakend())
        {
            return String.format("id(%s) pos(%s:%d:%d)",
                    id(), mChr[SE_START], orientation(true), position(true));
        }
        else
        {
            return String.format("id(%s) pos(%s:%d:%d -> %s:%d:%d)",
                    id(), mChr[SE_START], orientation(true), position(true),
                    mChr[SE_END], orientation(false), position(false));
        }
    }

    public final String posId(boolean useStart)
    {
        return String.format("%s: %s %s:%d:%d",
                id(), useStart ? "start" :"end", mChr[seIndex(useStart)], orientation(useStart), position(useStart));
    }

    public final ChromosomeArm arm(boolean isStart) { return mArm[seIndex(isStart)]; }
    public final String chrShort(boolean isStart) { return mChr[seIndex(isStart)]; }

    public void setChromosomalArms(final ChromosomeArm start, final ChromosomeArm end)
    {
        mArm[SE_START] = start;
        mArm[SE_END] = end;

        mBreakend[SE_START] = new SvBreakend(this, true);

        if(!isSglBreakend())
            mBreakend[SE_END] = new SvBreakend(this, false);
    }

    public final SvCluster getCluster() { return mCluster; }
    public void setCluster(final SvCluster cluster) { mCluster = cluster; }

    public final int length()
    {
        if(type() == BND || isSglBreakend())
            return 0;

        return abs(position(false) - position(true));
    }

    public static final String CR_DELIM = "-";

    public void addClusterReason(final ClusteringReason reason, final int otherId)
    {
        if(mClusterReason.contains(reason.toString()))
            return;

        mClusterReason = appendStr(mClusterReason, reason.toString(), SUBSET_DELIM);

        if(otherId > -1)
            mClusterReason += CR_DELIM + otherId;

        if(otherId == id())
        {
            LNX_LOGGER.warn("SV({}) reason({}) setting to own ID", id(), reason);
        }
    }

    public final String getClusterReason() { return mClusterReason; }
    public boolean hasClusterReason(ClusteringReason reason) { return mClusterReason.contains(reason.toString()); }

    public double jcn() { return mJcn; }

    public double copyNumberChange(boolean isStart)
    {
        // TEMP: precise DBs cause incorrect copy number change, so in this case use JCN
        if(mDbLink[seIndex(isStart)] != null && mDbLink[seIndex(isStart)].length() == 0)
            return mJcn;
        else
            return mCopyNumberChange[seIndex(isStart)];
    }

    public void setCopyNumberData(boolean isStart, double copyNumber, double copyNumberChange)
    {
        mCopyNumber[seIndex(isStart)] = copyNumber;
        mCopyNumberChange[seIndex(isStart)] = copyNumberChange;
    }

    public double getRoundedJcn(boolean enforceClonal)
    {
        double roundedJcn = round(jcn());
        return enforceClonal ? max(roundedJcn, 1) : roundedJcn;
    }

    public int getMaxAssembledBreakend()
    {
        return max(getAssembledLinkedPairs(true).size(), getAssembledLinkedPairs(false).size());
    }

    public int getImpliedJcn()
    {
        return max(getMaxAssembledBreakend(), (int) getRoundedJcn(true));
    }

    public int getNearestSvDistance() { return mNearestSvDistance; }
    public void setNearestSvDistance(int distance) { mNearestSvDistance = distance; }
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
            mLineElement[seIndex] = mLineElement[seIndex] + SUBSET_SPLIT + type;
    }

    public boolean isLineElement(boolean useStart)
    {
        return !mLineElement[seIndex(useStart)].equals(NO_LINE_ELEMENT);
    }

    public boolean inLineElement() { return isLineElement(true) || isLineElement(false); }
    public final String getLineElement(boolean useStart) { return mLineElement[seIndex(useStart)]; }

    public final List<SvLinkedPair> getLinkedPairs(boolean isStart)
    {
        return mTiLinks[seIndex(isStart)];
    }

    public final List<SvLinkedPair> getAssembledLinkedPairs(boolean isStart)
    {
        return mTiLinks[seIndex(isStart)].stream().filter(SvLinkedPair::isAssembled).collect(Collectors.toList());
    }

    public final SvLinkedPair getLinkedPair(boolean isStart)
    {
        return mTiLinks[seIndex(isStart)].isEmpty() ? null : mTiLinks[seIndex(isStart)].get(0);
    }

    public void addLinkedPair(final SvLinkedPair link, boolean isStart)
    {
        // add in order from shortest to longest
        final List<SvLinkedPair> links = mTiLinks[seIndex(isStart)];

        int index = 0;
        while(index < links.size())
        {
            final SvLinkedPair otherPair = links.get(index);

            if(otherPair.matches(link))
                return;

            if(link.length() < links.get(index).length())
                break;

            ++index;
        }

        links.add(index, link);
    }

    public final SvLinkedPair getDBLink(boolean isStart) { return mDbLink[seIndex(isStart)]; }
    public void setDBLink(final SvLinkedPair link, boolean isStart)
    {
        mDbLink[seIndex(isStart)] = link;
    }

    public final int getFoldbackId(boolean useStart)
    {
        if(mFoldbackBreakends[seIndex(useStart)] != null)
            return mFoldbackBreakends[seIndex(useStart)].getSV().id();
        else
            return -1;
    }

    public final SvBreakend getFoldbackBreakend(boolean isStart) { return mFoldbackBreakends[seIndex(isStart)]; }
    public int getFoldbackLength(boolean isStart) { return mFoldbackLength[seIndex(isStart)]; }
    public final String getFoldbackInfo(boolean isStart) { return mFoldbackInfo[seIndex(isStart)]; }

    public boolean isFoldback() { return mFoldbackBreakends[SE_START] != null || mFoldbackBreakends[SE_END] != null; }

    public boolean isSameSvFoldback() { return mFoldbackBreakends[SE_END] != null && mFoldbackBreakends[SE_END] == mBreakend[SE_START]; }

    public boolean isSingleBreakendFoldback()
    {
        if(mFoldbackBreakends[SE_START] != null && mFoldbackBreakends[SE_START] == mBreakend[SE_START])
            return true;
        else if(mFoldbackBreakends[SE_END] != null && mFoldbackBreakends[SE_END] == mBreakend[SE_END])
            return true;

        return false;
    }

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
        if(isInferredSgl())
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
        return !isSglBreakend() && !isLocal();
    }

    public final boolean isSimpleType()
    {
        return (type() == DEL || type() == DUP || type() == INS);
    }

   public String getAssemblyData(boolean isStart) { return mAssemblyData[seIndex(isStart)]; }

    // unit testing only
    public void setAssemblyData(boolean isStart, final String data)
    {
        mAssemblyData[seIndex(isStart)] = data;
        setAssemblyData(true);
    }

    public final List<String> getTIAssemblies(boolean isStart)
    {
        return mTIAssemblies[seIndex(isStart)];
    }

    public boolean isEquivBreakend()
    {
        return getAssemblyData(true).contains(ASSEMBLY_TYPE_EQV);
    }

    public final List<SglMapping> getSglMappings() { return mSglMappings; }

    public final List<GeneAnnotation> getGenesList(boolean isStart) { return mGenes[seIndex(isStart)]; }

    private static final int PRE_TRANSCRIPT_DISTANCE = 10000;

    public final String getGeneInBreakend(boolean isStart, boolean includeId)
    {
        // create a list of any genes which this breakend touches, but exclude the upstream distance used for fusions
        final List<GeneAnnotation> genesList = getGenesList(isStart).stream()
                .filter(x -> x.breakendWithinGene(PRE_TRANSCRIPT_DISTANCE))
                .collect(Collectors.toList());

        String genesStr = "";
        for(final GeneAnnotation gene : genesList)
        {
            String geneStr = includeId ? gene.StableId + ":" + gene.GeneName : gene.GeneName;
            genesStr = appendStr(genesStr, geneStr, ';');
        }

        return genesStr;
    }

    public boolean hasAssemblyLink(boolean isStart)
    {
        return mTiLinks[seIndex(isStart)].stream().anyMatch(SvLinkedPair::isAssembled);
    }

    private void setAssemblyData(boolean useExisting)
    {
        for(int se = SE_START; se <= SE_END; ++se)
        {
            List<String> tiAssemblies = mTIAssemblies[se];

            if(!useExisting)
            {
                mAssemblyData[se] = "";

                final String linkedByData = se == SE_START ? mSVData.startLinkedBy() : mSVData.endLinkedBy();

                if(!linkedByData.isEmpty() && !linkedByData.equals("."))
                {
                    mAssemblyData[se] = linkedByData.replaceAll(",", ";");
                }
            }

            if(!mAssemblyData[se].isEmpty())
            {
                String[] assemblyList = mAssemblyData[se].split(";");

                for(int i = 0; i < assemblyList.length; ++i)
                {
                    if(assemblyList[i].contains(ASSEMBLY_TYPE_TI))
                        tiAssemblies.add(assemblyList[i]);
                }
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

    public void setJcnRecalcData(double minJcn, double maxJcn)
    {
        if(maxJcn == 0)
        {
            // suggests an error in the calculation
            mJcnMax = mJcn;
        }
        else
        {
            mJcnMin = minJcn;
            mJcnMax = maxJcn;
            mJcn = (mJcnMin + mJcnMax) * 0.5;
        }

        mHasCalcJcn = true;
    }

    public boolean hasCalculatedJcn() { return mHasCalcJcn; }
    public double jcnMax() { return mHasCalcJcn ? mJcnMax : mJcn; }
    public double jcnMin() { return mHasCalcJcn ? mJcnMin : mJcn; }
    public double jcnUncertainty() { return mHasCalcJcn ? (mJcn - mJcnMin) * 0.5 : 0; }

    public double calcVaf(boolean isStart)
    {
        return copyNumberChange(isStart) / copyNumber(isStart);
    }

    public final SvCNData getCopyNumberData(boolean isStart, boolean isPrevious)
    {
        return isStart ? (isPrevious ? mCnDataPrevStart : mCnDataPostStart) : (isPrevious ? mCnDataPrevEnd : mCnDataPostEnd);
    }

    public boolean sglToCentromereOrTelomere()
    {
        if(!isSglBreakend() || isInferredSgl())
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

    public boolean sglToSatelliteRepeats()
    {
        if(!isSglBreakend() || isInferredSgl())
            return false;

        if(mSVData.insertSequenceRepeatClass().equals("Satellite/centr"))
            return true;

        final String repeatType = mSVData.insertSequenceRepeatType();

        return repeatType.equals("(CATTC)n") || repeatType.equals("(GAATG)n") || repeatType.equals("HSATII") || repeatType.equals("SAR");
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
}
