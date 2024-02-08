package com.hartwig.hmftools.linx.types;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.round;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.ITEM_DELIM;
import static com.hartwig.hmftools.common.utils.Strings.appendStr;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_PAIR;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.seIndex;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.BND;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.DEL;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.DUP;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.INF;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.INS;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.SGL;
import static com.hartwig.hmftools.linx.LinxConfig.LNX_LOGGER;
import static com.hartwig.hmftools.linx.LinxOutput.ITEM_DELIM_CHR;
import static com.hartwig.hmftools.linx.types.SglMapping.convertFromInsertSequenceAlignments;

import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.linx.gene.BreakendGeneData;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeFunctions;
import com.hartwig.hmftools.common.purple.ChromosomeArm;
import com.hartwig.hmftools.common.utils.sv.StartEndPair;
import com.hartwig.hmftools.common.sv.StructuralVariantData;
import com.hartwig.hmftools.common.sv.StructuralVariantType;
import com.hartwig.hmftools.linx.analysis.ClusteringReason;
import com.hartwig.hmftools.linx.annotators.LineElementType;
import com.hartwig.hmftools.linx.cn.SvCNData;

public class SvVarData
{
    // full set of DB fields
    private final StructuralVariantData mSVData;
    private final String[] mChr; // stripped of 'chr' for logging
    private final ChromosomeArm[] mArm;
    private final SvBreakend[] mBreakend;
    private final boolean[] mFragileSite;
    private final StartEndPair<Set<LineElementType>> mLineElements;

    private final String[] mAssemblyData;

    private SvCluster mCluster;
    private String mClusterReason;

    private final SvBreakend[] mFoldbackBreakends; // either the 2 breakends for this SV, or another SV's breakend
    private final int[] mFoldbackLength;
    private final String[] mFoldbackInfo;

    private int mNearestSvDistance;
    private String mNearestSvRelation;

    private final StartEndPair<List<LinkedPair>> mTiLinks; // start and end lists of inferred or assembled TIs

    private final DbPair[] mDbLink; // deletion bridge formed from this breakend to another
    private final StartEndPair<List<String>> mTIAssemblies;

    private final StartEndPair<List<BreakendGeneData>> mGenes;

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
    private SvVarData[] mLinkedSVs;
    private List<String> mAnnotationList;

    public static final String NONE_SEGMENT_INFERRED = "INFERRED";
    public static final String INF_SV_TYPE = "INF";
    public static final String SGL_CENTRO_SATELLITE = "Satellite/centr";
    public static final String SGL_TELO_SATELLITE = "Satellite/telo";

    public static final String ASSEMBLY_TYPE_TI = "asm";
    public static final String TRANSITIVE_TYPE_TI = "trs";
    public static final String ASSEMBLY_TYPE_EQV = "eqv";

    public static final String RELATION_TYPE_NEIGHBOUR = "NHBR";
    public static final String RELATION_TYPE_OVERLAP = "OVRL";

    public SvVarData(final StructuralVariantData svData)
    {
        mSVData = svData;

        mArm = new ChromosomeArm[SE_PAIR];
        mChr = new String[] { RefGenomeFunctions.stripChrPrefix(chromosome(true)), RefGenomeFunctions.stripChrPrefix(chromosome(false)) };

        mFragileSite = new boolean[SE_PAIR];
        mLineElements = new StartEndPair<>(Sets.newHashSet(), Sets.newHashSet());
        mBreakend = new SvBreakend[SE_PAIR];

        mNearestSvDistance = -1;
        mNearestSvRelation = "";

        mClusterReason = "";
        mCluster = null;

        mDbLink = new DbPair[SE_PAIR];
        mTiLinks = new StartEndPair<>(Lists.newArrayList(), Lists.newArrayList());

        mFoldbackBreakends = new SvBreakend[SE_PAIR];
        mFoldbackLength = new int[] {-1, -1};
        mFoldbackInfo = new String[] {"", ""};

        mGenes = new StartEndPair<>(Lists.newArrayList(), Lists.newArrayList());

        mAssemblyData = new String[SE_PAIR];
        mTIAssemblies = new StartEndPair<>(Lists.newArrayList(), Lists.newArrayList());

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
            convertFromInsertSequenceAlignments(mSglMappings, svData.insertSequenceAlignments(), orientation(true));
        }
        else
        {
            mSglMappings = null;
        }

        mLinkedSVs = null;
        mAnnotationList = null;
    }

    public int id() { return mSVData.id(); }
    public String idStr() { return String.valueOf(mSVData.id()); }

    public String toString() { return posId() + " " + typeStr(); }

    public StructuralVariantData getSvData() { return mSVData; }

    // for convenience
    public String chromosome(boolean isStart) { return isStart ? mSVData.startChromosome() : mSVData.endChromosome(); }
    public int position(boolean isStart) { return isStart ? mSVData.startPosition() : mSVData.endPosition(); }
    public byte orientation(boolean isStart){ return isStart ? mSVData.startOrientation() : mSVData.endOrientation(); }
    public double copyNumber(boolean isStart){ return mCopyNumber[seIndex(isStart)]; }
    public StructuralVariantType type() { return mSVData.type(); }

    public SvBreakend getBreakend(int seIndex) { return mBreakend[seIndex]; }
    public SvBreakend getBreakend(boolean isStart) { return mBreakend[seIndex(isStart)]; }
    public void setSglEndBreakend(final SvBreakend breakend) { mBreakend[SE_END] = breakend; }

    public boolean isSglBreakend() { return type() == SGL || type() == INF; }
    public boolean isInferredSgl() { return type() == INF; }

    public String posId()
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

    public String posId(boolean useStart)
    {
        return String.format("%s: %s %s:%d:%d",
                id(), useStart ? "start" :"end", mChr[seIndex(useStart)], orientation(useStart), position(useStart));
    }

    public ChromosomeArm arm(boolean isStart) { return mArm[seIndex(isStart)]; }
    public String chrShort(boolean isStart) { return mChr[seIndex(isStart)]; }

    public void setChromosomalArms(final ChromosomeArm start, final ChromosomeArm end)
    {
        mArm[SE_START] = start;
        mArm[SE_END] = end;

        mBreakend[SE_START] = new SvBreakend(this, true);

        if(!isSglBreakend())
            mBreakend[SE_END] = new SvBreakend(this, false);
    }

    public SvCluster getCluster() { return mCluster; }
    public void setCluster(final SvCluster cluster) { mCluster = cluster; }

    public int length()
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

        mClusterReason = appendStr(mClusterReason, reason.toString(), ITEM_DELIM_CHR);

        if(otherId > -1)
            mClusterReason += CR_DELIM + otherId;

        if(otherId == id())
        {
            LNX_LOGGER.warn("SV({}) reason({}) setting to own ID", id(), reason);
        }
    }

    public String getClusterReason() { return mClusterReason; }
    public boolean hasClusterReason(ClusteringReason reason) { return mClusterReason.contains(reason.toString()); }

    public double jcn() { return mJcn; }

    public double copyNumberChange(boolean isStart)
    {
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

    public void setFragileSites(boolean isFragileStart, boolean isFragileEnd)
    {
        mFragileSite[SE_START] = isFragileStart;
        mFragileSite[SE_END] = isFragileEnd;
    }

    public boolean isFragileSite(boolean useStart) { return mFragileSite[seIndex(useStart)]; }

    public void addLineElement(LineElementType type, boolean isStart) { mLineElements.get(isStart).add(type); }
    public boolean isLineElement(boolean useStart) { return !mLineElements.get(useStart).isEmpty(); }
    public boolean inLineElement() { return isLineElement(true) || isLineElement(false); }
    public boolean hasLineElement(LineElementType type, boolean isStart) { return mLineElements.get(isStart).contains(type); }
    public Set<LineElementType> getLineElement(boolean useStart) { return mLineElements.get(useStart); }

    public List<LinkedPair> getLinkedPairs(boolean isStart) { return mTiLinks.get(isStart); }

    public List<LinkedPair> getAssembledLinkedPairs(boolean isStart)
    {
        return mTiLinks.get(isStart).stream().filter(LinkedPair::isAssembled).collect(Collectors.toList());
    }

    public LinkedPair getLinkedPair(boolean isStart)
    {
        return mTiLinks.get(isStart).isEmpty() ? null : mTiLinks.get(isStart).get(0);
    }

    public void addLinkedPair(final LinkedPair link, boolean isStart)
    {
        // add in order from shortest to longest
        final List<LinkedPair> links = mTiLinks.get(isStart);

        int index = 0;
        while(index < links.size())
        {
            final LinkedPair otherPair = links.get(index);

            if(otherPair.matches(link))
                return;

            if(link.positionDistance() < links.get(index).positionDistance())
                break;

            ++index;
        }

        links.add(index, link);
    }

    public void clearClusteringData()
    {
        mClusterReason = "";
        mTiLinks.Start.clear();
        mTiLinks.End.clear();
    }

    public DbPair getDBLink(boolean isStart) { return mDbLink[seIndex(isStart)]; }
    public void setDBLink(final DbPair link, boolean isStart)
    {
        mDbLink[seIndex(isStart)] = link;
    }

    public int getFoldbackId(boolean useStart)
    {
        if(mFoldbackBreakends[seIndex(useStart)] != null)
            return mFoldbackBreakends[seIndex(useStart)].getSV().id();
        else
            return -1;
    }

    public SvBreakend getFoldbackBreakend(boolean isStart) { return mFoldbackBreakends[seIndex(isStart)]; }
    public int getFoldbackLength(boolean isStart) { return mFoldbackLength[seIndex(isStart)]; }
    public String getFoldbackInfo(boolean isStart) { return mFoldbackInfo[seIndex(isStart)]; }

    public boolean isFoldback() { return mFoldbackBreakends[SE_START] != null || mFoldbackBreakends[SE_END] != null; }

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

    public void setFoldbackLink(boolean isStart, final SvBreakend link, int length, String linkInfo)
    {
        mFoldbackBreakends[seIndex(isStart)] = link;
        mFoldbackLength[seIndex(isStart)] = length;
        mFoldbackInfo[seIndex(isStart)] = linkInfo;

        if(mCluster != null)
        {
            if(mFoldbackBreakends[SE_START] == null && mFoldbackBreakends[SE_END] == null)
            {
                mCluster.deregisterFoldback(this);
            }
            else
            {
                mCluster.registerFoldback(this);
            }
        }
    }

    public String typeStr()
    {
        if(isInferredSgl())
            return INF_SV_TYPE;
        else
            return type().toString();
    }

    public boolean isLocal()
    {
        // means that both ends are within the same chromosomal arm
        return chromosome(true).equals(chromosome(false)) && mArm[SE_START].equals(mArm[SE_END]);
    }

    public boolean isCrossArm()
    {
        return !isSglBreakend() && !isLocal();
    }

    public boolean isSimpleType()
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

    public List<String> getTIAssemblies(boolean isStart) { return mTIAssemblies.get(isStart); }

    public boolean isEquivBreakend()
    {
        return getAssemblyData(true).contains(ASSEMBLY_TYPE_EQV);
    }

    public void setLinkedSVs(final SvVarData var1, final SvVarData var2)
    {
        mLinkedSVs = new SvVarData[] {var1, var2};
    }

    public SvVarData[] getLinkedSVs() { return mLinkedSVs; }

    public List<BreakendGeneData> getGenesList(boolean isStart) { return mGenes.get(isStart); }

    public List<SglMapping> getSglMappings() { return mSglMappings; }

    private static final int PRE_TRANSCRIPT_DISTANCE = 10000;

    public String getGeneInBreakend(boolean isStart, boolean includeId)
    {
        // create a list of any genes which this breakend touches, but exclude the upstream distance used for fusions
        final List<BreakendGeneData> genesList = getGenesList(isStart).stream()
                .filter(x -> x.breakendWithinGene(PRE_TRANSCRIPT_DISTANCE))
                .collect(Collectors.toList());

        String genesStr = "";
        for(final BreakendGeneData gene : genesList)
        {
            String geneStr = includeId ? gene.geneId() + ":" + gene.geneName() : gene.geneName();
            genesStr = appendStr(genesStr, geneStr, ITEM_DELIM_CHR);
        }

        return genesStr;
    }

    public boolean hasAssemblyLink(boolean isStart)
    {
        return mTiLinks.get(isStart).stream().anyMatch(LinkedPair::isAssembled);
    }

    private void setAssemblyData(boolean useExisting)
    {
        for(int se = SE_START; se <= SE_END; ++se)
        {
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

                for(String assembly : assemblyList)
                {
                    if(assembly.contains(ASSEMBLY_TYPE_TI) || assembly.contains(TRANSITIVE_TYPE_TI))
                        mTIAssemblies.get(se).add(assembly);
                }
            }
        }
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

    public SvCNData getCopyNumberData(boolean isStart, boolean isPrevious)
    {
        return isStart ? (isPrevious ? mCnDataPrevStart : mCnDataPostStart) : (isPrevious ? mCnDataPrevEnd : mCnDataPostEnd);
    }

    public boolean sglToCentromereOrTelomere()
    {
        if(!isSglBreakend() || isInferredSgl())
            return false;

        if(mSVData.insertSequenceRepeatClass().equals(SGL_CENTRO_SATELLITE) || mSVData.insertSequenceRepeatClass().equals(SGL_TELO_SATELLITE))
            return true;

        if(mSVData.insertSequenceRepeatClass().equals("Simple_repeat"))
        {
            if(mSVData.insertSequenceRepeatType().equals("(CCCTAA)n") || mSVData.insertSequenceRepeatType().equals("(TTAGGG)n"))
                return true;
        }

        return false;
    }

    public boolean sglToSatelliteRepeats()
    {
        if(type() != SGL)
            return false;

        if(mSVData.insertSequenceRepeatClass().equals(SGL_CENTRO_SATELLITE))
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

    public List<String> getAnnotationList() { return mAnnotationList; }

    public void addAnnotation(final String annotation)
    {
        if(mAnnotationList == null)
            mAnnotationList = Lists.newArrayList();

        if(mAnnotationList.contains(annotation))
            return;

        mAnnotationList.add(annotation);
    }

    public boolean hasAnnotation(final String annotation) { return mAnnotationList != null && mAnnotationList.contains(annotation); }

    public String getAnnotations()
    {
        if(mAnnotationList == null)
            return "";

        return mAnnotationList.stream().collect (Collectors.joining (ITEM_DELIM));
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
