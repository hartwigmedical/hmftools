package com.hartwig.hmftools.svanalysis.types;

import static java.lang.Math.ceil;
import static java.lang.Math.floor;
import static java.lang.Math.max;
import static java.lang.Math.abs;
import static java.lang.Math.round;

import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.BND;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.typeAsInt;
import static com.hartwig.hmftools.svanalysis.analysis.ClusterAnalyser.SMALL_CLUSTER_SIZE;
import static com.hartwig.hmftools.svanalysis.analysis.SvClusteringMethods.DEFAULT_PROXIMITY_DISTANCE;
import static com.hartwig.hmftools.svanalysis.analysis.SvClusteringMethods.hasLowCNChangeSupport;
import static com.hartwig.hmftools.svanalysis.analysis.SvUtilities.addSvToChrBreakendMap;
import static com.hartwig.hmftools.svanalysis.analysis.SvUtilities.appendStr;
import static com.hartwig.hmftools.svanalysis.analysis.SvUtilities.calcConsistency;
import static com.hartwig.hmftools.svanalysis.analysis.SvUtilities.getSvTypesStr;
import static com.hartwig.hmftools.svanalysis.types.SvChain.CM_CHAIN_ENDS_AWAY;
import static com.hartwig.hmftools.svanalysis.types.SvChain.CM_DB;
import static com.hartwig.hmftools.svanalysis.types.SvChain.CM_SHORT_DB;
import static com.hartwig.hmftools.svanalysis.types.SvLinkedPair.ASSEMBLY_MATCH_INFER_ONLY;
import static com.hartwig.hmftools.svanalysis.types.SvLinkedPair.ASSEMBLY_MATCH_NONE;
import static com.hartwig.hmftools.svanalysis.types.SvLinkedPair.removedLinksWithSV;
import static com.hartwig.hmftools.svanalysis.types.SvVarData.SVI_END;
import static com.hartwig.hmftools.svanalysis.types.SvVarData.SVI_START;
import static com.hartwig.hmftools.svanalysis.types.SvVarData.isStart;
import static com.hartwig.hmftools.svanalysis.types.SvaConfig.SPECIFIC_CLUSTER_ID;
import static com.hartwig.hmftools.svanalysis.types.SvaConstants.SHORT_TI_LENGTH;
import static com.hartwig.hmftools.svanalysis.types.SvaConstants.SUBCLONAL_LOW_CNC_PERCENT;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantType;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

public class SvCluster
{
    private int mId;

    private int mConsistencyCount;
    private boolean mIsConsistent; // follows from telomere to centromere to telomore
    private String mDesc;
    private int[] mTypeCounts;

    private boolean mRequiresRecalc;
    private List<String> mAnnotationList;

    private List<SvVarData> mSVs;
    private List<SvVarData> mReplicatedSVs; // combined original and replicated SV
    private List<SvChain> mChains; // pairs of SVs linked into chains
    private List<SvLinkedPair> mLinkedPairs; // final set after chaining and linking
    private List<SvLinkedPair> mInferredLinkedPairs; // forming a templated insertion
    private List<SvLinkedPair> mAssemblyLinkedPairs; // TIs found during assembly
    private List<SvArmGroup> mArmGroups;
    private List<SvArmCluster> mArmClusters; // clusters of proximate SVs on an arm, currently only used for annotations
    private Map<String, List<SvBreakend>> mChrBreakendMap;
    private List<SvVarData> mUnchainedSVs; // includes replicated SVs
    private List<SvLOH> mLohEvents;
    private boolean mIsResolved;
    private String mResolvedType;
    private String mClusteringReasons;

    // for synthetic DELs and DUPs
    private long mSynDelDupTI;
    private long mSynDelDupLength;

    private boolean mHasReplicatedSVs;

    // cached lists of identified special cases
    private List<SvVarData> mLongDelDups;
    private List<SvVarData> mFoldbacks;
    private boolean mHasLinkingLineElements;
    private boolean mIsSubclonal;
    private List<SvVarData> mInversions;

    // state for SVs which link different arms or chromosomes
    private boolean mRecalcRemoteSVStatus;
    private List<SvVarData> mShortTIRemoteSVs;
    private List<SvVarData> mUnlinkedRemoteSVs;

    private double mMinCNChange;
    private double mMaxCNChange;
    private double mValidAllelePloidySegmentPerc;

    private int mOriginArms;
    private int mFragmentArms;

    public static String RESOLVED_TYPE_NONE = "NONE";
    public static String RESOLVED_TYPE_SIMPLE_SV = "SIMPLE";
    public static String RESOLVED_TYPE_RECIPROCAL_TRANS = "RECIP_TRANS";
    public static String RESOLVED_TYPE_RECIPROCAL_INV = "RECIP_INV";
    public static String RESOLVED_TYPE_DUP_BE = "DUP_BE";
    public static String RESOLVED_TYPE_POLY_G_C = "POLY_G_C";
    public static String RESOLVED_TYPE_LINE = "LINE";
    public static String RESOLVED_TYPE_DEL_INT_TI = "DEL_INT_TI";
    public static String RESOLVED_TYPE_DEL_EXT_TI = "DEL_EXT_TI";
    public static String RESOLVED_TYPE_DUP_INT_TI = "DUP_INT_TI";
    public static String RESOLVED_TYPE_DUP_EXT_TI = "DUP_EXT_TI";
    public static String RESOLVED_TYPE_SGL_PAIR_INS = "SGL_PAIR_INS";
    public static String RESOLVED_TYPE_SGL_PAIR_DEL = "SGL_PAIR_DEL";
    public static String RESOLVED_TYPE_SGL_PAIR_DUP = "SGL_PAIR_DUP";
    public static String RESOLVED_TYPE_SGL_PLUS_INCONSISTENT = "SGL_BND_INV";
    public static String RESOLVED_TYPE_COMPLEX = "COMPLEX";

    private static final Logger LOGGER = LogManager.getLogger(SvCluster.class);

    public SvCluster(final int clusterId)
    {
        mId = clusterId;
        mSVs = Lists.newArrayList();
        mReplicatedSVs = Lists.newArrayList();
        mArmGroups = Lists.newArrayList();
        mArmClusters = Lists.newArrayList();
        mTypeCounts = new int[StructuralVariantType.values().length];

        // annotation info
        mDesc = "";
        mConsistencyCount = 0;
        mIsConsistent = false;
        mIsResolved = false;
        mResolvedType = RESOLVED_TYPE_NONE;
        mSynDelDupTI = 0;
        mSynDelDupLength = 0;
        mRequiresRecalc = true;
        mAnnotationList = Lists.newArrayList();
        mChrBreakendMap = new HashMap();

        // chain data
        mLinkedPairs = Lists.newArrayList();
        mAssemblyLinkedPairs= Lists.newArrayList();
        mInferredLinkedPairs = Lists.newArrayList();
        mChains = Lists.newArrayList();
        mUnchainedSVs = Lists.newArrayList();
        mLohEvents = Lists.newArrayList();

        mLongDelDups = Lists.newArrayList();
        mFoldbacks = Lists.newArrayList();
        mHasLinkingLineElements = false;
        mInversions = Lists.newArrayList();
        mShortTIRemoteSVs = Lists.newArrayList();
        mUnlinkedRemoteSVs = Lists.newArrayList();
        mRecalcRemoteSVStatus = false;
        mIsSubclonal = false;

        mHasReplicatedSVs = false;

        mMinCNChange = 0;
        mMaxCNChange = 0;
        mValidAllelePloidySegmentPerc = 1.0;
        mClusteringReasons = "";

        mOriginArms = 0;
        mFragmentArms = 0;
    }

    public int id() { return mId; }

    public int getSvCount() { return mSVs.size(); }
    public int getSvCount(boolean includeReplicated) { return includeReplicated ? mReplicatedSVs.size() : mSVs.size(); }

    public final String getDesc() { return mDesc; }
    public final void setDesc(final String desc) { mDesc = desc; }

    public final List<SvVarData> getSVs() { return mSVs; }
    public final List<SvVarData> getSVs(boolean includeReplicated) { return includeReplicated ? mReplicatedSVs : mSVs; }
    public final SvVarData getSV(int index) { return index < mSVs.size() ? mSVs.get(index) : null; }

    public void addVariant(final SvVarData var)
    {
        if(mSVs.contains(var))
        {
            LOGGER.error("cluster({}) attempting to add SV again", mId, var.id());
            return;
        }

        if(!var.isReplicatedSv())
        {
            mSVs.add(var);
        }

        var.setCluster(this);

        mReplicatedSVs.add(var);
        mUnchainedSVs.add(var);

        mRequiresRecalc = true;

        if(!mHasLinkingLineElements)
        {
            mIsResolved = false;
            mResolvedType = RESOLVED_TYPE_NONE;
        }

        mSynDelDupTI = 0;
        mSynDelDupLength = 0;

        // isSpecificSV(var.id())

        if(var.isReplicatedSv())
        {
            mHasReplicatedSVs = true;
        }
        else
        {
            ++mTypeCounts[typeAsInt(var.type())];

            if (var.type() == BND || var.isCrossArm())
                mRecalcRemoteSVStatus = true;

            addSvToChrBreakendMap(var, mChrBreakendMap);

            // keep track of all SVs in their respective chromosomal arms
            for (int be = SVI_START; be <= SVI_END; ++be)
            {
                if (be == SVI_END && var.isNullBreakend())
                    continue;

                if (be == SVI_END && var.isLocal())
                    continue;

                boolean useStart = isStart(be);

                boolean groupFound = false;
                for (SvArmGroup armGroup : mArmGroups)
                {
                    if (armGroup.chromosome().equals(var.chromosome(useStart)) && armGroup.arm().equals(var.arm(useStart)))
                    {
                        armGroup.addVariant(var);
                        groupFound = true;
                        break;
                    }
                }

                if (!groupFound)
                {
                    SvArmGroup armGroup = new SvArmGroup(this, var.chromosome(useStart), var.arm(useStart));
                    armGroup.addVariant(var);
                    mArmGroups.add(armGroup);
                }
            }
        }
    }

    public void removeReplicatedSvs()
    {
        if(!mHasReplicatedSVs)
            return;

        int i = 0;
        while(i < mReplicatedSVs.size())
        {
            SvVarData var = mReplicatedSVs.get(i);

            if(var.isReplicatedSv())
            {
                removeReplicatedSv(var);
            }
            else
            {
                ++i;
            }
        }
    }

    public void removeReplicatedSv(final SvVarData var)
    {
        if(!var.isReplicatedSv())
            return;

        mReplicatedSVs.remove(var);
        mUnchainedSVs.remove(var);

        // remove from any cached links
        removedLinksWithSV(mAssemblyLinkedPairs, var);
        removedLinksWithSV(mInferredLinkedPairs, var);
        removedLinksWithSV(mLinkedPairs, var);

        // deregister from the original SV
        if(var.getReplicatedSv() != null)
        {
            int newReplicationCount = max(var.getReplicatedSv().getReplicatedCount() - 1, 0);
            var.getReplicatedSv().setReplicatedCount(newReplicationCount);
        }

        // retest cluster status again
        mHasReplicatedSVs = mReplicatedSVs.size() > mSVs.size();
        mRequiresRecalc = true;
    }

    public boolean isSyntheticSimpleType(boolean checkResolved)
    {
        if(checkResolved && !mIsResolved)
            return false;

        if(mResolvedType == RESOLVED_TYPE_DEL_EXT_TI || mResolvedType == RESOLVED_TYPE_DEL_INT_TI
        || mResolvedType == RESOLVED_TYPE_DUP_INT_TI || mResolvedType == RESOLVED_TYPE_DUP_EXT_TI)
        {
            return true;
        }

        if(mResolvedType == RESOLVED_TYPE_SGL_PAIR_DEL || mResolvedType == RESOLVED_TYPE_SGL_PAIR_DEL
        || mResolvedType == RESOLVED_TYPE_SGL_PAIR_DUP)
        {
            return true;
        }

        return false;
    }

    public List<SvArmGroup> getArmGroups() { return mArmGroups; }
    public Map<String, List<SvBreakend>> getChrBreakendMap() { return mChrBreakendMap; }

    public boolean hasReplicatedSVs() { return mHasReplicatedSVs; }

    public void addLohEvent(final SvLOH lohEvent)
    {
        if(!mLohEvents.contains(lohEvent))
            mLohEvents.add(lohEvent);
    }

    public final List<SvLOH> getLohEvents() { return mLohEvents; }

    public List<SvChain> getChains() { return mChains; }

    public void addChain(SvChain chain, boolean resetId)
    {
        if(resetId)
            chain.setId(mChains.size());

        mChains.add(chain);

        for(SvVarData var : chain.getSvList())
        {
            mUnchainedSVs.remove(var);
        }
    }

    public boolean isFullyChained(boolean requireConsistency)
    {
        if(!mUnchainedSVs.isEmpty() || mChains.isEmpty())
            return false;

        if(requireConsistency)
        {
            for (final SvChain chain : mChains)
            {
                if (!chain.isConsistent())
                    return false;
            }
        }

        return true;
    }

    public void dissolveLinksAndChains()
    {
        mUnchainedSVs.clear();
        mUnchainedSVs.addAll(mReplicatedSVs);
        mChains.clear();
        mLinkedPairs.clear();
        mInferredLinkedPairs.clear();

        mAssemblyLinkedPairs.clear();
        for(final SvVarData var : mSVs)
        {
            var.setAssemblyMatchType(ASSEMBLY_MATCH_NONE, true);
            var.setAssemblyMatchType(ASSEMBLY_MATCH_NONE, false);
        }
    }

    public List<SvVarData> getUnlinkedSVs() { return mUnchainedSVs; }

    public boolean isSimpleSingleSV() { return getSvCount() == 1 && isSimpleSVs(); }

    public boolean isSimpleSVs()
    {
        if(mSVs.size() > SMALL_CLUSTER_SIZE)
            return false;

        for(final SvVarData var : mSVs)
        {
            if(!var.isSimpleType())
                return false;
        }

        return true;
    }

    public final List<SvLinkedPair> getLinkedPairs() { return mLinkedPairs; }
    public final List<SvLinkedPair> getInferredLinkedPairs() { return mInferredLinkedPairs; }
    public final List<SvLinkedPair> getAssemblyLinkedPairs() { return mAssemblyLinkedPairs; }
    public void setInferredLinkedPairs(final List<SvLinkedPair> pairs) { mInferredLinkedPairs = pairs; }
    public void setAssemblyLinkedPairs(final List<SvLinkedPair> pairs) { mAssemblyLinkedPairs = pairs; }

    public void mergeOtherCluster(final SvCluster other)
    {
        mergeOtherCluster(other, true);
    }

    public void mergeOtherCluster(final SvCluster other, boolean logDetails)
    {
        if(other == this || other.id() == id())
        {
            LOGGER.error("attempting to merge same cluster({})", id());
            return;
        }

        // just add the other cluster's variants - no preservation of links or chains
        if(other.getSvCount() > getSvCount())
        {
            if(logDetails)
            {
                LOGGER.debug("cluster({} svs={}) merges in other cluster({} svs={}) and adopts new ID",
                        id(), getSvCount(), other.id(), other.getSvCount());
            }

            // maintain the id of the larger group
            mId = other.id();
        }
        else if(logDetails)
        {
            LOGGER.debug("cluster({} svs={}) merges in other cluster({} svs={})",
                    id(), getSvCount(), other.id(), other.getSvCount());
        }

        addVariantLists(other);

        other.getLohEvents().stream().forEach(this::addLohEvent);

        String[] clusterReasons = other.getClusteringReasons().split(";");
        for(int i = 0; i < clusterReasons.length; ++i)
        {
            addClusterReason(clusterReasons[i]);
        }

        if(other.hasLinkingLineElements())
        {
            // retain status as a LINE cluster
            markAsLine();
        }
    }

    private void addVariantLists(final SvCluster other)
    {
        for(final SvVarData var : other.getSVs(true))
        {
            addVariant(var);
        }

        mAssemblyLinkedPairs.addAll(other.getAssemblyLinkedPairs());
        mInferredLinkedPairs.addAll(other.getInferredLinkedPairs());
        mInversions.addAll(other.getInversions());
        mFoldbacks.addAll(other.getFoldbacks());
        mLongDelDups.addAll(other.getLongDelDups());
    }

    public void addClusterReason(final String reason)
    {
        if(!mClusteringReasons.contains(reason))
        {
            mClusteringReasons = appendStr(mClusteringReasons, reason, ';');
        }
    }

    public final String getClusteringReasons() { return mClusteringReasons; }

    public boolean isConsistent()
    {
        updateClusterDetails();
        return mIsConsistent;
    }

    public int getConsistencyCount()
    {
        updateClusterDetails();
        return mConsistencyCount;
    }

    public void setResolved(boolean toggle, final String type)
    {
        mIsResolved = toggle;
        mResolvedType = type;

        if(mDesc.isEmpty())
            mDesc = getClusterTypesAsString();
    }

    public boolean isResolved() { return mIsResolved; }
    public final String getResolvedType() { return mResolvedType; }

    public void setSynDelDupData(long length, long tiLength)
    {
        mSynDelDupLength = length;
        mSynDelDupTI = tiLength;
    }

    public long getSynDelDupTILength() { return mSynDelDupTI; }
    public long getSynDelDupLength() { return mSynDelDupLength; }

    private void updateClusterDetails()
    {
        if(!mRequiresRecalc)
            return;

        mConsistencyCount = calcConsistency(mReplicatedSVs);

        mIsConsistent = (mConsistencyCount == 0);

        mDesc = getClusterTypesAsString();

        setMinMaxCNChange();

        resetBreakendMapIndices();

        mRequiresRecalc = false;
    }

    public void logDetails()
    {
        updateClusterDetails();

        if(isSimpleSVs())
        {
            LOGGER.debug("cluster({}) simple svCount({}) desc({}) armCount({}) consistency({}) ",
                    id(), getSvCount(), getDesc(), getArmCount(), getConsistencyCount());
        }
        else
        {
            double chainedPerc = 1 - (getUnlinkedSVs().size()/mSVs.size());

            String otherInfo = "";

            if(!mFoldbacks.isEmpty())
            {
                otherInfo += String.format("foldbacks=%d", mFoldbacks.size());
            }

            if(!mLongDelDups.isEmpty())
            {
                otherInfo = appendStr(otherInfo, String.format("longDelDup=%d", mLongDelDups.size()), ' ');
            }

            if(!mInversions.isEmpty())
            {
                otherInfo = appendStr(otherInfo, String.format("inv=%d", mInversions.size()), ' ');
            }

            if(!mShortTIRemoteSVs.isEmpty())
            {
                otherInfo = appendStr(otherInfo, String.format("sti-bnd=%d", mShortTIRemoteSVs.size()), ' ');
            }

            if(!mUnlinkedRemoteSVs.isEmpty())
            {
                otherInfo = appendStr(otherInfo, String.format("unlnk-bnd=%d", mUnlinkedRemoteSVs.size()), ' ');
            }

            LOGGER.debug(String.format("cluster(%d) complex SVs(%d rep=%d) desc(%s res=%s) arms(%d) consis(%d) chains(%d perc=%.2f) replic(%s) %s",
                    id(), getSvCount(), getSvCount(true), getDesc(), mResolvedType,
                    getArmCount(), getConsistencyCount(),
                    mChains.size(), chainedPerc, mHasReplicatedSVs, otherInfo));
        }
    }

    public int getArmCount() { return mArmGroups.size(); }

    public final String getClusterTypesAsString()
    {
        if(mSVs.size() == 1)
        {
            return mSVs.get(0).typeStr();
        }

        return getSvTypesStr(mTypeCounts);
    }

    public int getTypeCount(StructuralVariantType type)
    {
        return mTypeCounts[typeAsInt(type)];
    }

    public final List<SvVarData> getLongDelDups() { return mLongDelDups; }
    public final List<SvVarData> getFoldbacks() { return mFoldbacks; }
    public final List<SvVarData> getInversions() { return mInversions; }

    public void registerFoldback(final SvVarData var)
    {
        if(!mFoldbacks.contains(var))
            mFoldbacks.add(var);
    }

    public void deregisterFoldback(final SvVarData var)
    {
        if(mFoldbacks.contains(var))
            mFoldbacks.remove(var);
    }

    public void registerInversion(final SvVarData var)
    {
        if(!mInversions.contains(var) && !var.isReplicatedSv())
            mInversions.add(var);
    }

    public void registerLongDelDup(final SvVarData var)
    {
        if(!mLongDelDups.contains(var) && !var.isReplicatedSv())
            mLongDelDups.add(var);
    }

    public void markAsLine()
    {
        mHasLinkingLineElements = true;
        setResolved(true, RESOLVED_TYPE_LINE);
    }

    public boolean hasLinkingLineElements() { return mHasLinkingLineElements; }

    public void markSubclonal()
    {
        long lowCNChangeSupportCount = mSVs.stream().filter(x -> hasLowCNChangeSupport(x)).count();
        mIsSubclonal = lowCNChangeSupportCount / (double)mSVs.size() > SUBCLONAL_LOW_CNC_PERCENT;
    }

    public boolean isSubclonal() { return mIsSubclonal; }

    public final List<SvVarData> getUnlinkedRemoteSVs() { return mUnlinkedRemoteSVs; }
    public final List<SvVarData> getShortTIRemoteSVs() { return mShortTIRemoteSVs; }

    public void setArmLinks()
    {
        if(!mRecalcRemoteSVStatus)
            return;

        // keep track of BND which are or aren't candidates for links between arms
        mShortTIRemoteSVs.clear();
        mUnlinkedRemoteSVs.clear();

        for (final SvChain chain : mChains)
        {
            // any pair of remote SVs which don't form a short TI are fair game
            for (final SvLinkedPair pair : chain.getLinkedPairs())
            {
                if (pair.first().isCrossArm() && pair.second().isCrossArm() && pair.length() <= SHORT_TI_LENGTH)
                {
                    if (!mShortTIRemoteSVs.contains(pair.first()))
                        mShortTIRemoteSVs.add(pair.first());

                    if (!mShortTIRemoteSVs.contains(pair.second()))
                        mShortTIRemoteSVs.add(pair.second());
                }
            }
        }

        mUnlinkedRemoteSVs = mSVs.stream()
                .filter(x -> x.isCrossArm())
                .filter(x -> !x.inLineElement())
                .filter(x -> !mShortTIRemoteSVs.contains(x))
                .collect(Collectors.toList());

        for (final SvArmGroup armGroup : mArmGroups)
        {
            armGroup.setBoundaries(mShortTIRemoteSVs);
        }

        mRecalcRemoteSVStatus = true;
    }

    private void setMinMaxCNChange()
    {
        if(mSVs.size() == 1)
        {
            mMinCNChange = mMaxCNChange = mSVs.get(0).getRoundedCNChange();
            return;
        }

        // establish the lowest copy number change, using calculated ploidy if present
        mMinCNChange = -1;
        mMaxCNChange = 0;

        int svCalcPloidyCount = 0;
        Map<Integer,Integer> ploidyFrequency = new HashMap();

        // isSpecificCluster(this);

        double tightestMinPloidy = 0;
        double tightestMaxPloidy = -1;
        int countHalfToOnePloidy = 0;
        double minCopyNumberChange = -1;

        for (final SvVarData var : mSVs)
        {
            double calcCNChange = var.getRoundedCNChange();

            if(calcCNChange <= 0)
            {
                LOGGER.debug("cluster({}) has SV({}) with zero effective CN change: {}",
                        mId, var.id(), String.format("start=%.2f end=%.2f p=%.2f pMax=%.2f pMin=%.2f",
                                var.copyNumberChange(true), var.copyNumberChange(false),
                                var.getSvData().ploidy(), var.ploidyMax(), var.ploidyMin()));
            }

            if (mMinCNChange < 0 || calcCNChange < mMinCNChange)
            {
                mMinCNChange = calcCNChange;
                minCopyNumberChange = calcCNChange;
            }

            mMaxCNChange = max(mMaxCNChange, calcCNChange);

            if(var.hasCalculatedPloidy())
            {
                ++svCalcPloidyCount;

                int minPloidyInt = (int)ceil(var.ploidyMin());
                int maxPloidyInt = (int)floor(var.ploidyMax());
                maxPloidyInt = max(minPloidyInt, maxPloidyInt);

                if(tightestMaxPloidy == -1 || var.ploidyMax() < tightestMaxPloidy)
                    tightestMaxPloidy = var.ploidyMax();

                tightestMinPloidy = max(var.ploidyMin(), tightestMinPloidy);

                if(var.ploidyMin() < 1 && var.ploidyMax() > 0.5)
                    ++countHalfToOnePloidy;

                for(int i = minPloidyInt; i <= maxPloidyInt; ++i)
                {
                    Integer svCount = ploidyFrequency.get(i);
                    if(svCount == null)
                        ploidyFrequency.put(i, 1);
                    else
                        ploidyFrequency.put(i, svCount+1);
                }
            }
        }

        if(svCalcPloidyCount > 0)
        {
            mMinCNChange = -1;
            mMaxCNChange = 0;

            for (Map.Entry<Integer, Integer> entry : ploidyFrequency.entrySet())
            {
                int ploidy = entry.getKey();
                int svCount = entry.getValue();

                if (svCount == svCalcPloidyCount)
                {
                    // all SVs can settle on the same ploidy value, so take this
                    mMaxCNChange = ploidy;
                    mMinCNChange = ploidy;
                    break;
                }

                if (ploidy > 0 && (mMinCNChange < 0 || ploidy < mMinCNChange))
                     mMinCNChange = ploidy;

                mMinCNChange = max(mMinCNChange, minCopyNumberChange);

                mMaxCNChange = max(mMaxCNChange, ploidy);
            }

            if(mMinCNChange < mMaxCNChange)
            {
                if (tightestMaxPloidy > tightestMinPloidy && tightestMaxPloidy - tightestMinPloidy < 1)
                {
                    // if all SVs cover the same value but it's not an integer, still consider them uniform
                    mMinCNChange = 1;
                    mMaxCNChange = 1;
                }
                else if (countHalfToOnePloidy == svCalcPloidyCount)
                {
                    mMinCNChange = 1;
                    mMaxCNChange = 1;
                }
            }

        }
    }

    public double getMaxCNChange() { return mMaxCNChange; }
    public double getMinCNChange() { return mMinCNChange; }

    public boolean hasVariedCopyNumber()
    {
        if(mRequiresRecalc)
            updateClusterDetails();

        if(mSVs.size() == 1)
            return false;

        return (mMaxCNChange > mMinCNChange && mMinCNChange >= 0);
    }

    public double getValidAllelePloidySegmentPerc() { return mValidAllelePloidySegmentPerc; }
    public void setValidAllelePloidySegmentPerc(double percent) { mValidAllelePloidySegmentPerc = percent; }

    private void resetBreakendMapIndices()
    {
        for (Map.Entry<String, List<SvBreakend>> entry : mChrBreakendMap.entrySet())
        {
            List<SvBreakend> breakendList = entry.getValue();

            for (int i = 0; i < breakendList.size(); ++i)
            {
                final SvBreakend breakend = breakendList.get(i);
                breakend.setClusterChrPosIndex(i);
            }
        }
    }

    public void cacheLinkedPairs()
    {
        // moves assembly and inferred linked pairs which are used in chains to a set of 'final' linked pairs
        mLinkedPairs.clear();

        // add all chained links
        for (final SvChain chain : mChains)
        {
            mLinkedPairs.addAll(chain.getLinkedPairs());
        }

        for(SvVarData var : mReplicatedSVs)
        {
            for(int be = SVI_START; be <= SVI_END; ++be)
            {
                boolean useStart = isStart(be);

                if(var.isNullBreakend() && !useStart)
                    continue;

                if (var.getLinkedPair(useStart) != null && var.getLinkedPair(useStart).isInferred())
                    var.setLinkedPair(null, useStart);
            }
        }

        // mark the resultant set of inferred links - the assembly links will have already been marked
        for (SvLinkedPair pair : mLinkedPairs)
        {
            if(pair.isInferred())
            {
                pair.first().setAssemblyMatchType(ASSEMBLY_MATCH_INFER_ONLY, pair.firstLinkOnStart());
                pair.second().setAssemblyMatchType(ASSEMBLY_MATCH_INFER_ONLY, pair.secondLinkOnStart());
            }

            pair.first().setLinkedPair(pair, pair.firstLinkOnStart());
            pair.second().setLinkedPair(pair, pair.secondLinkOnStart());
        }
    }

    public final SvChain findChain(final SvVarData var)
    {
        for(final SvChain chain : mChains)
        {
            if(chain.getSvIndex(var) >= 0)
                return chain;
        }

        return null;
    }

    public final List<SvChain> findChains(final SvVarData var)
    {
        return mChains.stream().filter(x -> x.hasSV(var, true)).collect(Collectors.toList());
    }

    public final SvChain findSameChainForSVs(SvVarData var1, SvVarData var2)
    {
        List<SvChain> chains1 = findChains(var1);
        List<SvChain> chains2 = findChains(var2);

        for(SvChain chain1 : chains1)
        {
            for(SvChain chain2 : chains2)
            {
                if(chain1 == chain2)
                    return chain1;
            }
        }

        return null;
    }

    public int getChainId(final SvVarData var)
    {
        final SvChain chain = findChain(var);

        if(chain != null)
            return chain.id();

        // otherwise set an id based on index in the unchained variants list
        for(int i = 0; i < mUnchainedSVs.size(); ++i)
        {
            final SvVarData unchainedSv = mUnchainedSVs.get(i);

            if(unchainedSv.equals(var, true))
                return mChains.size() + i + 1;
        }

        return var.dbId();
    }

    public final List<SvArmCluster> getArmClusters() { return mArmClusters; }

    public void buildArmClusters()
    {
        isSpecificCluster(this);

        for (Map.Entry<String, List<SvBreakend>> entry : mChrBreakendMap.entrySet())
        {
            List<SvBreakend> breakendList = entry.getValue();

            SvArmCluster prevArmCluster = null;

            for (int i = 0; i < breakendList.size(); ++i)
            {
                final SvBreakend breakend = breakendList.get(i);
                SvVarData var = breakend.getSV();

                // ensure that a pair of foldback breakends are put into the same arm cluster
                if(var.isFoldback() && var.getFoldbackBreakend(breakend.usesStart()) != null)
                {
                    SvBreakend otherFoldbackBreakend = var.getFoldbackBreakend(breakend.usesStart());
                    SvArmCluster existingAC = findArmCluster(otherFoldbackBreakend);

                    if(existingAC != null)
                    {
                        existingAC.addBreakend(breakend);
                        continue;
                    }
                }

                // first test the previous arm cluster
                if(prevArmCluster != null)
                {
                    if(breakend.arm() == prevArmCluster.arm() && breakend.position() - prevArmCluster.posEnd() <= DEFAULT_PROXIMITY_DISTANCE)
                    {
                        prevArmCluster.addBreakend(breakend);
                        continue;
                    }

                    // prevArmCluster = null;
                }

                boolean groupFound = false;

                for(final SvArmCluster armCluster : mArmClusters)
                {
                    if(!breakend.chromosome().equals(armCluster.chromosome()) || breakend.arm() != armCluster.arm())
                        continue;

                    // test whether position is within range
                    if(breakend.position() >= armCluster.posStart() - DEFAULT_PROXIMITY_DISTANCE
                    && breakend.position() <= armCluster.posEnd() + DEFAULT_PROXIMITY_DISTANCE)
                    {
                        armCluster.addBreakend(breakend);
                        groupFound = true;
                        prevArmCluster = armCluster;
                        break;
                    }
                }

                if(!groupFound)
                {
                    SvArmCluster armCluster = new SvArmCluster(mArmClusters.size(), this, breakend.chromosome(), breakend.arm());
                    armCluster.addBreakend(breakend);
                    mArmClusters.add(armCluster);
                    prevArmCluster = armCluster;
                }
            }
        }

        // isSpecificCluster(this);
        mArmClusters.forEach(x -> x.setFeatures());
    }

    public SvArmCluster findArmCluster(final SvBreakend breakend)
    {
        for(final SvArmCluster armCluster : mArmClusters)
        {
            if(armCluster.getBreakends().contains(breakend))
                return armCluster;
        }

        return null;
    }

    public void setArmData(int origins, int fragments) { mOriginArms = origins; mFragmentArms = fragments; }
    public int getOriginArms() { return mOriginArms; }
    public int getFragmentArms() { return mFragmentArms; }

    public int[] getLinkMetrics()
    {
        int[] chainData = new int[CM_CHAIN_ENDS_AWAY +1];
        mChains.stream().forEach(x -> x.extractChainMetrics(chainData));

        for(SvVarData var : mSVs)
        {
            for(int be = SVI_START; be <= SVI_END; ++be)
            {
                if (be == SVI_END && var.isNullBreakend())
                    continue;

                boolean useStart = isStart(be);
                SvLinkedPair dbLink = var.getDBLink(useStart);

                // only take matches on the lower breakend to avoid double-counting DBs
                if(dbLink != null && dbLink.getBreakend(true).getSV() == var && dbLink.getOtherSV(var).getCluster() == this)
                {
                    ++chainData[CM_DB];

                    if(dbLink.length() <= 100)
                        ++chainData[CM_SHORT_DB];
                }
            }
        }

        return chainData;
    }

    public static boolean isSpecificCluster(final SvCluster cluster)
    {
        if(cluster.id() == SPECIFIC_CLUSTER_ID)
            return true;

        return false;
    }

    // private static int SPECIFIC_CLUSTER_ID_2 = 29;
    private static int SPECIFIC_CLUSTER_ID_2 = -1;

    public static boolean areSpecificClusters(final SvCluster cluster1, final SvCluster cluster2)
    {
        if((cluster1.id() == SPECIFIC_CLUSTER_ID && cluster2.id() == SPECIFIC_CLUSTER_ID_2)
        || (cluster2.id() == SPECIFIC_CLUSTER_ID && cluster1.id() == SPECIFIC_CLUSTER_ID_2))
            return true;

        return false;
    }

    public static String CLUSTER_ANNONTATION_DM = "DM";
    public static String CLUSTER_ANNONTATION_CT = "CT";

    public final List<String> getAnnotationList() { return mAnnotationList; }
    public final void addAnnotation(final String annotation)
    {
        if(mAnnotationList.contains(annotation))
            return;

        mAnnotationList.add(annotation);
    }

    public String getAnnotations() { return mAnnotationList.stream().collect (Collectors.joining (";")); }

}
