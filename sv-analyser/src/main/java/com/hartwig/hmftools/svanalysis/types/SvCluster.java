package com.hartwig.hmftools.svanalysis.types;

import static java.lang.Math.max;
import static java.lang.Math.abs;

import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.BND;
import static com.hartwig.hmftools.svanalysis.analysis.ClusterAnalyser.SHORT_TI_LENGTH;
import static com.hartwig.hmftools.svanalysis.analysis.ClusterAnalyser.SMALL_CLUSTER_SIZE;
import static com.hartwig.hmftools.svanalysis.analysis.SvUtilities.calcConsistency;
import static com.hartwig.hmftools.svanalysis.types.SvLinkedPair.ASSEMBLY_MATCH_INFER_ONLY;
import static com.hartwig.hmftools.svanalysis.types.SvLinkedPair.removedLinksWithSV;
import static com.hartwig.hmftools.svanalysis.types.SvVarData.SVI_END;
import static com.hartwig.hmftools.svanalysis.types.SvVarData.SVI_START;
import static com.hartwig.hmftools.svanalysis.types.SvVarData.isStart;
import static com.hartwig.hmftools.svanalysis.types.SvLinkedPair.findLinkedPair;

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
    private Map<String, Integer> mTypeCountMap;

    private boolean mRequiresRecalc;
    private List<String> mAnnotationList;

    private List<SvVarData> mSVs;
    private List<SvChain> mChains; // pairs of SVs linked into chains
    private List<SvLinkedPair> mLinkedPairs; // final set after chaining and linking
    private List<SvLinkedPair> mInferredLinkedPairs; // forming a TI or DB
    private List<SvLinkedPair> mAssemblyLinkedPairs; // TIs found during assembly
    private List<SvVarData> mSpanningSVs; // having 2 duplicate (matching) BEs
    private List<SvArmGroup> mArmGroups;
    private List<SvVarData> mUnchainedSVs;
    private boolean mIsResolved;
    private String mResolvedType;

    // specific cluster type data

    // for synthetic DELs and DUPs
    private long mSynDelDupTI;
    private long mSynDelDupLength;

    private List<SvCluster> mSubClusters;
    private Map<String, List<SvCNData>> mChrCNDataMap;
    private boolean mHasReplicatedSVs;

    // cached lists of identified special cases
    private List<SvVarData> mLongDelDups;
    private List<SvVarData> mFoldbacks;
    private boolean mHasLinkingLineElements;
    private List<SvVarData> mInversions;

    // state for SVs which link different arms or chromosomes
    private boolean mRecalcRemoteSVStatus;
    private List<SvVarData> mShortTIRemoteSVs;
    private List<SvVarData> mUnlinkedRemoteSVs;

    private int mMinCopyNumber;
    private int mMaxCopyNumber;

    public static String RESOLVED_TYPE_SIMPLE_SV = "SimpleSV";
    public static String RESOLVED_TYPE_RECIPROCAL_TRANS = "RecipTrans";

    public static String RESOLVED_TYPE_NONE = "None";
    public static String RESOLVED_LOW_QUALITY = "LowQual";
    public static String RESOLVED_TYPE_LINE = "Line";
    public static String RESOLVED_TYPE_DEL_INT_TI = "DEL_Int_TI";
    public static String RESOLVED_TYPE_DEL_EXT_TI = "DEL_Ext_TI";
    public static String RESOLVED_TYPE_DUP_INT_TI = "DUP_Int_TI";
    public static String RESOLVED_TYPE_DUP_EXT_TI = "DUP_Ext_TI";
    public static String RESOLVED_TYPE_SGL_PAIR_INS = "SglPair_INS";
    public static String RESOLVED_TYPE_SGL_PAIR_DEL = "SglPair_DEL";
    public static String RESOLVED_TYPE_SGL_PAIR_DUP = "SglPair_DUP";
    public static String RESOLVED_TYPE_SGL_PLUS_INCONSISTENT = "SglPlusInc";

    public static String RESOLVED_TYPE_SIMPLE_CHAIN = "SimpleChain";
    public static String RESOLVED_TYPE_COMPLEX_CHAIN = "ComplexChain";

    private static final Logger LOGGER = LogManager.getLogger(SvCluster.class);

    public SvCluster(final int clusterId)
    {
        mId = clusterId;
        mSVs = Lists.newArrayList();
        mArmGroups = Lists.newArrayList();
        mTypeCountMap = new HashMap();

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

        // chain data
        mLinkedPairs = Lists.newArrayList();
        mAssemblyLinkedPairs= Lists.newArrayList();
        mInferredLinkedPairs = Lists.newArrayList();
        mSpanningSVs = Lists.newArrayList();
        mChains = Lists.newArrayList();
        mUnchainedSVs = Lists.newArrayList();

        mLongDelDups = Lists.newArrayList();
        mFoldbacks = Lists.newArrayList();
        mHasLinkingLineElements = false;
        mInversions = Lists.newArrayList();
        mShortTIRemoteSVs = Lists.newArrayList();
        mUnlinkedRemoteSVs = Lists.newArrayList();
        mRecalcRemoteSVStatus = false;

        mSubClusters = Lists.newArrayList();
        mChrCNDataMap = null;
        mHasReplicatedSVs = false;

        mMinCopyNumber = 0;
        mMaxCopyNumber = 0;
    }

    public int id() { return mId; }

    public int getCount() { return mSVs.size(); }

    public final String getDesc() { return mDesc; }
    public final void setDesc(final String desc) { mDesc = desc; }

    public List<SvVarData> getSVs() { return mSVs; }
    public final SvVarData getSV(int index) { return index < mSVs.size() ? mSVs.get(index) : null; }

    // private static String LOG_SPECIFIC_VAR_ID = "83211";
    private static String LOG_SPECIFIC_VAR_ID = "";

    public void addVariant(final SvVarData var)
    {
        if(mSVs.contains(var))
        {
            LOGGER.error("cluster({}) attempting to add SV again", mId, var.id());
            return;
        }

        mSVs.add(var);
        var.setCluster(this);

        mUnchainedSVs.add(var);
        mRequiresRecalc = true;

        if(!mHasLinkingLineElements)
        {
            mIsResolved = false;
            mResolvedType = RESOLVED_TYPE_NONE;
        }

        mSynDelDupTI = 0;
        mSynDelDupLength = 0;

        if(!LOG_SPECIFIC_VAR_ID.isEmpty() && var.id().equals(LOG_SPECIFIC_VAR_ID))
        {
            LOGGER.debug("spec var");
        }

        if(var.isReplicatedSv())
        {
            mHasReplicatedSVs = true;
        }
        else
        {
            String typeStr = var.typeStr();
            if(mTypeCountMap.containsKey(typeStr))
            {
                mTypeCountMap.replace(typeStr, mTypeCountMap.get(typeStr) + 1);
            }
            else
            {
                mTypeCountMap.put(typeStr, 1);
            }

            if(var.type() == BND || var.isCrossArm())
                mRecalcRemoteSVStatus = true;
        }

        // keep track of all SVs in their respective chromosomal arms
        for(int be = SVI_START; be <= SVI_END; ++be)
        {
            if(be == SVI_END && var.isNullBreakend())
                continue;

            if(be == SVI_END && var.isLocal())
                continue;

            boolean useStart = isStart(be);

            boolean groupFound = false;
            for(SvArmGroup armGroup : mArmGroups)
            {
                if(armGroup.chromosome().equals(var.chromosome(useStart)) && armGroup.arm().equals(var.arm(useStart)))
                {
                    armGroup.addVariant(var);
                    groupFound = true;
                    break;
                }
            }

            if(!groupFound)
            {
                SvArmGroup armGroup = new SvArmGroup(this, var.chromosome(useStart), var.arm(useStart));
                armGroup.addVariant(var);
                mArmGroups.add(armGroup);
            }
        }
    }

    public void removeReplicatedSv(final SvVarData var)
    {
        if(!var.isReplicatedSv())
            return;

        mSVs.remove(var);
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

        for(int be = SVI_START; be <= SVI_END; ++be)
        {
            if (be == SVI_END && var.isNullBreakend())
                continue;

            if (be == SVI_END && var.isLocal())
                continue;

            boolean useStart = isStart(be);

            for (SvArmGroup armGroup : mArmGroups)
            {
                if(armGroup.chromosome().equals(var.chromosome(useStart)) && armGroup.arm().equals(var.arm(useStart)))
                {
                    mArmGroups.remove(var);
                    break;
                }
            }
        }

        // retest cluster status again
        mHasReplicatedSVs = false;
        for(final SvVarData sv : mSVs)
        {
            if (sv.isReplicatedSv())
            {
                mHasReplicatedSVs = true;
                break;
            }
        }

        mRequiresRecalc = true;
    }

    public boolean isSyntheticSimpleType()
    {
        if(!mIsResolved)
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

    public boolean hasArmGroup(final String chr, final String arm) { return getArmGroup(chr, arm) != null; }

    public SvArmGroup getArmGroup(final String chr, final String arm)
    {
        for(SvArmGroup armGroup : mArmGroups)
        {
            if (armGroup.chromosome().equals(chr) && armGroup.arm().equals(arm))
                return armGroup;
        }

        return null;
    }

    public void setChrCNData(Map<String, List<SvCNData>> map) { mChrCNDataMap = map; }
    public final Map<String, List<SvCNData>> getChrCNData() { return mChrCNDataMap; }

    public int getUniqueSvCount()
    {
        if(!mHasReplicatedSVs)
            return mSVs.size();

        int count = 0;

        for(final SvVarData var : mSVs)
        {
            if(!var.isReplicatedSv())
                ++count;
        }

        return count;
    }

    public boolean hasReplicatedSVs() { return mHasReplicatedSVs; }

    public List<SvChain> getChains() { return mChains; }

    public void addChain(SvChain chain)
    {
        chain.setId(mChains.size());
        mChains.add(chain);

        for(SvVarData var : chain.getSvList())
        {
            mUnchainedSVs.remove(var);
        }
    }

    public boolean isFullyChained()
    {
        if(!mUnchainedSVs.isEmpty() || mChains.isEmpty())
            return false;

        for (final SvChain chain : mChains)
        {
            if (!chain.isConsistent())
                return false;
        }

        return true;
    }

    public List<SvVarData> getUnlinkedSVs() { return mUnchainedSVs; }

    public boolean isSimpleSingleSV() { return isSimpleSVs(); }

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
    public final List<SvVarData> getSpanningSVs() { return mSpanningSVs; }
    public void setInferredLinkedPairs(final List<SvLinkedPair> pairs) { mInferredLinkedPairs = pairs; }
    public void setAssemblyLinkedPairs(final List<SvLinkedPair> pairs) { mAssemblyLinkedPairs = pairs; }
    public void setSpanningSVs(final List<SvVarData> svList) { mSpanningSVs = svList; }

    public void mergeOtherCluster(final SvCluster other)
    {
        // just add the other cluster's variants - no preservation of links or chains
        if(other.getCount() > getCount())
        {
            LOGGER.debug("cluster({} svs={}) merges in other cluster({} svs={})",
                    other.id(), other.getCount(), id(), getCount());

            // maintain the id of the larger group
            mId = other.id();
        }
        else
        {
            LOGGER.debug("cluster({} svs={}) merges in other cluster({} svs={})",
                    id(), getCount(), other.id(), other.getCount());
        }

        addVariantLists(other);

        if(other.isFullyChained())
        {
            for (SvChain chain : other.getChains())
            {
                addChain(chain);
            }
        }

    }
    public void addSubCluster(SvCluster cluster)
    {
        if(mSubClusters.contains(cluster))
            return;

        if(cluster.hasSubClusters())
        {
            for(SvCluster subCluster : cluster.getSubClusters())
            {
                addSubCluster(subCluster);
            }

            return;
        }
        else
        {
            mSubClusters.add(cluster);
        }

        // merge the second cluster into the first
        addVariantLists(cluster);

        for(SvChain chain : cluster.getChains())
        {
            addChain(chain);
        }
    }

    private void addVariantLists(final SvCluster other)
    {
        for(final SvVarData var : other.getSVs())
        {
            addVariant(var);
        }

        mAssemblyLinkedPairs.addAll(other.getAssemblyLinkedPairs());
        mInferredLinkedPairs.addAll(other.getInferredLinkedPairs());
        mInversions.addAll(other.getInversions());
        mFoldbacks.addAll(other.getFoldbacks());
        mLongDelDups.addAll(other.getLongDelDups());
    }

    public int getMaxChainCount()
    {
        int maxCount = 0;
        for(final SvChain chain : mChains)
        {
            maxCount = max(maxCount, chain.getSvCount());
        }

        return maxCount;
    }

    public List<SvCluster> getSubClusters() { return mSubClusters; }
    public boolean hasSubClusters() { return !mSubClusters.isEmpty(); }

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

    public void setResolved(boolean toggle, final String type) { mIsResolved = toggle; mResolvedType = type; }
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

        mConsistencyCount = calcConsistency(mSVs);

        mIsConsistent = (mConsistencyCount == 0);

        mDesc = getClusterTypesAsString();

        setMinMaxCopyNumber();

        mRequiresRecalc = false;
    }

    public void logDetails()
    {
        updateClusterDetails();

        if(isSimpleSingleSV())
        {
            LOGGER.debug("cluster({}) simple svCount({}) desc({}) armCount({}) consistency({}) ",
                    id(), getCount(), getDesc(), getArmCount(), getConsistencyCount());
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
                if(!otherInfo.isEmpty())
                    otherInfo += " ";

                otherInfo += String.format("longDelDup=%d", mLongDelDups.size());
            }

            if(!mInversions.isEmpty())
            {
                if(!otherInfo.isEmpty())
                    otherInfo += " ";

                otherInfo += String.format("inv=%d", mInversions.size());
            }

            if(!mShortTIRemoteSVs.isEmpty())
            {
                if (!otherInfo.isEmpty())
                    otherInfo += " ";

                otherInfo += String.format("sti-bnd=%d", mShortTIRemoteSVs.size());
            }

            if(!mUnlinkedRemoteSVs.isEmpty())
            {
                if (!otherInfo.isEmpty())
                    otherInfo += " ";

                otherInfo += String.format("unlnk-bnd=%d", mUnlinkedRemoteSVs.size());
            }

            LOGGER.debug(String.format("cluster(%d) complex SVs(%d rep=%d) desc(%s res=%s) arms(%d) consis(%d) chains(%d perc=%.2f) replic(%s) %s",
                    id(), getUniqueSvCount(), getCount(), getDesc(), mResolvedType,
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

        // the following map-based naming convention leads
        // to a predictable ordering of types: INV, CRS, BND, DEL and DUP
        String clusterTypeStr = "";

        for(Map.Entry<String, Integer> entry : mTypeCountMap.entrySet())
        {
            if(!clusterTypeStr.isEmpty())
            {
                clusterTypeStr += "_";
            }

            clusterTypeStr += entry.getKey() + "=" + entry.getValue();

        }

        return clusterTypeStr;
    }

    public int getTypeCount(StructuralVariantType type)
    {
        Integer typeCount = mTypeCountMap.get(type.toString());
        return typeCount != null ? typeCount : 0;
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
                .filter(x -> !x.isReplicatedSv())
                .filter(x -> !mShortTIRemoteSVs.contains(x))
                .collect(Collectors.toList());

        for (final SvArmGroup armGroup : mArmGroups)
        {
            armGroup.setBoundaries(mShortTIRemoteSVs);
        }

        mRecalcRemoteSVStatus = true;
    }

    private void setMinMaxCopyNumber()
    {
        // first establish the lowest copy number change
        mMinCopyNumber = -1;
        mMaxCopyNumber = -1;

        for (final SvVarData var : mSVs)
        {
            int calcCopyNumber = var.impliedCopyNumber(true);

            if (mMinCopyNumber < 0 || calcCopyNumber < mMinCopyNumber)
                mMinCopyNumber = calcCopyNumber;

            mMaxCopyNumber = max(mMaxCopyNumber, calcCopyNumber);
        }
    }

    public int getMaxCopyNumber() { return mMaxCopyNumber; }
    public int getMinCopyNumber() { return mMinCopyNumber; }

    public boolean hasVariedCopyNumber()
    {
        if(mRequiresRecalc)
            updateClusterDetails();

        return (mMaxCopyNumber > mMinCopyNumber && mMinCopyNumber > 0);
    }

    public void cacheLinkedPairs()
    {
        // moves assembly and inferred linked pairs which are used in chains to a set of 'final' linked pairs
        // other potential inferred linked pairs are skipped
        List<SvLinkedPair> linkedPairs;

        if(isFullyChained())
        {
            linkedPairs = mChains.get(0).getLinkedPairs();
        }
        else
        {
            linkedPairs = Lists.newArrayList();

            // add all chained links
            for (final SvChain chain : mChains)
            {
                linkedPairs.addAll(chain.getLinkedPairs());
            }

            // any any unchained assembly links
            for (final SvLinkedPair pair : mAssemblyLinkedPairs)
            {
                if (!linkedPairs.contains(pair))
                {
                    linkedPairs.add(pair);
                }
            }
        }

        // mark the resultant set of inferred links - the assembly links will have already been marked
        for (SvLinkedPair pair : linkedPairs)
        {
            if(pair.isInferred())
            {
                pair.first().setAssemblyMatchType(ASSEMBLY_MATCH_INFER_ONLY, pair.firstLinkOnStart());
                pair.second().setAssemblyMatchType(ASSEMBLY_MATCH_INFER_ONLY, pair.secondLinkOnStart());
            }
        }

        mLinkedPairs = linkedPairs;
    }

    public final SvLinkedPair getLinkedPair(final SvVarData var, boolean useStart)
    {
        return findLinkedPair(mLinkedPairs, var, useStart);
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

    public final SvChain findChain(final SvLinkedPair pair)
    {
        for(final SvChain chain : mChains)
        {
            if(chain.hasLinkedPair(pair))
                return chain;
        }

        return null;
    }

    private static int SPECIFIC_CLUSTER_ID = -1;
    // private static int SPECIFIC_CLUSTER_ID = 400;

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
