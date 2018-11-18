package com.hartwig.hmftools.svanalysis.types;

import static java.lang.Math.max;
import static java.lang.Math.abs;

import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.BND;
import static com.hartwig.hmftools.svanalysis.analysis.ClusterAnalyser.SHORT_TI_LENGTH;
import static com.hartwig.hmftools.svanalysis.analysis.ClusterAnalyser.SMALL_CLUSTER_SIZE;
import static com.hartwig.hmftools.svanalysis.analysis.SvUtilities.PERMITED_DUP_BE_DISTANCE;
import static com.hartwig.hmftools.svanalysis.analysis.SvUtilities.calcConsistency;
import static com.hartwig.hmftools.svanalysis.analysis.SvUtilities.variantMatchesBreakend;
import static com.hartwig.hmftools.svanalysis.types.SvLinkedPair.ASSEMBLY_MATCH_INFER_ONLY;
import static com.hartwig.hmftools.svanalysis.types.SvVarData.SVI_END;
import static com.hartwig.hmftools.svanalysis.types.SvVarData.SVI_START;
import static com.hartwig.hmftools.svanalysis.types.SvVarData.isStart;
import static com.hartwig.hmftools.svanalysis.types.SvLinkedPair.findLinkedPair;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantType;
import com.hartwig.hmftools.svanalysis.analysis.SvUtilities;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

public class SvCluster
{
    private int mClusterId;

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
    private boolean mIsFullyChained;
    private List<SvVarData> mUnchainedSVs;
    private boolean mIsResolved;
    private String mResolvedType;
    private long mLengthOverride;

    private List<SvCluster> mSubClusters;
    private Map<String, List<SvCNData>> mChrCNDataMap;
    private boolean mHasReplicatedSVs;

    // cached lists of identified special cases
    private List<SvVarData> mLongDelDups;
    private List<SvVarData> mFoldbacks;
    private List<SvVarData> mLineElements;
    private List<SvVarData> mInversions;
    private List<SvVarData> mUnlinkedBnds;

    private int mMinCopyNumber;
    private int mMaxCopyNumber;

    public static String RESOLVED_TYPE_SIMPLE_SV = "SimpleSV";
    public static String RESOLVED_TYPE_RECIPROCAL_TRANS = "RecipTrans";
    // public static String RESOLVED_TYPE_RECIPROCAL_INV = "RecipInv";

    public static String RESOLVED_TYPE_NONE = "None";
    public static String RESOLVED_LOW_QUALITY = "LowQual";
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
        mClusterId = clusterId;
        mSVs = Lists.newArrayList();
        mArmGroups = Lists.newArrayList();
        mTypeCountMap = new HashMap();

        // annotation info
        mDesc = "";
        mConsistencyCount = 0;
        mIsConsistent = false;
        mIsResolved = false;
        mResolvedType = RESOLVED_TYPE_NONE;
        mLengthOverride = 0;
        mRequiresRecalc = true;
        mAnnotationList = Lists.newArrayList();

        // chain data
        mLinkedPairs = Lists.newArrayList();
        mAssemblyLinkedPairs= Lists.newArrayList();
        mInferredLinkedPairs = Lists.newArrayList();
        mSpanningSVs = Lists.newArrayList();
        mChains = Lists.newArrayList();
        mIsFullyChained = false;
        mUnchainedSVs = Lists.newArrayList();

        mLongDelDups = Lists.newArrayList();
        mFoldbacks = Lists.newArrayList();
        mLineElements = Lists.newArrayList();
        mInversions = Lists.newArrayList();
        mUnlinkedBnds = Lists.newArrayList();

        mSubClusters = Lists.newArrayList();
        mChrCNDataMap = null;
        mHasReplicatedSVs = false;

        mMinCopyNumber = 0;
        mMaxCopyNumber = 0;
    }

    public int getId() { return mClusterId; }

    public int getCount() { return mSVs.size(); }

    public final String getDesc() { return mDesc; }
    public final void setDesc(final String desc) { mDesc = desc; }

    public List<SvVarData> getSVs() { return mSVs; }

    // private static String LOG_SPECIFIC_VAR_ID = "83211";
    private static String LOG_SPECIFIC_VAR_ID = "";

    public void addVariant(final SvVarData var)
    {
        if(mSVs.contains(var))
        {
            LOGGER.error("cluster({}) attempting to add SV again", mClusterId, var.id());
            return;
        }

        mSVs.add(var);
        var.setCluster(this);

        if(var.inLineElement())
        {
            mLineElements.add(var);
        }
        else if(var.type() == BND && !var.isReplicatedSv())
        {
            mUnlinkedBnds.add(var);
        }

        mUnchainedSVs.add(var);
        mRequiresRecalc = true;
        mIsResolved = false;
        mResolvedType = RESOLVED_TYPE_NONE;
        mLengthOverride = 0;

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
        mChains.add(chain);

        for(SvVarData var : chain.getSvList())
        {
            mUnchainedSVs.remove(var);
        }
    }

    public void setIsFullyChained(boolean toggle) { mIsFullyChained = toggle; mUnchainedSVs.clear(); }
    public boolean isFullyChained() { return mIsFullyChained; }

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

    // private static int SPECIFIC_CLUSTER_ID = 158;
    public void mergeOtherCluster(final SvCluster other)
    {
        /*
        if(mClusterId == SPECIFIC_CLUSTER_ID)
        {
            LOGGER.debug("spec cluster");
        }
        */

        // just add the other cluster's variants - no preservation of links or chains
        if(other.getCount() > getCount())
        {
            LOGGER.debug("cluster({} svs={}) merges in other cluster({} svs={})",
                    other.getId(), other.getCount(), getId(), getCount());

            // maintain the id of the larger group
            mClusterId = other.getId();
        }
        else
        {
            LOGGER.debug("cluster({} svs={}) merges in other cluster({} svs={})",
                    getId(), getCount(), other.getId(), other.getCount());
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
    public void setLengthOverride(long length) { mLengthOverride = length; }
    public long getLengthOverride() { return mLengthOverride; }

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

    public boolean hasConsistentArms()
    {
        for(final SvArmGroup armGroup : mArmGroups)
        {
            if(calcConsistency(armGroup.getSVs()) != 0)
                return false;
        }

        return true;
    }

    public void logDetails()
    {
        updateClusterDetails();

        if(isSimpleSingleSV())
        {
            LOGGER.debug("cluster({}) simple svCount({}) desc({}) armCount({}) consistency({}) ",
                    getId(), getCount(), getDesc(), getChromosomalArmCount(), getConsistencyCount());
        }
        else
        {
            double chainedPerc = 1 - (getUnlinkedSVs().size()/mSVs.size());

            String otherInfo = "";

            if(!mFoldbacks.isEmpty())
            {
                otherInfo += String.format("foldbacks=%d", mFoldbacks.size());
            }

            if(!mLineElements.isEmpty())
            {
                if(!otherInfo.isEmpty())
                    otherInfo += " ";

                otherInfo += String.format("line=%d", mLineElements.size());
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

            if(!mUnlinkedBnds.isEmpty())
            {
                if (!otherInfo.isEmpty())
                    otherInfo += " ";

                otherInfo += String.format("bnd=%d", mUnlinkedBnds.size());
            }

            LOGGER.debug(String.format("cluster(%d) complex SVs(%d rep=%d) desc(%s res=%s) arms(%d) consis(%d) chains(%d perc=%.2f) replic(%s) %s",
                    getId(), getUniqueSvCount(), getCount(), getDesc(), mResolvedType,
                    getChromosomalArmCount(), getConsistencyCount(),
                    mChains.size(), chainedPerc, mHasReplicatedSVs, otherInfo));
        }
    }

    public int getChromosomalArmCount() { return mArmGroups.size(); }

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
        if(!mInversions.contains(var))
            mInversions.add(var);
    }

    public void registerLongDelDup(final SvVarData var)
    {
        if(!mLongDelDups.contains(var))
            mLongDelDups.add(var);
    }

    public boolean hasLinkingLineElements()
    {
        for (final SvVarData var : mLineElements)
        {
            if(var.isTranslocation() && var.inLineElement())
                return true;
        }

        return false;
    }

    public final List<SvVarData> getUnlinkedBnds() { return mUnlinkedBnds; }

    public void setUnlinkedBnds()
    {
        mUnlinkedBnds.clear();

        for(final SvVarData var : mUnchainedSVs)
        {
            if(var.type() == BND && !var.inLineElement() && !var.isReplicatedSv())
            {
                mUnlinkedBnds.add(var);
            }
        }

        for(final SvChain chain : mChains)
        {
            // any BND which doesn't form a short TI is fair game
            for(final SvLinkedPair pair : chain.getLinkedPairs())
            {
                if(pair.first().type() == BND && pair.second().type() == BND)
                {
                    if(pair.length() <= SHORT_TI_LENGTH)
                        continue;
                }

                if(pair.first().type() == BND && !pair.first().inLineElement() && !mUnlinkedBnds.contains(pair.first()))
                    mUnlinkedBnds.add(pair.first());

                if(pair.second().type() == BND && !pair.second().inLineElement()  && !mUnlinkedBnds.contains(pair.second()))
                    mUnlinkedBnds.add(pair.second());
            }
        }
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
        List<SvLinkedPair> linkedPairs;

        if(isFullyChained())
        {
            linkedPairs = mChains.get(0).getLinkedPairs();
        }
        else
        {
            linkedPairs = Lists.newArrayList();

            // add all chained links
            for(final SvChain chain : mChains)
            {
                linkedPairs.addAll(chain.getLinkedPairs());
            }

            // any any unchained assembly links
            for(final SvLinkedPair pair : mAssemblyLinkedPairs)
            {
                if(!linkedPairs.contains(pair))
                    linkedPairs.add(pair);
            }

            // finally add any other potential inferred links which don't clash with existing links
            for(final SvLinkedPair pair : mInferredLinkedPairs)
            {
                if(linkedPairs.contains(pair))
                    continue;

                boolean hasClash = false;
                for(final SvLinkedPair existingLink : linkedPairs)
                {
                    if (existingLink.hasLinkClash(pair))
                    {
                        hasClash = true;
                        break;
                    }
                }

                if(!hasClash)
                    linkedPairs.add(pair);
            }

            // mark the resultant set of inferred links
            for (SvLinkedPair pair : linkedPairs)
            {
                if(pair.isInferred())
                {
                    pair.first().setAssemblyMatchType(ASSEMBLY_MATCH_INFER_ONLY, pair.firstLinkOnStart());
                    pair.second().setAssemblyMatchType(ASSEMBLY_MATCH_INFER_ONLY, pair.secondLinkOnStart());
                }
            }
        }

        /*
        if(!linkedPairs.isEmpty())
        {
            LOGGER.debug("cluster({}: {} count={}) has {} linked pairs",
                    mClusterId, mDesc, getUniqueSvCount(), linkedPairs.size());
        }
        */

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

    public final List<String> getAnnotationList() { return mAnnotationList; }
    public final void addAnnotation(final String annotation)
    {
        if(mAnnotationList.contains(annotation))
            return;

        mAnnotationList.add(annotation);
    }

    public String getAnnotations() { return mAnnotationList.stream ().collect (Collectors.joining (";")); }

}
