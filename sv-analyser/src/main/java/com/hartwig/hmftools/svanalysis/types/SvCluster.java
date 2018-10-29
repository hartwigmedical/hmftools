package com.hartwig.hmftools.svanalysis.types;

import static java.lang.Math.max;
import static java.lang.Math.abs;

import static com.hartwig.hmftools.svanalysis.analysis.ClusterAnalyser.SMALL_CLUSTER_SIZE;
import static com.hartwig.hmftools.svanalysis.analysis.SvUtilities.PERMITED_DUP_BE_DISTANCE;
import static com.hartwig.hmftools.svanalysis.analysis.SvUtilities.calcConsistency;
import static com.hartwig.hmftools.svanalysis.analysis.SvUtilities.variantMatchesBreakend;
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
    private boolean mRequiresRecalc;
    private List<String> mAnnotationList;

    private List<SvVarData> mSVs;
    private List<SvBreakend> mUniqueBreakends; // for duplicate BE searches
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

    private int mMinCopyNumber;
    private int mMaxCopyNumber;


    public static String RESOLVED_TYPE_SIMPLE_SV = "SimpleSV";
    public static String RESOLVED_TYPE_RECIPROCAL_TRANS = "RecipTrans";
    // public static String RESOLVED_TYPE_RECIPROCAL_INV = "RecipInv";

    public static String RESOLVED_TYPE_NONE = "None";
    public static String RESOLVED_TYPE_DEL_INT_TI = "DEL_Int_TI";
    public static String RESOLVED_TYPE_DEL_EXT_TI = "DEL_Ext_TI";
    public static String RESOLVED_TYPE_DUP_INT_TI = "DUP_Int_TI";
    public static String RESOLVED_TYPE_DUP_EXT_TI = "DUP_Ext_TI";

    public static String RESOLVED_TYPE_SIMPLE_CHAIN = "SimpleChain";
    public static String RESOLVED_TYPE_COMPLEX_CHAIN = "ComplexChain";

    private static final Logger LOGGER = LogManager.getLogger(SvCluster.class);

    public SvCluster(final int clusterId)
    {
        mClusterId = clusterId;
        mSVs = Lists.newArrayList();
        mUniqueBreakends = Lists.newArrayList();
        mArmGroups = Lists.newArrayList();

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

    public void addVariant(final SvVarData var)
    {
        mSVs.add(var);
        mUnchainedSVs.add(var);
        mRequiresRecalc = true;
        mIsResolved = false;
        mResolvedType = RESOLVED_TYPE_NONE;
        mLengthOverride = 0;

        if(var.isReplicatedSv())
            mHasReplicatedSVs = true;

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

    public void mergeOtherCluster(final SvCluster other)
    {
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

        for(final SvVarData var : other.getSVs())
        {
            addVariant(var);
        }
    }

    public void removeReplicatedSvs()
    {
        if(!mHasReplicatedSVs)
            return;

        int index = 0;
        while (index < mSVs.size())
        {
            mSVs.get(index).setReplicatedCount(0);
            if (mSVs.get(index).isReplicatedSv())
                mSVs.remove(index);
            else
                ++index;
        }

        mHasReplicatedSVs = false;
        updateClusterDetails();
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

    public final String linkingChromosome(boolean useStart)
    {
        if(isFullyChained())
        {
            return mChains.get(0).openChromosome(useStart);
        }
        else if(!mSubClusters.isEmpty())
        {
            return mSubClusters.get(0).linkingChromosome(useStart);
        }

        String chr = "";
        for(final SvVarData var: mSVs)
        {
            if(var.isNullBreakend() && !useStart)
                return "";

            if(chr.isEmpty())
                chr = var.chromosome(useStart);
            else if(!chr.equals(var.chromosome(useStart)))
                return "";
        }

        return chr;
    }

    public final String linkingArm(boolean useStart)
    {
        if(isFullyChained())
        {
            return mChains.get(0).openArm(useStart);
        }
        else if(!mSubClusters.isEmpty())
        {
            return mSubClusters.get(0).linkingArm(useStart);
        }
        // take the arm from the SVs if they're all the same
        String arm = "";
        for(final SvVarData var: mSVs)
        {
            if(var.isNullBreakend() && !useStart)
                return "";

            if(arm.isEmpty())
                arm = var.arm(useStart);
            else if(!arm.equals(var.arm(useStart)))
                return "";
        }

        return arm;
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
    public void setLinkedPairs(final List<SvLinkedPair> pairs) { mLinkedPairs = pairs; }
    public void setInferredLinkedPairs(final List<SvLinkedPair> pairs) { mInferredLinkedPairs = pairs; }
    public void setAssemblyLinkedPairs(final List<SvLinkedPair> pairs) { mAssemblyLinkedPairs = pairs; }
    public void setSpanningSVs(final List<SvVarData> svList) { mSpanningSVs = svList; }

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
        for(final SvVarData var : cluster.getSVs())
        {
            addVariant(var);
        }

        mAssemblyLinkedPairs.addAll(cluster.getAssemblyLinkedPairs());
        mInferredLinkedPairs.addAll(cluster.getInferredLinkedPairs());

        for(SvChain chain : cluster.getChains())
        {
            addChain(chain);
        }
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

        // for now, just link to this
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
            LOGGER.info("cluster({}) simple svCount({}) desc({}) armCount({}) consistency({}) ",
                    getId(), getCount(), getDesc(), getChromosomalArmCount(), getConsistencyCount());
        }
        else
        {
            double chainedPerc = 1 - (getUnlinkedSVs().size()/mSVs.size());
            LOGGER.info(String.format("cluster(%d) complex SVs(%d rep=%d) desc(%s) arms(%d) consistency(%d) chains(%d perc=%.2f) replic(%s)",
                    getId(), getUniqueSvCount(), getCount(), getDesc(), getChromosomalArmCount(), getConsistencyCount(),
                    mChains.size(), chainedPerc, mHasReplicatedSVs));

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
        Map<String, Integer> typeMap = new HashMap<>();

        for(final SvVarData var : mSVs)
        {
            if(var.isReplicatedSv())
                continue;

            String typeStr = var.typeStr();
            if(typeMap.containsKey(typeStr))
            {
                typeMap.replace(typeStr, typeMap.get(typeStr) + 1);
            }
            else
            {
                typeMap.put(typeStr, 1);
            }
        }

        for(Map.Entry<String, Integer> entry : typeMap.entrySet())
        {
            if(!clusterTypeStr.isEmpty())
            {
                clusterTypeStr += "_";
            }

            clusterTypeStr += entry.getKey() + "=" + entry.getValue();

        }

        return clusterTypeStr;
    }

    public boolean hasLinkingLineElements()
    {
        for (final SvVarData var : mSVs)
        {
            if(var.isTranslocation() && (var.isLineElement(true) || var.isLineElement(false)))
                return true;
        }

        return false;
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

    public void setUniqueBreakends()
    {
        // group any matching BEs (same position and orientation)
        for (final SvVarData var : mSVs)
        {
            if(var.type() == StructuralVariantType.INS)
                continue;

            if(var.isReplicatedSv())
                continue;

            addVariantToUniqueBreakends(var);
        }

        // cache this against the SV
        for (SvVarData var : mSVs)
        {
            if(var.type() == StructuralVariantType.INS)
                continue;

            for(final SvBreakend breakend : mUniqueBreakends)
            {
                if(variantMatchesBreakend(var, breakend, true, PERMITED_DUP_BE_DISTANCE) && breakend.getCount() > 1)
                    var.setIsDupBEStart(true);

                if(variantMatchesBreakend(var, breakend, false, PERMITED_DUP_BE_DISTANCE) && breakend.getCount() > 1)
                {
                    if(var.type() == StructuralVariantType.INV && var.position(false) - var.position(true) <= PERMITED_DUP_BE_DISTANCE)
                    {
                        // avoid setting both end of an INV as duplicate if they match
                        continue;
                    }

                    var.setIsDupBEEnd(true);
                }
            }
        }
    }

    private void addVariantToUniqueBreakends(final SvVarData var)
    {
        for(int be = SVI_START; be <= SVI_END; ++be)
        {
            boolean useStart = isStart(be);

            if(var.type() == StructuralVariantType.INV && var.position(false) - var.position(true) <= PERMITED_DUP_BE_DISTANCE)
            {
                // avoid setting both end of an INV as duplicate if they match
                if(!useStart)
                    continue;
            }

            // ignore the null breakends
            if(var.isNullBreakend() && !useStart)
                continue;

            boolean found = false;
            for(SvBreakend breakend : mUniqueBreakends)
            {
                if (variantMatchesBreakend(var, breakend, useStart, PERMITED_DUP_BE_DISTANCE))
                {
                    breakend.addToCount(1);
                    found = true;
                    break;
                }
            }

            if(!found)
            {
                // add a new entry
                mUniqueBreakends.add(new SvBreakend(var.chromosome(useStart), var.position(useStart), var.orientation(useStart)));
            }
        }
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

    public static final SvCluster findCluster(final SvVarData var, final List<SvCluster> clusters)
    {
        for(final SvCluster cluster : clusters)
        {
            if(cluster.getSVs().contains(var))
                return cluster;
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
