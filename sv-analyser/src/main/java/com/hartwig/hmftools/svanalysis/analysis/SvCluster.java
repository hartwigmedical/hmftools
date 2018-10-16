package com.hartwig.hmftools.svanalysis.analysis;

import static java.lang.Math.max;
import static java.lang.Math.abs;

import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.BND;
import static com.hartwig.hmftools.svanalysis.analysis.ClusterAnalyser.SMALL_CLUSTER_SIZE;
import static com.hartwig.hmftools.svanalysis.analysis.SvUtilities.PERMITED_DUP_BE_DISTANCE;
import static com.hartwig.hmftools.svanalysis.analysis.SvUtilities.variantMatchesBreakend;
import static com.hartwig.hmftools.svanalysis.annotators.LineElementAnnotator.NO_LINE_ELEMENT;
import static com.hartwig.hmftools.svanalysis.types.SvClusterData.SVI_END;
import static com.hartwig.hmftools.svanalysis.types.SvClusterData.SVI_START;
import static com.hartwig.hmftools.svanalysis.types.SvClusterData.isStart;
import static com.hartwig.hmftools.svanalysis.types.SvLinkedPair.findLinkedPair;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantType;
import com.hartwig.hmftools.svanalysis.types.SvArmGroup;
import com.hartwig.hmftools.svanalysis.types.SvBreakend;
import com.hartwig.hmftools.svanalysis.types.SvCNData;
import com.hartwig.hmftools.svanalysis.types.SvChain;
import com.hartwig.hmftools.svanalysis.types.SvClusterData;
import com.hartwig.hmftools.svanalysis.types.SvLinkedPair;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

public class SvCluster
{
    final SvUtilities mUtils;

    private int mClusterId;

    private int mConsistencyCount;
    private boolean mIsConsistent; // follows from telomere to centromere to telomore
    private String mDesc;
    private boolean mRequiresRecalc;
    private List<String> mAnnotationList;

    private List<SvClusterData> mSVs;
    private List<SvBreakend> mUniqueBreakends; // for duplicate BE searches
    private List<SvChain> mChains; // pairs of SVs linked into chains
    private List<SvLinkedPair> mLinkedPairs; // final set after chaining and linking
    private List<SvLinkedPair> mInferredLinkedPairs; // forming a TI or DB
    private List<SvLinkedPair> mAssemblyLinkedPairs; // TIs found during assembly
    private List<SvClusterData> mSpanningSVs; // having 2 duplicate (matching) BEs
    private List<SvArmGroup> mArmGroups;
    private boolean mIsFullyChained;
    private List<SvClusterData> mUnchainedSVs;

    private List<SvCluster> mSubClusters;
    private Map<String, List<SvCNData>> mChrCNDataMap;
    private boolean mHasReplicatedSVs;

    private int mMinCopyNumber;
    private int mMaxCopyNumber;

    private static final Logger LOGGER = LogManager.getLogger(SvCluster.class);

    public SvCluster(final int clusterId, final SvUtilities utils)
    {
        mClusterId = clusterId;
        mUtils = utils;
        mSVs = Lists.newArrayList();
        mUniqueBreakends = Lists.newArrayList();
        mArmGroups = Lists.newArrayList();

        // annotation info
        mConsistencyCount = 0;
        mIsConsistent = false;
        mDesc = "";
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

    public List<SvClusterData> getSVs() { return mSVs; }

    public void addVariant(final SvClusterData var)
    {
        mSVs.add(var);
        mUnchainedSVs.add(var);
        mRequiresRecalc = true;

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

    public boolean isSingleArm() { return mArmGroups.size() == 1; }

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
        for(final SvClusterData var: mSVs)
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
        for(final SvClusterData var: mSVs)
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

        for(final SvClusterData var : mSVs)
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

        for(SvClusterData var : chain.getSvList())
        {
            mUnchainedSVs.remove(var);
        }
    }

    public void setIsFullyChained(boolean toggle) { mIsFullyChained = toggle; mUnchainedSVs.clear(); }
    public boolean isFullyChained() { return mIsFullyChained; }

    public List<SvClusterData> getUnlinkedSVs() { return mUnchainedSVs; }

    public boolean isSimpleSingleSV() { return isSimpleSVs(); }

    public boolean isSimpleSVs()
    {
        if(mSVs.size() > SMALL_CLUSTER_SIZE)
            return false;

        for(final SvClusterData var : mSVs)
        {
            if(!var.isSimpleType())
                return false;
        }

        return true;
    }

    public final List<SvLinkedPair> getLinkedPairs() { return mLinkedPairs; }
    public final List<SvLinkedPair> getInferredLinkedPairs() { return mInferredLinkedPairs; }
    public final List<SvLinkedPair> getAssemblyLinkedPairs() { return mAssemblyLinkedPairs; }
    public final List<SvClusterData> getSpanningSVs() { return mSpanningSVs; }
    public void setLinkedPairs(final List<SvLinkedPair> pairs) { mLinkedPairs = pairs; }
    public void setInferredLinkedPairs(final List<SvLinkedPair> pairs) { mInferredLinkedPairs = pairs; }
    public void setAssemblyLinkedPairs(final List<SvLinkedPair> pairs) { mAssemblyLinkedPairs = pairs; }
    public void setSpanningSVs(final List<SvClusterData> svList) { mSpanningSVs = svList; }

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
        for(final SvClusterData var : cluster.getSVs())
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

    private void updateClusterDetails()
    {
        if(!mRequiresRecalc)
            return;

        mConsistencyCount = mUtils.calcConsistency(mSVs);

        // for now, just link to this
        mIsConsistent = (mConsistencyCount == 0);

        mDesc = getClusterTypesAsString();

        setMinMaxCopyNumber();

        mRequiresRecalc = false;
    }

    public void logDetails()
    {
        LOGGER.debug("cluster({}) svCount({}) desc({}) armCount({}) consistent({} count={}) chains({}) {}",
                getId(), getCount(), getDesc(), getChromosomalArmCount(), isConsistent(), getConsistencyCount(),
                mChains.size(), isSimpleSingleSV() ? "simple" : (mIsFullyChained ? "full-chained" : ""));
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

        for(final SvClusterData var : mSVs)
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
        for (final SvClusterData var : mSVs)
        {
            if(var.type() == BND && (var.isLineElement(true) || var.isLineElement(false)))
                return true;
        }

        return false;
    }

    private void setMinMaxCopyNumber()
    {
        // first establish the lowest copy number change
        mMinCopyNumber = -1;
        mMaxCopyNumber = -1;

        for (final SvClusterData var : mSVs)
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
        for (final SvClusterData var : mSVs)
        {
            if(var.type() == StructuralVariantType.INS)
                continue;

            if(var.isReplicatedSv())
                continue;

            addVariantToUniqueBreakends(var);
        }

        // cache this against the SV
        for (SvClusterData var : mSVs)
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

    private void addVariantToUniqueBreakends(final SvClusterData var)
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

    public final SvLinkedPair getLinkedPair(final SvClusterData var, boolean useStart)
    {
        return findLinkedPair(mLinkedPairs, var, useStart);
    }

    public final SvChain findChain(final SvClusterData var)
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

    public static final SvCluster findCluster(final SvClusterData var, final List<SvCluster> clusters)
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
