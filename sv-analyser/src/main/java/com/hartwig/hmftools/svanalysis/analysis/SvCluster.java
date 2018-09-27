package com.hartwig.hmftools.svanalysis.analysis;

import static java.lang.Math.abs;

import static com.hartwig.hmftools.svanalysis.analysis.SvUtilities.PERMITED_DUP_BE_DISTANCE;
import static com.hartwig.hmftools.svanalysis.analysis.SvUtilities.variantMatchesBreakend;
import static com.hartwig.hmftools.svanalysis.annotators.LineElementAnnotator.NO_LINE_ELEMENT;
import static com.hartwig.hmftools.svanalysis.types.SvClusterData.SVI_END;
import static com.hartwig.hmftools.svanalysis.types.SvClusterData.SVI_START;
import static com.hartwig.hmftools.svanalysis.types.SvClusterData.isStart;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantType;
import com.hartwig.hmftools.svanalysis.types.SvArmGroup;
import com.hartwig.hmftools.svanalysis.types.SvBreakend;
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
    private List<String> mAnnotationList;

    private List<SvClusterData> mSVs; // SVs within close promximity of each other
    private List<SvBreakend> mUniqueBreakends; // for duplicate BE searches
    private List<SvChain> mChains; // pairs of SVs linked into chains
    private List<SvLinkedPair> mLinkedPairs; // forming a TI or DB
    private List<SvClusterData> mSpanningSVs; // having 2 duplicate (matching) BEs
    private List<SvArmGroup> mArmGroups;

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
        mAnnotationList = Lists.newArrayList();

        // chain data
        mLinkedPairs = Lists.newArrayList();
        mSpanningSVs = Lists.newArrayList();
        mChains = Lists.newArrayList();
    }

    public int getId() { return mClusterId; }

    public int getCount() { return mSVs.size(); }

    public final String getDesc() { return mDesc; }
    public final void setDesc(final String desc) { mDesc = desc; }

    public List<SvClusterData> getSVs() { return mSVs; }

    public void addVariant(final SvClusterData var)
    {
        mSVs.add(var);

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

    public List<SvChain> getChains() { return mChains; }
    public void addChain(SvChain chain) { mChains.add(chain); }

    public final List<SvLinkedPair> getLinkedPairs() { return mLinkedPairs; }
    public final List<SvClusterData> getSpanningSVs() { return mSpanningSVs; }
    public void setLinkedPairs(final List<SvLinkedPair> pairs) { mLinkedPairs = pairs; }
    public void setSpanningSVs(final List<SvClusterData> svList) { mSpanningSVs = svList; }

    public boolean isConsistent() { return mIsConsistent; }

    public int getConsistencyCount()
    {
        return mConsistencyCount;
    }
    public double getConsistencyScore() { return abs(mConsistencyCount/(double)getCount()); }

    public void setConsistencyCount()
    {
        mConsistencyCount = mUtils.calcConsistency(mSVs);

        // for now, just link to this
        mIsConsistent = (mConsistencyCount == 0);
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

    public int getLineElementCount()
    {
        int count = 0;
        for (final SvClusterData var : mSVs)
        {
            if(!var.isStartLineElement().equals(NO_LINE_ELEMENT) || !var.isEndLineElement().equals(NO_LINE_ELEMENT))
                ++count;
        }

        return count;
    }

    public void setUniqueBreakends()
    {
        // group any matching BEs (same position and orientation)
        for (final SvClusterData var : mSVs)
        {
            if(var.type() == StructuralVariantType.INS)
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
        for(int i = 0; i < 2; ++i)
        {
            boolean useStart = (i == 0);

            if(var.type() == StructuralVariantType.INV && var.position(false) - var.position(true) <= PERMITED_DUP_BE_DISTANCE)
            {
                // avoid setting both end of an INV as duplicate if they match
                if(i == 1)
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

    public int getDuplicateBECount()
    {
        int count = 0;
        for(final SvBreakend breakend : mUniqueBreakends)
        {
            if(breakend.getCount() > 1)
                count += breakend.getCount();
        }

        return count;
    }

    public int getDuplicateBESiteCount()
    {
        int count = 0;
        for(final SvBreakend breakend : mUniqueBreakends)
        {
            if(breakend.getCount() > 1)
                ++count;
        }

        return count;
    }

    public final SvLinkedPair findLinkedPair(final SvClusterData var, boolean useStart)
    {
        return findLinkedPair(mLinkedPairs, var, useStart);
    }

    public static final SvLinkedPair findLinkedPair(final List<SvLinkedPair> linkedPairs, final SvClusterData var, boolean useStart)
    {
        for(final SvLinkedPair pair : linkedPairs)
        {
            if(var.equals(pair.first()) && useStart == pair.firstLinkOnStart())
                return pair;

            if(var.equals(pair.second()) && useStart == pair.secondLinkOnStart())
                return pair;
        }

        return null;
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

    public final List<String> getAnnotationList() { return mAnnotationList; }
    public final void addAnnotation(final String annotation)
    {
        if(mAnnotationList.contains(annotation))
            return;

        mAnnotationList.add(annotation);
    }

    public String getAnnotations() { return mAnnotationList.stream ().collect (Collectors.joining (";")); }

}
