package com.hartwig.hmftools.svanalysis.analysis;

import static java.lang.Math.abs;

import static com.hartwig.hmftools.svanalysis.analysis.SvUtilities.PERMITED_DUP_BE_DISTANCE;
import static com.hartwig.hmftools.svanalysis.annotators.FragileSiteAnnotator.NO_FS;
import static com.hartwig.hmftools.svanalysis.annotators.LineElementAnnotator.NO_LINE_ELEMENT;
import static com.hartwig.hmftools.svanalysis.analysis.SvUtilities.DOUBLE_STRANDED_BREAK;
import static com.hartwig.hmftools.svanalysis.analysis.SvUtilities.REPLICATION_EVENT;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.svanalysis.types.SvBreakend;
import com.hartwig.hmftools.svanalysis.types.SvChain;
import com.hartwig.hmftools.svanalysis.types.SvClusterData;
import com.hartwig.hmftools.svanalysis.types.SvFootprint;
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
    private int mNextFootprintId;

    private int mConsistencyCount;
    private boolean mIsConsistent;
    private String mDesc;
    private List<String> mAnnotationList;
    private int mChromosomeArmCount;

    private List<SvClusterData> mClusteredSVs;
    private List<SvBreakend> mUniqueBreakends;
    private List<SvChain> mChains;
    private List<SvLinkedPair> mLinkedPairs;
    private List<SvClusterData> mSpanningSVs;

    private static final Logger LOGGER = LogManager.getLogger(SvCluster.class);

    public SvCluster(final int clusterId, final SvUtilities utils)
    {
        mNextFootprintId = 0;
        mClusterId = clusterId;
        mUtils = utils;
        mClusteredSVs = Lists.newArrayList();
        mUniqueBreakends = Lists.newArrayList();

        // annotation info
        mConsistencyCount = 0;
        mIsConsistent = false;
        mDesc = "";
        mAnnotationList = Lists.newArrayList();
        mChromosomeArmCount = 0;

        // chain data
        mLinkedPairs = Lists.newArrayList();
        mSpanningSVs = Lists.newArrayList();
        mChains = Lists.newArrayList();
    }

    public int getId() { return mClusterId; }

    public int getCount() { return mClusteredSVs.size(); }
    public List<SvClusterData> getSVs() { return mClusteredSVs; }

    public void addVariant(final SvClusterData variant)
    {
        mClusteredSVs.add(variant);
    }

    public List<SvChain> getChains() { return mChains; }

    public void addChain(SvChain chain)
    {
        mChains.add(chain);
    }

    public final List<SvLinkedPair> getLinkedPairs() { return mLinkedPairs; }
    public final List<SvClusterData> getSpanningSVs() { return mSpanningSVs; }
    public void setLinkedPairs(final List<SvLinkedPair> pairs) { mLinkedPairs = pairs; }
    public void setSpanningSVs(final List<SvClusterData> svList) { mSpanningSVs = svList; }

    public final String getDesc() { return mDesc; }
    public final void setDesc(final String desc) { mDesc = desc; }

    public void setIsConsistent(boolean toggle) { mIsConsistent = toggle; }

    public boolean isConsistent() { return mIsConsistent; }

    public int getConsistencyCount()
    {
        return mConsistencyCount;
    }

    public void setConsistencyCount()
    {
        mConsistencyCount = mUtils.calcConsistency(mClusteredSVs);

        // for now, just link to this
        mIsConsistent = (mConsistencyCount == 0);
    }

    public int getChromosomalArmCount() { return mChromosomeArmCount; }

    public void setChromosomalArmCount()
    {
        // unique arm count
        List<String> chrArmlist = Lists.newArrayList();

        for(final SvClusterData var : mClusteredSVs)
        {
            String chrArmStart = var.chromosome(true) + "_" + var.arm(true);
            String chrArmEnd = var.chromosome(false) + "_" + var.arm(false);

            if(!chrArmlist.contains(chrArmStart))
                chrArmlist.add(chrArmStart);

            if(!chrArmlist.contains(chrArmEnd))
                chrArmlist.add(chrArmEnd);
        }

        mChromosomeArmCount = chrArmlist.size();
    }

    public final String getClusterTypesAsString()
    {
        if(mClusteredSVs.size() == 1)
        {
            return mClusteredSVs.get(0).typeStr();
        }

        // the following map-based naming convention leads
        // to a predictable ordering of types: INV, CRS, BND, DEL and DUP
        String clusterTypeStr = "";
        Map<String, Integer> typeMap = new HashMap<>();

        for(final SvClusterData var : mClusteredSVs) {

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

    public int getFragileSiteCount() {
        int count = 0;
        for (final SvClusterData var : mClusteredSVs) {
            if(var.isStartFragileSite() != NO_FS || var.isEndFragileSite() != NO_FS)
                ++count;
        }

        return count;
    }

    public int getLineElementCount() {
        int count = 0;
        for (final SvClusterData var : mClusteredSVs) {
            if(var.isStartLineElement() != NO_LINE_ELEMENT || var.isEndLineElement() != NO_LINE_ELEMENT)
                ++count;
        }

        return count;
    }

    public void setUniqueBreakends()
    {
        // group any matching BEs (same position and orientation)
        for (final SvClusterData var : mClusteredSVs) {
            addVariantToUniqueBreakends(var);
        }

        // cache this against the SV
        for (SvClusterData var : mClusteredSVs) {

            for(final SvBreakend breakend : mUniqueBreakends)
            {
                if(variantMatchesBreakend(var, breakend, true, PERMITED_DUP_BE_DISTANCE) && breakend.getCount() > 1)
                    var.setIsDupBEStart(true);

                if(variantMatchesBreakend(var, breakend, false, PERMITED_DUP_BE_DISTANCE) && breakend.getCount() > 1)
                    var.setIsDupBEEnd(true);
            }
        }
    }

    private void addVariantToUniqueBreakends(final SvClusterData var)
    {
        for(int i = 0; i < 2; ++i)
        {
            boolean useStart = (i == 0);

            boolean found = false;
            for(SvBreakend breakend : mUniqueBreakends) {
                if (variantMatchesBreakend(var, breakend, useStart, PERMITED_DUP_BE_DISTANCE))
                {
                    breakend.addToCount(1);
                    found = true;
                    break;
                }
            }

            if(!found) {
                // add a new entry
                mUniqueBreakends.add(new SvBreakend(var.chromosome(useStart), var.position(useStart), var.orientation(useStart)));
            }
        }
    }

    public boolean variantMatchesBreakend(final SvClusterData var, final SvBreakend breakend, boolean useStart, int permittedDist)
    {
        return breakend.chromosome().equals(var.chromosome(useStart))
            && abs(breakend.position() - var.position(useStart)) <= permittedDist
            && breakend.orientation() == var.orientation(useStart);
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
        for(final SvLinkedPair pair : mLinkedPairs)
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
