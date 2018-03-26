package com.hartwig.hmftools.svannotation.analysis;

import static com.hartwig.hmftools.svannotation.FragileSiteAnnotator.NO_FS;
import static com.hartwig.hmftools.svannotation.LineElementAnnotator.NO_LINE_ELEMENT;
import static com.hartwig.hmftools.svannotation.analysis.SvUtilities.DOUBLE_STRANDED_BREAK;
import static com.hartwig.hmftools.svannotation.analysis.SvUtilities.REPLICATION_EVENT;

import com.google.common.collect.Lists;

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
    private int mDSBCount ;
    private int mTempInsertCount;
    private String mDSBLengths;
    private String mTempInsertLengths;
    private int mChromosomeArmCount;

    private List<SvClusterData> mClusteredSVs;
    private List<SvBreakend> mUniqueBreakends;
    private List<SvFootprint> mFootprints;
    private List<SvCluster> mSubClusters;

    private SvClusterData mSpanningSV; // the SV whose positions define this group

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
        mAnnotationList= Lists.newArrayList();
        mDSBCount = 0;
        mDSBLengths = "";
        mTempInsertCount = 0;
        mTempInsertLengths = "";
        mChromosomeArmCount = 0;

        // subclustering and linkages
        mFootprints = Lists.newArrayList();
        mSubClusters = Lists.newArrayList();

        mSpanningSV = null;
    }

    public int getClusterId() { return mClusterId; }

    private int getNextFootprintId() { return mNextFootprintId++; }

    public List<SvClusterData> getSVs() { return mClusteredSVs; }
    public List<SvFootprint> getFootprints() { return mFootprints; }
    public List<SvCluster> getSubClusters() { return mSubClusters; }

    public void addVariant(final SvClusterData variant)
    {
        mClusteredSVs.add(variant);

        // addToFootprint(variant);
    }

    public final String getDesc() { return mDesc; }
    public final void setDesc(final String desc) { mDesc = desc; }
    public final List<String> getAnnotationList() { return mAnnotationList; }
    public final void addAnnotation(final String annotation)
    {
        if(mAnnotationList.contains(annotation))
            return;

        mAnnotationList.add(annotation);
    }

    public String getAnnotations() { return mAnnotationList.stream ().collect (Collectors.joining (";")); }

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
        LOGGER.debug("cluster({}) has {} unique arms", mClusterId, mChromosomeArmCount);
    }

    public final boolean isReplicationEvent()
    {
        return mAnnotationList.contains(REPLICATION_EVENT);
    }

    public final int getDSBCount() { return mDSBCount; }
    public void setDSBCount(int count) { mDSBCount = count; }

    public final String getDSBLengths() { return mDSBLengths; }
    public final void setDSBLengths(final String lens) { mDSBLengths = lens; }

    public final int getTICount() { return mTempInsertCount; }
    public void setTICount(int count) { mTempInsertCount = count; }

    public final String getTempInsertLengths() { return mTempInsertLengths; }
    public final void setTempInsertLengths(final String lens) { mTempInsertLengths = lens; }

    public final boolean isDSBEvent()
    {
        return mAnnotationList.contains(DOUBLE_STRANDED_BREAK);
    }

    public final String getClusterTypesAsString()
    {
        if(mClusteredSVs.size() == 1)
        {
            return mClusteredSVs.get(0).typeStr();
        }
        else if(mClusteredSVs.size() == 2)
        {
            return mClusteredSVs.get(0).typeStr() + "_" + mClusteredSVs.get(1).typeStr();
        }

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
    }

    private void addVariantToUniqueBreakends(final SvClusterData var)
    {
        for(int i = 0; i < 2; ++i)
        {
            boolean useStart = (i == 0);

            boolean found = false;
            for(SvBreakend breakend : mUniqueBreakends) {
                if (variantMatchesBreakend(var, breakend, useStart))
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

    private boolean variantMatchesBreakend(final SvClusterData var, final SvBreakend breakend, boolean useStart)
    {
        return breakend.chromosome().equals(var.chromosome(useStart))
            && breakend.position() == var.position(useStart)
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

    public boolean isSvDuplicateBE(final SvClusterData var, boolean useStart)
    {
        for(final SvBreakend breakend : mUniqueBreakends)
        {
            if(variantMatchesBreakend(var, breakend, useStart))
                return breakend.getCount() > 1;
        }

        return false;
    }


    public void addSubCluster(SvCluster subCluster)
    {
        mSubClusters.add(subCluster);
    }

    public SvCluster findClosestCluster(String chromosome, long position)
    {
        if(mSpanningSV == null)
            return null;

        if(!mUtils.isWithin(mSpanningSV, chromosome, position))
        {
            return null;
        }

        if(mSubClusters.isEmpty())
            return this;

        for(SvCluster subCluster : mSubClusters)
        {
            SvCluster matchingCluster = subCluster.findClosestCluster(chromosome, position);
            if(matchingCluster != null)
                return matchingCluster;
        }

        // no more precise match in the sub-clusters than this cluster
        return this;
    }

    private void addToFootprint(final SvClusterData variant)
    {
        // also add to relevant footprint
        boolean footprintFound = false;
        for(SvFootprint footprint : mFootprints)
        {
            for(final SvClusterData otherVar : footprint.getSVs())
            {
                if(mUtils.areVariantsLinkedByDistance(variant, otherVar))
                {
                    footprint.addVariant(variant);
                    LOGGER.debug("cluster({}) footprint({}) added variant({}), count({})", mClusterId, footprint.getFootprintId(), variant.id(), footprint.getSVs().size());
                    footprintFound = true;
                    break;

                }
            }
        }

        if(!footprintFound)
        {
            // add a new one
            SvFootprint newFootprint = new SvFootprint(getNextFootprintId());

            LOGGER.debug("cluster({}) new footprint({}) added variant({}),", mClusterId, newFootprint.getFootprintId(), variant.id());
            newFootprint.addVariant(variant);

            mFootprints.add(newFootprint);
        }
    }

    public void setSpanningSV(SvClusterData spanningSV)
    {
        addVariant(spanningSV);
        mSpanningSV = spanningSV;
        mSpanningSV.setStartCluster(this);
        mSpanningSV.setEndCluster(this);

        LOGGER.debug("cluster({}) set spanning variant({}) with subSVsCount({})", mClusterId, spanningSV.id(), spanningSV.getSubSVs().size());

        // and add any of its sub-SVs
        for(SvClusterData subVar : spanningSV.getSubSVs())
        {
            addVariant(subVar);
        }
    }

    public SvClusterData getSpanningSV() { return mSpanningSV; }

    public long getSpan() { return mSpanningSV != null ? mSpanningSV.getSpan() : -1; }

    public boolean isIsolatedSingleSV()
    {
        // has no other overlapping SVs or wholy contained SVs
        return (mClusteredSVs.size() == 1);
    }

}
