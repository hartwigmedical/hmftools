package com.hartwig.hmftools.svannotation.analysis;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.svannotation.analysis.SvClusterData;
import com.hartwig.hmftools.svannotation.analysis.SvFootprint;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import java.util.List;

public class SvCluster
{
    final SvUtilities mUtils;

    private int mClusterId;
    private int mNextFootprintId;

    private List<SvClusterData> mClusteredSVs;
    private List<SvFootprint> mFootprints;
    private List<SvCluster> mSubClusters;

    private SvClusterData mSpanningSV; // the SV whose positions define this group

    // private List<SvClusterData> mCrossingSVs; // SVs which span more than a single sub-cluster

    // links to a parent cluster if applicable
    //private SvCluster mParentClusterStart;
    //private SvCluster mParentClusterEnd;

    private static final Logger LOGGER = LogManager.getLogger(SvCluster.class);

    public SvCluster(final int clusterId, final SvUtilities utils)
    {
        mNextFootprintId = 0;
        mClusterId = clusterId;
        mUtils = utils;
        mClusteredSVs = Lists.newArrayList();
        mFootprints = Lists.newArrayList();
        mSubClusters = Lists.newArrayList();
        // mSpanningSVs = Lists.newArrayList();

        mSpanningSV = null;
        //mParentClusterStart = null;
        //mParentClusterEnd = null;
    }

    public int getClusterId() { return mClusterId; }

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

    private int getNextFootprintId() { return mNextFootprintId++; }

    public List<SvClusterData> getSVs() { return mClusteredSVs; }
    public List<SvFootprint> getFootprints() { return mFootprints; }
    public List<SvCluster> getSubClusters() { return mSubClusters; }

    public void addVariant(final SvClusterData variant)
    {
        mClusteredSVs.add(variant);

        addToFootprint(variant);
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

        else
        {
            for(SvCluster subCluster : mSubClusters)
            {
                SvCluster matchingCluster = subCluster.findClosestCluster(chromosome, position);
                if(matchingCluster != null)
                    return matchingCluster;
            }

            // no more precise match in the sub-clusters than this cluster
            return this;
        }
    }

    private void addToFootprint(final SvClusterData variant)
    {
        // also add to relevant footprint
        boolean footprintFound = false;
        for(SvFootprint footprint : mFootprints)
        {
            for(final SvClusterData otherVar : footprint.getSVs())
            {
                if(mUtils.areVariantsLinked(variant, otherVar))
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
}
