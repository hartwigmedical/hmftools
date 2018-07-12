package com.hartwig.hmftools.svanalysis.analysis;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.svanalysis.types.SvClusterData;
import com.hartwig.hmftools.svanalysis.types.SvPathData;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class SvPathFinder {


    final SvClusteringConfig mConfig;
    final SvUtilities mUtils;

    List<SvClusterData> mAllVariants;

    HashMap<String, SvPathData> mLinkedSvPaths; // keyed by chromosomal arm if self-contained
    HashMap<String, List<SvClusterData>> mChrArmVariants;

    private static final Logger LOGGER = LogManager.getLogger(ClusterAnalyser.class);

    public static String MIXED_LINKED_PATH = "Mixed";

    public SvPathFinder(final SvClusteringConfig config, final SvUtilities utils)
    {
        mConfig = config;
        mUtils = utils;
        mAllVariants = Lists.newArrayList();
        mLinkedSvPaths = new HashMap();
        mChrArmVariants = new HashMap();
    }

    public void findPaths(final String sampleId, final List<SvClusterData> allVariants)
    {
        mAllVariants = allVariants;

        // assigning all SVs to arms makes the search for paths faster
        assignToChrArms();

        for(Map.Entry<String, List<SvClusterData>> entry : mChrArmVariants.entrySet())
        {
            List<SvClusterData> chrArmVariants = entry.getValue();




            LOGGER.debug("chrArm({}) has {} variants", entry.getKey(), entry.getValue().size());
        }


//        List<SvClusterData> unassignedVariants = Lists.newArrayList(allVariants);
//
//        SvPathData interArmPath = null;
//
//
//        // assign each variant once to a cluster using proximity as a test
//        int currentIndex = 0;
//
//        while(currentIndex < unassignedVariants.size()) {
//
//            SvClusterData currentVar = unassignedVariants.get(currentIndex);
//
//            String chrArmStart = mUtils.getVariantChrArm(currentVar,true);
//            String chrArmEnd = mUtils.getVariantChrArm(currentVar,false);
//
//            SvPathData pathData = null;
//
//            if(chrArmStart.equals(chrArmEnd))
//            {
//                pathData = mLinkedSvPaths.get(chrArmEnd);
//
//                if(pathData == null) {
//                    pathData = new SvPathData(getNextPathId(), chrArmEnd, mUtils);
//                    mLinkedSvPaths.put(chrArmStart, pathData);
//                }
//            }
//            else
//            {
//                if(interArmPath == null)
//                {
//                    interArmPath = new SvPathData(getNextPathId(), MIXED_LINKED_PATH, mUtils);
//                    mLinkedSvPaths.put(MIXED_LINKED_PATH, interArmPath);
//                }
//
//                pathData = interArmPath;
//            }
//
//            pathData.addVariant(currentVar);
//
////            LOGGER.debug("creating new path({}) with next variant({}:{})", newCluster.getId(), currentIndex, currentVar.id());
//
//            unassignedVariants.remove(currentIndex); // index will remain the same and so point to the next item
//
//            // findLinkedVariants(pathData, currentVar, currentIndex, unassignedVariants);
//
//            LOGGER.debug("linkedPathCount({}) remainingVariants({})", mLinkedSvPaths.size(), unassignedVariants.size());
//        }

    }

    private void assignToChrArms()
    {
        // first assign SVs to chromosomal arms, assigning anything cross-arm to its own collection
        // further more, if a chr-arm has a BE in it, move all SVs in that chr-arm to a '_M' (ie mixed) instance of it
        // to separate the search for independent vs cross-arm paths

        for(SvClusterData var : mAllVariants)
        {
            String chrArmStart = mUtils.getVariantChrArm(var,true);
            String chrArmEnd = mUtils.getVariantChrArm(var,false);

            List<SvClusterData> chrArmSvList = null;
            List<SvClusterData> interArmList = null;

            if(chrArmStart.equals(chrArmEnd))
            {
                // first search for a Mixed chr-arm
                String mixedChrArm = chrArmEnd + "_M";
                chrArmSvList = mChrArmVariants.get(mixedChrArm);

                if(chrArmSvList != null) {
                    // a mixed chr-arm already exists, so use this
                }
                else {

                    chrArmSvList = mChrArmVariants.get(chrArmEnd);

                    if (chrArmSvList == null) {
                        chrArmSvList = Lists.newArrayList();
                        mChrArmVariants.put(chrArmStart, chrArmSvList);
                    }
                }
            }
            else
            {
                if(interArmList == null)
                {
                    interArmList = Lists.newArrayList();
                    mChrArmVariants.put(MIXED_LINKED_PATH, interArmList);
                }

                chrArmSvList = interArmList;
            }

            chrArmSvList.add(var);
        }

        // log all chr-arms and their counts
        for(Map.Entry<String, List<SvClusterData>> entry : mChrArmVariants.entrySet())
        {
            LOGGER.debug("chrArm({}) has {} variants", entry.getKey(), entry.getValue().size());
        }
    }

    /*
    private void findLinkedVariants(SvPathData pathData, SvClusterData currentVar, List<SvClusterData> unassignedVariants)
    {
        // since the SVs are order from lowest chromosome and position, the first SV in a new path starts with the start position
        // and so looks for a linked end using its end position
        boolean curUnlinkedBEIsStart = false;

        while (!unassignedVariants.isEmpty()) {



            // compare with all other SVs in this cluster
            boolean matched = false;
            for (SvClusterData otherVar : cluster.getSVs())
            {
                // test each possible linkage
                if (!mUtils.areVariantsLinkedByDistance(currentVar, otherVar))
                {
                    //LOGGER.debug("non-linked SVs: v1({}) and v2({})", currentVar.posId(), otherVar.posId());
                    continue;
                }

                pathData.addVariant(currentVar);
                LOGGER.debug("pathData({}) add matched variant({}), totalCount({})",
                        pathData.getId(), currentVar.id(), pathData.getLinkedVariants().size());

                matched = true;
                break;
            }

            if(matched)
            {
                unassignedVariants.remove(currentIndex);

                // as soon as a new SV is added to this cluster, need to start checking from the beginning again
                currentIndex = 0;
            }
            else
            {
                ++currentIndex;
            }
        }

    }
    */

    private int getNextPathId()
    {
        return mLinkedSvPaths.size();
    }



}
