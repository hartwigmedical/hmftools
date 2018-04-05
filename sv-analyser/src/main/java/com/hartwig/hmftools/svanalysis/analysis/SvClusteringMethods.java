package com.hartwig.hmftools.svanalysis.analysis;

import static com.hartwig.hmftools.svanalysis.analysis.SvUtilities.SV_GROUP_ENCLOSED;
import static com.hartwig.hmftools.svanalysis.analysis.SvUtilities.SV_GROUP_ENCLOSING;
import static com.hartwig.hmftools.svanalysis.analysis.SvUtilities.SV_GROUP_NEIGHBOURS;
import static com.hartwig.hmftools.svanalysis.analysis.SvUtilities.SV_GROUP_OVERLAP;

import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.common.variant.structural.StructuralVariantType;
import com.hartwig.hmftools.svanalysis.types.SvClusterData;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class SvClusteringMethods {

    private static final Logger LOGGER = LogManager.getLogger(SvClusteringMethods.class);

    private final SvUtilities mUtils;

    private Map<String, Integer> mChrArmSvCount;
    private Map<String, Double> mChrArmSvExpected;
    private Map<String, Double> mChrArmSvRate;
    private double mMedianChrArmRate;

    private static double REF_BASE_LENGTH = 10000000.0;

    public SvClusteringMethods(final SvUtilities clusteringUtils)
    {
        mUtils = clusteringUtils;
        mChrArmSvCount = new HashMap();
        mChrArmSvExpected = new HashMap();
        mChrArmSvRate = new HashMap();
        mMedianChrArmRate = 0;
    }

    public void findLinkedSVsByDistance(SvCluster cluster, List<SvClusterData> unassignedVariants) {
        // look for any other SVs which form part of this cluster based on proximity
        int currentIndex = 0;

        while (currentIndex < unassignedVariants.size()) {
            SvClusterData currentVar = unassignedVariants.get(currentIndex);

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

                cluster.addVariant(currentVar);
                LOGGER.debug("cluster({}) add matched variant({}), totalCount({})", cluster.getClusterId(), currentVar.id(), cluster.getSVs().size());

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

    public void setClosestSVData(List<SvClusterData> allVariants, List<SvCluster> clusters) {

        for (SvCluster cluster : clusters) {

            // logically the closest SV will already be in the same cluster, otherwise use the full set
            boolean useClusterSVs = cluster.getSVs().size() > 1;

            List<SvClusterData> variants = useClusterSVs ? cluster.getSVs() : allVariants;

            for (int index1 = 0 ; index1 < variants.size(); ++index1) {

                SvClusterData var1 = variants.get(index1);

                int shortestLength = var1.getNearestSVLength(); // use a starting value if set from a previous linked variant
                String shortestLinkageType = var1.getNearestSVLinkType();

                // start at the next variant to avoid rechecking pairs
                for (int index2 = index1+1 ; index2 < variants.size(); ++index2) {

                    SvClusterData var2 = variants.get(index2);

                    int proximity = mUtils.getShortestProximity(var1, var2);

                    if(proximity < 0)
                        continue;

                    String linkageType = "";

                    // apply some specific tests
                    if (mUtils.isWithin(var1, var2)) {
                        linkageType = SV_GROUP_ENCLOSING;
                    } else if (mUtils.isWithin(var2, var1)) {
                        linkageType = SV_GROUP_ENCLOSED;
                    } else if (mUtils.isOverlapping(var1, var2)) {
                        linkageType = SV_GROUP_OVERLAP;
                    } else {
                        linkageType = SV_GROUP_NEIGHBOURS;
                    }

                    // some validity checks
                    // NOTE: this may then lead to an inconsistency with clustering by base distance alone
                    if(mUtils.areTypePair(var1, var2, StructuralVariantType.DEL, StructuralVariantType.DEL)
                    && linkageType != SV_GROUP_NEIGHBOURS)
                    {
                        continue;
                    }

                    if (shortestLength == -1 || proximity < shortestLength) {
                        // record this shortest linkage info
                        shortestLength = proximity;
                        shortestLinkageType = linkageType;
                    }

                    // update the second variant if it's also now nearer
                    if(var2.getNearestSVLength() == -1 || proximity < var2.getNearestSVLength())
                    {
                        var2.setNearestSVLength(proximity);

                        if(shortestLinkageType == SV_GROUP_ENCLOSING)
                            var2.setNearestSVLinkType(SV_GROUP_ENCLOSED);
                        else if(shortestLinkageType == SV_GROUP_ENCLOSED)
                            var2.setNearestSVLinkType(SV_GROUP_ENCLOSING);
                        else
                            var2.setNearestSVLinkType(shortestLinkageType);
                    }
                }

                var1.setNearestSVLength(shortestLength);
                var1.setNearestSVLinkType(shortestLinkageType);

                LOGGER.debug("var({}) shortestLength({}) and type({})", var1.id(), shortestLength, shortestLinkageType);
            }
        }
    }

    public void setClosestLinkedSVData(List<SvClusterData> allVariants) {

        for (int index1 = 0 ; index1 < allVariants.size(); ++index1) {

            SvClusterData var1 = allVariants.get(index1);

            int shortestTILength = var1.getNearestTILength(); // use a starting value if set from a previous linked variant
            int shortestDBLength = var1.getNearestDBLength();

            for (int index2 = index1+1 ; index2 < allVariants.size(); ++index2) {

                SvClusterData var2 = allVariants.get(index2);

                int pairTILength = -1;
                int pairDBLength = -1;

                // compare each BE for this pair
                for(int i = 0; i < 2; ++i) {

                    boolean v1Start = (i == 0);

                    for (int j = 0; j < 2; ++j) {

                        boolean v2Start = (j == 0);

                        if (mUtils.areLinkedSection(var1, var2, v1Start, v2Start)) {

                            // test for templated insertion
                            int tiLength = mUtils.getProximity(var1, var2, v1Start, v2Start);

                            if (tiLength >= 0 && (pairTILength < 0 || tiLength < pairTILength)) {
                                // new shortest TI
                                pairTILength = tiLength;
                            }
                        }

                        if(mUtils.areSectionBreak(var1, var2, v1Start, v2Start))
                        {
                            // then for a deletion bridge
                            int dbLength = mUtils.getProximity(var1, var2, v1Start, v2Start);

                            if (dbLength >= 0 && (pairDBLength < 0 || dbLength < pairDBLength)) {
                                pairDBLength = dbLength;
                            }
                        }
                    }
                }

                if (pairTILength >= 0 && (shortestTILength < 0 || pairTILength < shortestTILength)) {
                    shortestTILength = pairTILength;
                }

                if (pairDBLength >= 0 && (shortestDBLength < 0 || pairDBLength < shortestDBLength)) {
                    shortestDBLength = pairDBLength;
                }

                // update the second variant if it's also now nearer (this is an optimisation)
                if(pairTILength >= 0 && (var2.getNearestTILength() == -1 || pairTILength < var2.getNearestTILength()))
                {
                    var2.setNearestTILength(pairTILength);
                }

                if(pairDBLength >= 0 && (var2.getNearestDBLength() == -1 || pairDBLength < var2.getNearestDBLength()))
                {
                    var2.setNearestDBLength(pairDBLength);
                }
            }

            var1.setNearestTILength(shortestTILength);
            var1.setNearestDBLength(shortestDBLength);

            LOGGER.debug("var({}) shortestTILength({}) shortestDBLength({})", var1.id(), shortestTILength, shortestDBLength);
        }
    }

    public void setChromosomalArmStats(final List<SvClusterData> allVariants)
    {
        mChrArmSvCount.clear();
        mChrArmSvExpected.clear();
        mChrArmSvRate.clear();
        mMedianChrArmRate = 0;

        // form a map of unique arm to SV count
        for(final SvClusterData var : allVariants)
        {
            String chrArmStart = mUtils.getVariantChrArm(var,true);
            String chrArmEnd = mUtils.getVariantChrArm(var,false);

            if(mChrArmSvCount.containsKey(chrArmStart))
            {
                mChrArmSvCount.replace(chrArmStart, mChrArmSvCount.get(chrArmStart) + 1);
            }
            else
            {
                mChrArmSvCount.put(chrArmStart, 1);
            }

            if(mChrArmSvCount.containsKey(chrArmEnd))
            {
                mChrArmSvCount.replace(chrArmEnd, mChrArmSvCount.get(chrArmEnd) + 1);
            }
            else
            {
                mChrArmSvCount.put(chrArmEnd, 1);
            }
        }

        // now normalise these counts by expressing as an expected

        // now determine the background rate by taking the median value from amongst the arms
        // factoring in the arms which have no Q (14-16, 21-22) and excluding the X & Ys
        for(Map.Entry<String, Integer> entry : mChrArmSvCount.entrySet()) {

            final String chrArm = entry.getKey();
            final String chromosome = mUtils.getChrFromChrArm(chrArm);
            final String arm = mUtils.getArmFromChrArm(chrArm);

            long chrArmLength = mUtils.getChromosomalArmLength(chromosome, arm);
            int svCount = entry.getValue();
            double ratePerLength = svCount / (chrArmLength / REF_BASE_LENGTH); // the factor isn't important

            mChrArmSvRate.put(chrArm, ratePerLength);
            LOGGER.debug("chrArm({}) ratePerMill({}) from count({}) length({})", chrArm, ratePerLength, svCount, chrArmLength);
        }

        mChrArmSvRate = sortByValue(mChrArmSvRate, false);

        mMedianChrArmRate = 0;
        int chrArmIndex = 0;
        for(Map.Entry<String, Double> entry : mChrArmSvRate.entrySet())
        {
            LOGGER.debug("chrArm({}: {}) svRate({})", chrArmIndex, entry.getKey(), entry.getValue());

            if(chrArmIndex == 20)
                mMedianChrArmRate = entry.getValue();

           ++chrArmIndex;
        }

        LOGGER.debug("median SV rate({})", mMedianChrArmRate);

        // now create another map of expected SV count per arm using the median rate
        for(Map.Entry<String, Double> entry : mChrArmSvRate.entrySet()) {

            final String chrArm = entry.getKey();
            final String chromosome = mUtils.getChrFromChrArm(chrArm);
            final String arm = mUtils.getArmFromChrArm(chrArm);

            long chrArmLength = mUtils.getChromosomalArmLength(chromosome, arm);
            double expectedSvCount = (int)Math.round((chrArmLength / REF_BASE_LENGTH) * mMedianChrArmRate);
            LOGGER.debug("chrArm({}) expectedSvCount({}) vs actual({})", chrArm, expectedSvCount, mChrArmSvCount.get(chrArm));

            mChrArmSvExpected.put(chrArm, expectedSvCount);
        }
    }

    public String getChrArmData(final SvClusterData var)
    {
        String chrArmStart = mUtils.getVariantChrArm(var,true);
        String chrArmEnd = mUtils.getVariantChrArm(var,false);

        // report Start SV count : Expected SV Count : End SV Count : Expected SV Count
        return String.format("%d,%.2f,%d,%.2f",
                mChrArmSvCount.get(chrArmStart), mChrArmSvExpected.get(chrArmStart),
                mChrArmSvCount.get(chrArmEnd), mChrArmSvExpected.get(chrArmEnd));
    }

    private static Map<String, Double> sortByValue(Map<String, Double> unsortMap, final boolean order)
    {
        List<Map.Entry<String, Double>> list = new LinkedList<Map.Entry<String, Double>>(unsortMap.entrySet());

        Collections.sort(list, new Comparator<Map.Entry<String, Double>>()
        {
            public int compare(Map.Entry<String, Double> o1,
                    Map.Entry<String, Double> o2)
            {
                if (order)
                {
                    return o1.getValue().compareTo(o2.getValue());
                }
                else
                {
                    return o2.getValue().compareTo(o1.getValue());

                }
            }
        });

        // Maintaining insertion order with the help of LinkedList
        Map<String, Double> sortedMap = new LinkedHashMap<String, Double>();
        for (Map.Entry<String, Double> entry : list)
        {
            sortedMap.put(entry.getKey(), entry.getValue());
        }

        return sortedMap;
    }


    // these methods are unused for now..

    /*

    public void runClusteringByLocals() {
        // first find outer-most clusters at CRMS level and identify others which will span CRMS

        // LOGGER.debug("sample({}) isolating inter-chromosomal variants", mSampleId);

        // first remove cross-chromosomal from initial consideration
        // List<SvClusterData> crossChromosomal = findCrossChromosomalSVs();

        LOGGER.debug("sample({}) clustering contained variants", mSampleId);

        createOutermostChromosomalClusters();

        LOGGER.debug("sample({}) assigning sub-clusters", mSampleId);

        // now need to assign and create sub-clusters within these outer-most clusters
        for(SvCluster cluster : mClusters) {
            assignSubClusters(cluster);
        }

        LOGGER.debug("sample({}) matching {} overlapping variants", mSampleId, mUnassignedVariants.size());

        // final create links for the cross-chromosomal SVs and other overlapping SVs
        assignCrossClusterSVs();
    }

    private List<SvClusterData> findCrossChromosomalSVs()
    {
        // first remove cross-chromosomal from initial consideration
        List<SvClusterData> crossChromosomal = Lists.newArrayList();

        int currentIndex = 0;
        while (currentIndex < mUnassignedVariants.size()) {
            SvClusterData currentVar = mUnassignedVariants.get(currentIndex);

            // first skip over cross-CRMS for now
            if (!currentVar.isLocal()) {
                crossChromosomal.add(currentVar);

                LOGGER.debug("cross-CRMS SV: {}", currentVar.posId());

                mUnassignedVariants.remove(currentIndex); // index will point at next
                continue;
            } else {
                ++currentIndex;
            }
        }

        return crossChromosomal;

    }

    private void createOutermostChromosomalClusters()
    {
        int currentIndex = 0;

        List<SvClusterData> subSVsToRemove = Lists.newArrayList();

        while(currentIndex < mUnassignedVariants.size()) {
            SvClusterData currentVar = mUnassignedVariants.get(currentIndex);

            if(!currentVar.isLocal())
            {
                LOGGER.debug("skipping inter-chromosomal SV: {}", currentVar.posId());
                ++currentIndex;
                continue;
            }

            boolean spansOtherSVs = false;
            List<SvClusterData> subSVs = Lists.newArrayList();

            for (SvClusterData otherVar : mUnassignedVariants)
            {
                if(otherVar.id().equals(currentVar.id()))
                {
                    continue;
                }

                if(mClusteringUtils.isLocalOverlap(currentVar, otherVar))
                {
                    LOGGER.debug("local overlap SVs: v1({}) and v2({})", currentVar.posId(), otherVar.posId());
                    spansOtherSVs = true;
                }

                if (mClusteringUtils.isWithin(currentVar, otherVar))
                {
                    if(currentVar.addSubSV(otherVar)) {
                        LOGGER.debug("wholy contained: outer({}) and inner({})", currentVar.posId(), otherVar.posId());
                    }

                }
                else if (mClusteringUtils.isWithin(otherVar, currentVar))
                {
                    if(otherVar.addSubSV(currentVar)) {
                        LOGGER.debug("wholy contained: outer({}) and inner({})", otherVar.posId(), currentVar.posId());
                    }
                }
            }

            if(!currentVar.isSubSV() && (currentVar.hasSubSVs() || !spansOtherSVs))
            {
                LOGGER.debug("adding outer-most CRMS SV: {}", currentVar.posId());

                // make a new cluster
                SvCluster newCluster = new SvCluster(getNextClusterId(), mClusteringUtils);

                newCluster.setSpanningSV(currentVar);
                mClusters.add(newCluster);

                mUnassignedVariants.remove(currentIndex); // index will point at next
            }
            else
            {
                if(currentVar.isSubSV())
                    subSVsToRemove.add(currentVar);

                ++currentIndex;
            }
        }

        for(SvClusterData variant : subSVsToRemove)
        {
            // remove from unassigned
            mUnassignedVariants.remove(variant);
        }
    }

    private void assignSubClusters(SvCluster cluster)
    {
        LOGGER.debug("cluster({}) sub-clustering {} variants", cluster.getClusterId(), cluster.getSVs().size());

        // create sub-clusters for any SVs which have sub SVs
        for(SvClusterData currentVar : cluster.getSVs())
        {
            if(currentVar == cluster.getSpanningSV())
                continue;

            if(!currentVar.hasSubSVs())
            {
                if(currentVar.areClustersSet())
                {
                    // keep with the smallest cluster
                    long clusterSpan = currentVar.getStartCluster().getSpan();
                    long newSpan = cluster.getSpan();

                    if(newSpan == -1 || newSpan > clusterSpan)
                    {
                        continue;
                    }

                }

                LOGGER.debug("variant({}) assigned to cluster({}) span({} -> {})",
                        currentVar.posId(), cluster.getClusterId(), cluster.getSpanningSV().position(true), cluster.getSpanningSV().position(false));

                currentVar.setStartCluster(cluster);
                currentVar.setEndCluster(cluster);
                continue;
            }
            else {
                // also check that this variant is a sub-variant of another variant at this level
                boolean skipClustering = false;
                for(SvClusterData other : cluster.getSVs())
                {
                    if(currentVar.equals(other) || other == cluster.getSpanningSV())
                        continue;

                    if(other.getSubSVs().contains(currentVar))
                    {
                        LOGGER.debug("skipping sub-cluster for variant({}) since part of other({})", currentVar.id(), other.id());
                        skipClustering = true;
                        break;
                    }
                }

                if(skipClustering)
                    continue;
            }

            // create a cluster, set its parent cluster and then call iteratively
            SvCluster newCluster = new SvCluster(getNextClusterId(), mClusteringUtils);

            LOGGER.debug("cluster({}) creating new subCluster({})", cluster.getClusterId(), newCluster.getClusterId());

            // add the spanning SV and its sub-SVs
            newCluster.setSpanningSV(currentVar);
            cluster.addSubCluster(newCluster);

            // and call recursively to create lower-level clusters
            assignSubClusters(newCluster);
        }
    }

    private void assignCrossClusterSVs()
    {
        // for all remaining SVs, including those which cross chromosomes,
        // assign them to the most precise cluster (ie that with the smallest positional range)
        for(SvClusterData currentVar : mUnassignedVariants) {

            for(int i = 0; i < 2; ++i) {

                boolean isStart = (i == 0);
                String isStartStr = isStart ? "start" : "end";

                // first start position, the end
                boolean found = false;
                for (SvCluster cluster : mClusters) {
                    SvCluster matchedCluster = cluster.findClosestCluster(currentVar.chromosome(isStart), currentVar.position(isStart));
                    if (matchedCluster != null) {
                        LOGGER.debug("variant: id({}) {}({}:{}) matched with cluster({}) on spanningVariant({})",
                                currentVar.id(), isStartStr, currentVar.chromosome(isStart), currentVar.position(isStart),
                                matchedCluster.getClusterId(), matchedCluster.getSpanningSV().posId());

                        if(isStart)
                            currentVar.setStartCluster(matchedCluster);
                        else
                            currentVar.setEndCluster(matchedCluster);

                        found = true;
                        break;
                    }
                }

                if (!found) {
                    LOGGER.debug("variant: id({}) {}({}:{}) unmatched",
                            currentVar.id(), isStartStr, currentVar.chromosome(isStart), currentVar.position(isStart));
                }
            }
        }
    }

    */
}
