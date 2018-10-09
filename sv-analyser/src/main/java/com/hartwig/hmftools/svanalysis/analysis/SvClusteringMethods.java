package com.hartwig.hmftools.svanalysis.analysis;

import static java.lang.Math.abs;
import static java.lang.Math.round;

import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.DEL;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.DUP;
import static com.hartwig.hmftools.svanalysis.analysis.SvUtilities.CHROMOSOME_ARM_Q;
import static com.hartwig.hmftools.svanalysis.analysis.SvUtilities.SV_GROUP_ENCLOSED;
import static com.hartwig.hmftools.svanalysis.analysis.SvUtilities.SV_GROUP_ENCLOSING;
import static com.hartwig.hmftools.svanalysis.analysis.SvUtilities.SV_GROUP_NEIGHBOURS;
import static com.hartwig.hmftools.svanalysis.analysis.SvUtilities.SV_GROUP_OVERLAP;
import static com.hartwig.hmftools.svanalysis.analysis.SvUtilities.isOverlapping;
import static com.hartwig.hmftools.svanalysis.annotators.LineElementAnnotator.NO_LINE_ELEMENT;
import static com.hartwig.hmftools.svanalysis.types.SvCNData.CN_SEG_TELOMERE;
import static com.hartwig.hmftools.svanalysis.types.SvClusterData.RELATION_TYPE_NEIGHBOUR;
import static com.hartwig.hmftools.svanalysis.types.SvClusterData.RELATION_TYPE_OVERLAP;

import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantType;
import com.hartwig.hmftools.svanalysis.types.SvArmGroup;
import com.hartwig.hmftools.svanalysis.types.SvBreakend;
import com.hartwig.hmftools.svanalysis.types.SvCNData;
import com.hartwig.hmftools.svanalysis.types.SvClusterData;
import com.hartwig.hmftools.svanalysis.types.SvFootprint;
import com.hartwig.hmftools.svanalysis.types.SvLOH;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class SvClusteringMethods {

    private static final Logger LOGGER = LogManager.getLogger(SvClusteringMethods.class);

    private final SvUtilities mUtils;

    private Map<String, Integer> mChrArmSvCount;
    private Map<String, Double> mChrArmSvExpected;
    private Map<String, Double> mChrArmSvRate;
    private double mMedianChrArmRate;
    private int mNextClusterId;

    private Map<String, List<SvBreakend>> mChrBreakendMap; // every breakend on a chromosome, ordered by asending position
    private Map<String, Integer> mChrCopyNumberMap; // copy number for whole chromosomes if clear
    private Map<String, List<SvCNData>> mChrCNDataMap; // copy number segments recreated from SVs
    private Map<String, List<SvLOH>> mSampleLohData;

    private static final double REF_BASE_LENGTH = 10000000D;

    public SvClusteringMethods(final SvUtilities clusteringUtils)
    {
        mUtils = clusteringUtils;
        mChrArmSvCount = Maps.newHashMap();
        mChrArmSvExpected = Maps.newHashMap();
        mChrArmSvRate = Maps.newHashMap();
        mMedianChrArmRate = 0;
        mNextClusterId = 0;

        mChrBreakendMap = new HashMap();
        mChrCopyNumberMap = new HashMap();
        mChrCNDataMap = new HashMap();
        mSampleLohData = null;
    }

    public Map<String, List<SvBreakend>> getChrBreakendMap() { return mChrBreakendMap; }
    public Map<String, Integer> getChrCopyNumberMap() { return mChrCopyNumberMap; }
    public int getNextClusterId() { return mNextClusterId++; }
    public void setSampleLohData(final Map<String, List<SvLOH>> data) { mSampleLohData = data; }

    public void clusterByBaseDistance(List<SvClusterData> allVariants, List<SvCluster> clusters)
    {
        mNextClusterId = 0;

        List<SvClusterData> unassignedVariants = Lists.newArrayList(allVariants);

        // assign each variant once to a cluster using proximity as a test
        int currentIndex = 0;

        while(currentIndex < unassignedVariants.size())
        {
            SvClusterData currentVar = unassignedVariants.get(currentIndex);

            // make a new cluster
            SvCluster newCluster = new SvCluster(getNextClusterId(), mUtils);

            // first remove the current SV from consideration
            newCluster.addVariant(currentVar);
            unassignedVariants.remove(currentIndex); // index will remain the same and so point to the next item

            // and then search for all other linked ones
            findLinkedSVsByDistance(newCluster, unassignedVariants);

            clusters.add(newCluster);
        }
    }

    private void findLinkedSVsByDistance(SvCluster cluster, List<SvClusterData> unassignedVariants)
    {
        // look for any other SVs which form part of this cluster based on proximity
        int currentIndex = 0;

        while (currentIndex < unassignedVariants.size())
        {
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
                // LOGGER.debug("cluster({}) add matched variant({}), totalCount({})", cluster.getId(), currentVar.id(), cluster.getSVs().size());

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

    public void mergeClusters(final String sampleId, List<SvCluster> clusters)
    {
        int initClusterCount = clusters.size();

        int iterations = 0;

        // the merge must be run a few times since as clusters grow, more single SVs and other clusters
        // will then fall within the bounds of the new larger clusters
        boolean foundMerges = mergeOnCommonArmLinks(clusters);

        while(foundMerges && iterations < 5)
        {
            foundMerges = mergeOnCommonArmLinks(clusters);
            ++iterations;
        }

        mergeOnLOHEvents(sampleId, clusters);

        if(clusters.size() < initClusterCount)
        {
            LOGGER.info("reduced cluster count({} -> {}) iterations({})", initClusterCount, clusters.size(), iterations);
        }
    }

    private void mergeOnLOHEvents(final String sampleId, List<SvCluster> clusters)
    {
        if(mSampleLohData == null)
            return;

        // first extract all the SVs from the LOH events
        final List<SvLOH> lohList = mSampleLohData.get(sampleId);

        if(lohList == null)
            return;

        for(final SvLOH lohEvent : lohList)
        {
            if(lohEvent.StartSV.isEmpty() || lohEvent.EndSV.isEmpty() || lohEvent.StartSV.equals("0") || lohEvent.EndSV.equals("0"))
                continue;

            SvCluster startCluster = findClusterFromVariantId(lohEvent.StartSV, clusters);
            SvCluster endCluster = findClusterFromVariantId(lohEvent.EndSV, clusters);

            if(startCluster == null)
            {
                LOGGER.error("sample({}) start varId({}) not found in any cluster", sampleId, lohEvent.StartSV);
                continue;
            }
            else if(endCluster == null)
            {
                LOGGER.error("sample({}) end varId({}) not found in any cluster", sampleId, lohEvent.EndSV);
                continue;
            }

            if(startCluster == endCluster)
                continue;

            LOGGER.debug("cluster({} svs={}) merges in other cluster({} svs={}) on LOH event(sv1={} sv2={} len={})",
                    startCluster.getId(), startCluster.getCount(), endCluster.getId(), endCluster.getCount(),
                    lohEvent.StartSV, lohEvent.EndSV, lohEvent.Length);

            // merge the second cluster into the first
            for(final SvClusterData var : endCluster.getSVs())
            {
                startCluster.addVariant(var);
            }

            clusters.remove(endCluster);
        }
    }

    private SvCluster findClusterFromVariantId(final String varId, List<SvCluster> clusters)
    {
        for(SvCluster cluster : clusters)
        {
            for(final SvClusterData var : cluster.getSVs())
            {
                if(var.id().equals(varId))
                    return cluster;
            }
        }

        return null;
    }

    private boolean mergeOnCommonArmLinks(List<SvCluster> clusters)
    {
        // merge any clusters which have the same inter-arm links and return true if merges were found
        int initClusterCount = clusters.size();

        int index1 = 0;
        while(index1 < clusters.size())
        {
            SvCluster cluster1 = clusters.get(index1);

            List<SvArmGroup> armGroups1 = cluster1.getArmGroups();

            if (armGroups1.size() == 1)
            {
                ++index1;
                continue;
            }

            if(cluster1.getCount() == 2 && cluster1.isConsistent())
            {
                // skip balanced inversion & translocations
                ++index1;
                continue;
            }

            int index2 = index1 + 1;
            while(index2 < clusters.size())
            {
                SvCluster cluster2 = clusters.get(index2);
                List<SvArmGroup> armGroups2 = cluster2.getArmGroups();

                if(armGroups2.size() == 1)
                {
                    ++index2;
                    continue;
                }

                if(cluster2.getCount() == 2 && cluster2.isConsistent())
                {
                    // skip balanced inversion & translocations
                    ++index2;
                    continue;
                }

                int sharedArmCount = 0;

                for(SvArmGroup armGroup1 : armGroups1)
                {
                    for(SvArmGroup armGroup2 : armGroups2)
                    {
                        if(armGroup1.matches(armGroup2))
                        {
                            ++sharedArmCount;

                            if(sharedArmCount >= 2)
                                break;
                        }
                    }

                    if(sharedArmCount >= 2)
                        break;
                }

                if(sharedArmCount < 2)
                {
                    ++index2;
                    continue;
                }

                LOGGER.debug("cluster({} svs={}) merges in other cluster({} svs={}) on common arms",
                        cluster1.getId(), cluster1.getCount(), cluster2.getId(), cluster2.getCount());

                // merge the second cluster into the first
                for(final SvClusterData var : cluster2.getSVs())
                {
                    cluster1.addVariant(var);
                }

                clusters.remove(index2);
            }

            ++index1;
        }

        return clusters.size() < initClusterCount;
    }

    private boolean mergeOnOverlaps(List<SvCluster> clusters)
    {
        // merge any overlapping complex clusters and return true if merges were found
        int initClusterCount = clusters.size();

        int index1 = 0;
        while(index1 < clusters.size())
        {
            SvCluster cluster1 = clusters.get(index1);

            List<SvArmGroup> armGroups1 = cluster1.getArmGroups();

            int index2 = index1 + 1;
            while(index2 < clusters.size())
            {
                SvCluster cluster2 = clusters.get(index2);
                List<SvArmGroup> armGroups2 = cluster2.getArmGroups();

                boolean canMergeArms = false;

                for(SvArmGroup armGroup1 : armGroups1)
                {
                    for(SvArmGroup armGroup2 : armGroups2)
                    {
                        if (canMergeArmGroups(armGroup1, armGroup2))
                        {
                            LOGGER.debug("arm({} svs={}) overlaps with arm({} svs={})",
                                    armGroup1.posId(), armGroup1.getCount(), armGroup2.posId(), armGroup2.getCount());

                            canMergeArms = true;
                            break;
                        }
                    }

                    if(canMergeArms)
                        break;
                }

                if(!canMergeArms)
                {
                    ++index2;
                    continue;
                }

                LOGGER.debug("cluster({} svs={}) merges in other cluster({} svs={})",
                        cluster1.getId(), cluster1.getCount(), cluster2.getId(), cluster2.getCount());

                // merge the second cluster into the first
                for(final SvClusterData var : cluster2.getSVs())
                {
                    cluster1.addVariant(var);
                }

                clusters.remove(index2);
            }

            ++index1;
        }

        return clusters.size() < initClusterCount;
    }

    private static int COMPLEX_CLUSTER_COUNT = 4;

    private boolean canMergeArmGroups(final SvArmGroup group1, final SvArmGroup group2)
    {
        if(!group1.chromosome().equals(group2.chromosome()) || !group1.arm().equals(group2.arm()))
            return false;

        // don't merge single-variant clusters which have no overlap with the other cluster
        final List<SvClusterData> list1 = group1.getSVs();
        final List<SvClusterData> list2 = group2.getSVs();
        int count1 = list1.size();
        int count2 = list2.size();

        final SvClusterData var1 = count1 == 1 ? list1.get(0) : null;
        final SvClusterData var2 = count2 == 1 ? list2.get(0) : null;

        boolean logDetails = false;

        /*
        // if(group1.chromosome().equals("15") && group1.arm().equals("Q") && (var1 != null || var2 != null))
        if(group1.chromosome().equals("15") && group1.arm().equals("Q")
        && ((var1 != null && var1.id().equals("14232")) || (var2 != null && var2.id().equals("14232"))))
        {
            logDetails = true;
        }
        */

        if(var1 != null && var1.getNearestSvRelation() == RELATION_TYPE_NEIGHBOUR)
            return false;
        else if(var2 != null && var2.getNearestSvRelation() == RELATION_TYPE_NEIGHBOUR)
            return false;

        // don't merge simple DUPs and DELs
        if(var1 != null && var2 != null)
        {
            if((var1.type() == DUP || var1.type() == DEL) && (var2.type() == DUP || var2.type() == DEL))
            {
                if(logDetails)
                    LOGGER.debug("vars({} and {}) are simple SVs", var1.posId(), var2.posId());

                return false;
            }
        }

        // check for overlapping SVs based on the outer start and end positions
        if((group1.hasEndSet() && group1.posEnd() < group2.posStart()) || (group2.hasEndSet() && group1.posStart() > group2.posEnd()))
        {
            if(logDetails)
                LOGGER.debug("groups({} and {}) outside range", group1.posId(), group2.posId());

            return false;
        }

        // there is either overlap for full enclosure for these 2 groups

        if(count1 >= COMPLEX_CLUSTER_COUNT || count2 >= COMPLEX_CLUSTER_COUNT)
            return true;

        // if no Svs in a small cluster overlap with any of the SVs in the other cluster, don't merge
        if(hasOverlappingSVs(list1, list2))
            return true;

        if(logDetails)
            LOGGER.debug("no overlapping vars for groups({} and {}) outside range", group1.posId(), group2.posId());

        return false;
    }

    private boolean hasOverlappingSVs(final List<SvClusterData> list1, final List<SvClusterData> list2)
    {
        for(final SvClusterData var1 : list1)
        {
            for(final SvClusterData var2 : list2)
            {
                if(isOverlapping(var1, var2))
                    return true;
            }
        }

        return false;
    }

    public void populateChromosomeBreakendMap(List<SvClusterData> allVariants)
    {
        mChrBreakendMap.clear();

        // add each SV's breakends to a map keyed by chromosome, with the breakends in order of position lowest to highest
        for(final SvClusterData var : allVariants)
        {
            // add each breakend in turn
            for(int i = 0; i < 2 ; ++i)
            {
                boolean useStart = (i == 0);

                if(!useStart && var.isNullBreakend())
                    continue;

                final String chr = var.chromosome(useStart);
                long position = var.position(useStart);

                if(!mChrBreakendMap.containsKey(chr))
                {
                    List<SvBreakend> breakendList = Lists.newArrayList();
                    breakendList.add(new SvBreakend(var, useStart));
                    mChrBreakendMap.put(chr, breakendList);
                    continue;
                }

                // otherwise add the variant in order by ascending position
                List<SvBreakend> breakendList = mChrBreakendMap.get(chr);

                int index = 0;
                for(;index < breakendList.size(); ++index)
                {
                    final SvBreakend breakend = breakendList.get(index);

                    if(position < breakend.position())
                        break;
                }

                breakendList.add(index, new SvBreakend(var, useStart));
            }
        }
    }

    private static double MAX_COPY_NUMBER_INTEGER_DIFF = 0.2;

    public void calcCopyNumberData(final String sampleId)
    {
        // look for duplication on each chromosome
        for (final Map.Entry<String, List<SvBreakend>> entry : mChrBreakendMap.entrySet())
        {
            final String chromosome = entry.getKey();
            final List<SvBreakend> breakendList = entry.getValue();

            int minCopyNumber = -1;
            boolean isValid = true;

            for (final SvBreakend breakend : breakendList)
            {
                final SvClusterData var = breakend.getSV();

                double copyNumber = var.copyNumber(breakend.usesStart());
                int copyNumRounded = (int)round(copyNumber);
                if(abs(copyNumber - copyNumRounded) > MAX_COPY_NUMBER_INTEGER_DIFF)
                {
                    isValid = false;
                    break;
                }

                int chromatidCopyNumber = copyNumRounded / 2;

                if(minCopyNumber < 0 || minCopyNumber > chromatidCopyNumber)
                {
                    minCopyNumber = chromatidCopyNumber;
                }
            }

            if(minCopyNumber > 1 && isValid)
            {
                LOGGER.debug("sample({}) chromosome({}) has copyNumber({})", sampleId, chromosome, minCopyNumber);
                mChrCopyNumberMap.put(chromosome, minCopyNumber);
            }
            else
            {
                mChrCopyNumberMap.put(chromosome, 1);
            }
        }
    }

    public void createCopyNumberSegments()
    {
        mChrCNDataMap.clear();

        int cnId = 0;

        for(Map.Entry<String, List<SvBreakend>> entry : mChrBreakendMap.entrySet())
        {
            final String chromosome = entry.getKey();

            List<SvBreakend> breakendList = entry.getValue();

            List<SvCNData> copyNumberList = Lists.newArrayList();

            for (int i = 0; i < breakendList.size(); ++i)
            {
                final SvBreakend breakend = breakendList.get(i);
                final SvClusterData var = breakend.getSV();

                if(i == 0)
                {
                    // add a section for the telomere - note: copy number is unknown from just the SVs
                    copyNumberList.add(
                            new SvCNData(cnId++, chromosome, 0, breakend.position() - 1, CN_SEG_TELOMERE, var.type().toString(),
                                    0, 0, 0, 2, ""));
                }

                double copyNumber = var.copyNumber(breakend.usesStart());

                long nextPosition = 0;
                String nextType;

                if(i == breakendList.size() - 1)
                {
                    // last entry - add the other telomere
                    if(!mUtils.CHROMOSOME_LENGTHS.containsKey(chromosome))
                    {
                        // LOGGER.debug("missing chromosome({})", chromosome);
                        continue;
                    }

                    nextPosition = mUtils.CHROMOSOME_LENGTHS.get(chromosome);
                    nextType = CN_SEG_TELOMERE;
                }
                else
                {
                    final SvBreakend nextBreakend = breakendList.get(i+1);
                    final SvClusterData nextVar = nextBreakend.getSV();
                    nextType = nextVar.type().toString();
                    nextPosition = nextBreakend.position() - 1;
                }

                copyNumberList.add(
                        new SvCNData(cnId++, chromosome, breakend.position(), nextPosition, var.type().toString(), nextType, 0, 0, 0, copyNumber, ""));
            }

            mChrCNDataMap.put(chromosome, copyNumberList);
        }
    }

    public void annotateNearestSvData()
    {
        // mark each SV's nearest other SV and its relationship - neighbouring or overlapping
        for(Map.Entry<String, List<SvBreakend>> entry : mChrBreakendMap.entrySet())
        {
            List<SvBreakend> breakendList = entry.getValue();

            int breakendCount = breakendList.size();

            for(int i = 0; i < breakendList.size(); ++i)
            {
                final SvBreakend breakend = breakendList.get(i);
                SvClusterData var = breakend.getSV();

                /*
                if(var.id().equals("52520"))
                {
                    LOGGER.debug("s");
                }
                */

                SvBreakend prevBreakend = (i > 0) ? breakendList.get(i - 1) : null;
                SvBreakend nextBreakend = (i < breakendCount-1) ? breakendList.get(i + 1) : null;

                long closestDistance = -1;
                if(prevBreakend != null && prevBreakend.getSV() != var)
                {
                    long distance = breakend.position() - prevBreakend.position();
                    closestDistance = distance;
                }

                if(nextBreakend != null && nextBreakend.getSV() != var)
                {
                    long distance = nextBreakend.position() - breakend.position();
                    if(closestDistance < 0 || distance < closestDistance)
                        closestDistance = distance;
                }

                if(closestDistance >= 0 && (var.getNearestSvDistance() == -1 || closestDistance < var.getNearestSvDistance()))
                    var.setNearestSvDistance(closestDistance);

                String relationType = "";
                if((prevBreakend != null && prevBreakend.getSV() == var) || (nextBreakend != null && nextBreakend.getSV() == var))
                    relationType = RELATION_TYPE_NEIGHBOUR;
                else
                    relationType = RELATION_TYPE_OVERLAP;

                var.setNearestSvRelation(relationType);
            }
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

            // ensure an entry exists
            if (!mChrArmSvCount.containsKey(chrArmStart))
            {
                mChrArmSvCount.put(chrArmStart, 0);
            }

            // exclude LINE elements from back-ground rates
            if(var.isStartLineElement().equals(NO_LINE_ELEMENT))
            {
                mChrArmSvCount.replace(chrArmStart, mChrArmSvCount.get(chrArmStart) + 1);
            }

            if(!var.isNullBreakend())
            {
                if (!chrArmStart.equals(chrArmEnd) && !mChrArmSvCount.containsKey(chrArmEnd))
                {
                    mChrArmSvCount.put(chrArmEnd, 0);
                }

                if (var.isEndLineElement().equals(NO_LINE_ELEMENT))
                {
                    mChrArmSvCount.replace(chrArmEnd, mChrArmSvCount.get(chrArmEnd) + 1);
                }
            }
        }

        // now determine the background rate by taking the median value from amongst the arms
        // factoring in the arms which have no Q (14-16, 21-22) and excluding the X & Ys
        for(Map.Entry<String, Integer> entry : mChrArmSvCount.entrySet())
        {
            final String chrArm = entry.getKey();
            final String chromosome = mUtils.getChrFromChrArm(chrArm);
            final String arm = mUtils.getArmFromChrArm(chrArm);

            long chrArmLength = mUtils.getChromosomalArmLength(chromosome, arm);
            int svCount = entry.getValue();
            double ratePerLength = svCount / (chrArmLength / REF_BASE_LENGTH); // the factor isn't important

            mChrArmSvRate.put(chrArm, ratePerLength);
            // LOGGER.debug("chrArm({}) ratePerMill({}) from count({}) length({})", chrArm, ratePerLength, svCount, chrArmLength);
        }

        mChrArmSvRate = sortByValue(mChrArmSvRate, false);

        mMedianChrArmRate = 0;
        int chrArmIndex = 0;
        for(Map.Entry<String, Double> entry : mChrArmSvRate.entrySet())
        {
            // LOGGER.debug("chrArm({}: {}) svRate({})", chrArmIndex, entry.getKey(), entry.getValue());

            if(chrArmIndex == 20)
                mMedianChrArmRate = entry.getValue();

           ++chrArmIndex;
        }

        LOGGER.debug(String.format("median SV rate(%.2f)", mMedianChrArmRate));

        // now create another map of expected SV count per arm using the median rate
        for(Map.Entry<String, Double> entry : mChrArmSvRate.entrySet())
        {
            final String chrArm = entry.getKey();
            final String chromosome = mUtils.getChrFromChrArm(chrArm);
            final String arm = mUtils.getArmFromChrArm(chrArm);

            long chrArmLength = mUtils.getChromosomalArmLength(chromosome, arm);
            double expectedSvCount = (int) round((chrArmLength / REF_BASE_LENGTH) * mMedianChrArmRate);
            LOGGER.debug("chrArm({}) expectedSvCount({}) vs actual({})", chrArm, expectedSvCount, mChrArmSvCount.get(chrArm));

            mChrArmSvExpected.put(chrArm, expectedSvCount);
        }
    }

    public String getChrArmData(final SvClusterData var)
    {
        String chrArmStart = mUtils.getVariantChrArm(var,true);

        boolean hasEnd = !var.isNullBreakend();
        String chrArmEnd = hasEnd ? mUtils.getVariantChrArm(var,false) : "";

        // report Start SV count : Expected SV Count : End SV Count : Expected SV Count
        return String.format("%d,%.2f,%d,%.2f",
                mChrArmSvCount.get(chrArmStart), mChrArmSvExpected.get(chrArmStart),
                hasEnd ? mChrArmSvCount.get(chrArmEnd) : 0, hasEnd ? mChrArmSvExpected.get(chrArmEnd) : 0.0);
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

}
