package com.hartwig.hmftools.svanalysis.analysis;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.round;

import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.BND;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.DEL;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.DUP;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.INS;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.INV;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.SGL;
import static com.hartwig.hmftools.svanalysis.analysis.ClusterAnalyser.SMALL_CLUSTER_SIZE;
import static com.hartwig.hmftools.svanalysis.analysis.LinkFinder.MIN_TEMPLATED_INSERTION_LENGTH;
import static com.hartwig.hmftools.svanalysis.analysis.LinkFinder.areLinkedSection;
import static com.hartwig.hmftools.svanalysis.analysis.LinkFinder.areSectionBreak;
import static com.hartwig.hmftools.svanalysis.analysis.SvUtilities.CHROMOSOME_ARM_P;
import static com.hartwig.hmftools.svanalysis.analysis.SvUtilities.CHROMOSOME_ARM_Q;
import static com.hartwig.hmftools.svanalysis.analysis.SvUtilities.calcConsistency;
import static com.hartwig.hmftools.svanalysis.analysis.SvUtilities.isLocalOverlap;
import static com.hartwig.hmftools.svanalysis.analysis.SvUtilities.makeChrArmStr;
import static com.hartwig.hmftools.svanalysis.types.SvCluster.RESOLVED_LOW_QUALITY;
import static com.hartwig.hmftools.svanalysis.types.SvCluster.RESOLVED_TYPE_DEL_EXT_TI;
import static com.hartwig.hmftools.svanalysis.types.SvCluster.RESOLVED_TYPE_DEL_INT_TI;
import static com.hartwig.hmftools.svanalysis.types.SvCluster.RESOLVED_TYPE_DUP_EXT_TI;
import static com.hartwig.hmftools.svanalysis.types.SvCluster.RESOLVED_TYPE_DUP_INT_TI;
import static com.hartwig.hmftools.svanalysis.types.SvCluster.RESOLVED_TYPE_RECIPROCAL_TRANS;
import static com.hartwig.hmftools.svanalysis.types.SvCluster.RESOLVED_TYPE_SIMPLE_INS;
import static com.hartwig.hmftools.svanalysis.types.SvCluster.findCluster;
import static com.hartwig.hmftools.svanalysis.analysis.SvUtilities.calcTypeCount;
import static com.hartwig.hmftools.svanalysis.types.SvCNData.CN_SEG_TELOMERE;
import static com.hartwig.hmftools.svanalysis.types.SvLinkedPair.LINK_TYPE_DB;
import static com.hartwig.hmftools.svanalysis.types.SvLinkedPair.LINK_TYPE_TI;
import static com.hartwig.hmftools.svanalysis.types.SvVarData.ASSEMBLY_TYPE_EQV;
import static com.hartwig.hmftools.svanalysis.types.SvVarData.RELATION_TYPE_NEIGHBOUR;
import static com.hartwig.hmftools.svanalysis.types.SvVarData.RELATION_TYPE_OVERLAP;
import static com.hartwig.hmftools.svanalysis.types.SvVarData.haveLinkedAssemblies;

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
import com.hartwig.hmftools.svanalysis.types.SvChain;
import com.hartwig.hmftools.svanalysis.types.SvCluster;
import com.hartwig.hmftools.svanalysis.types.SvLinkedPair;
import com.hartwig.hmftools.svanalysis.types.SvVarData;
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

    private long mDelDupCutoffLength;

    private boolean mFoldbackColsLogged;
    private boolean mInversionPairColsLogged;

    private static double REF_BASE_LENGTH = 10000000D;
    public static int MAX_SIMPLE_DUP_DEL_CUTOFF = 5000000;
    public static int MIN_SIMPLE_DUP_DEL_CUTOFF = 100000;

    public static String CLUSTER_REASON_PROXIMITY = "Proximity";
    public static String CLUSTER_REASON_LOH = "LOH";
    public static String CLUSTER_REASON_COMMON_ARMS = "CommonArms";
    public static String CLUSTER_REASON_FOLDBACKS = "Foldback";
    public static String CLUSTER_REASON_SOLO_SINGLE = "SoloSingle";
    public static String CLUSTER_REASON_INV_OVERLAP = "InvOverlap";
    public static String CLUSTER_REASON_LONG_DEL_DUP = "LongDelDup";

    public SvClusteringMethods(final SvUtilities clusteringUtils)
    {
        mUtils = clusteringUtils;
        mChrArmSvCount = Maps.newHashMap();
        mChrArmSvExpected = Maps.newHashMap();
        mChrArmSvRate = Maps.newHashMap();
        mMedianChrArmRate = 0;
        mNextClusterId = 0;

        mDelDupCutoffLength = 0;

        mChrBreakendMap = new HashMap();
        mChrCopyNumberMap = new HashMap();
        mChrCNDataMap = new HashMap();
        mSampleLohData = null;
        mFoldbackColsLogged = false;
        mInversionPairColsLogged = false;
    }

    public Map<String, List<SvBreakend>> getChrBreakendMap() { return mChrBreakendMap; }
    public int getNextClusterId() { return mNextClusterId++; }
    public void setSampleLohData(final Map<String, List<SvLOH>> data) { mSampleLohData = data; }

    public void clusterByBaseDistance(List<SvVarData> allVariants, List<SvCluster> clusters)
    {
        mNextClusterId = 0;

        List<SvVarData> unassignedVariants = Lists.newArrayList(allVariants);

        // assign each variant once to a cluster using proximity as a test
        int currentIndex = 0;

        while(currentIndex < unassignedVariants.size())
        {
            SvVarData currentVar = unassignedVariants.get(currentIndex);

            // make a new cluster
            SvCluster newCluster = new SvCluster(getNextClusterId());

            // first remove the current SV from consideration
            newCluster.addVariant(currentVar);
            unassignedVariants.remove(currentIndex); // index will remain the same and so point to the next item

            // exceptions to proximity clustering
            if(isEquivSingleBreakend(currentVar))
            {
                newCluster.setResolved(true, RESOLVED_LOW_QUALITY);
                clusters.add(newCluster);
                continue;
            }

            // and then search for all other linked ones
            findLinkedSVsByDistance(newCluster, unassignedVariants, true);

            // check for invalid clusters - currently just overlapping DELs
            if(newCluster.getCount() > 1 && newCluster.getCount() == newCluster.getTypeCount(DEL))
            {
                addDelGroupClusters(clusters, newCluster);
            }
            else
            {
                clusters.add(newCluster);
            }
        }
    }

    private void addDelGroupClusters(List<SvCluster> clusters, SvCluster delGroupCluster)
    {
        // only cluster if proximate and not overlapping
        int currentIndex = 0;

        List<SvVarData> delSVs = Lists.newArrayList();
        delSVs.addAll(delGroupCluster.getSVs());

        while(currentIndex < delSVs.size())
        {
            SvVarData currentVar = delSVs.get(currentIndex);

            // make a new cluster
            SvCluster newCluster = new SvCluster(getNextClusterId());

            // first remove the current SV from consideration
            newCluster.addVariant(currentVar);
            delSVs.remove(currentIndex); // index will remain the same and so point to the next item

            // and then search for all other linked ones
            findLinkedSVsByDistance(newCluster, delSVs, false);

            clusters.add(newCluster);
        }
    }

    private void findLinkedSVsByDistance(SvCluster cluster, List<SvVarData> unassignedVariants, boolean allowOverlaps)
    {
        // look for any other SVs which form part of this cluster based on proximity
        int currentIndex = 0;

        while (currentIndex < unassignedVariants.size())
        {
            SvVarData currentVar = unassignedVariants.get(currentIndex);

            if(isEquivSingleBreakend(currentVar))
            {
                ++currentIndex;
                continue;
            }

            // compare with all other SVs in this cluster
            boolean matched = false;
            for (SvVarData otherVar : cluster.getSVs())
            {
                    // test each possible linkage
                if (!mUtils.areVariantsLinkedByDistance(currentVar, otherVar))
                    continue;

                if(!allowOverlaps && !(currentVar.position(false) < otherVar.position(true) || otherVar.position(false) < currentVar.position(true)))
                    continue;

                cluster.addVariant(currentVar);
                currentVar.addClusterReason(CLUSTER_REASON_PROXIMITY, otherVar.id());

                if(otherVar.getClusterReason().isEmpty())
                    otherVar.addClusterReason(CLUSTER_REASON_PROXIMITY, currentVar.id());

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
        // first apply replication rules since this can affect consistency
        for(SvCluster cluster : clusters)
        {
            // check for sub-clonal / low-copy-number supported variants
            if(cluster.getCount() == 1 && isLowQualityVariant(cluster.getSVs().get(0)))
            {
                cluster.setResolved(true, RESOLVED_LOW_QUALITY);
                continue;
            }

            applyCopyNumberReplication(sampleId, cluster);
        }

        int initClusterCount = clusters.size();

        int iterations = 0;

        // the merge must be run a few times since as clusters grow, more single SVs and other clusters
        // will then fall within the bounds of the new larger clusters
        boolean foundMerges = true;

        while(foundMerges && iterations < 5)
        {
            foundMerges = mergeOnCommonArmLinks(clusters);

            if(mergeOnUnresolvedSVs(clusters, true, false))
                foundMerges = true;

            if(mergeLongDupsDels(clusters))
                foundMerges = true;

            ++iterations;
        }

        mergeOnLOHEvents(sampleId, clusters);

        if(clusters.size() < initClusterCount)
        {
            LOGGER.debug("reduced cluster count({} -> {}) iterations({})", initClusterCount, clusters.size(), iterations);
        }
    }

    public static void applyCopyNumberReplication(final String sampleId, SvCluster cluster)
    {
        // use the relative copy number change to replicate some SVs within a cluster
        if(!cluster.hasVariedCopyNumber())
            return;

        // first establish the lowest copy number change
        int minCopyNumber = cluster.getMinCopyNumber();

        // replicate the SVs which have a higher copy number than their peers
        int clusterCount = cluster.getCount();

        for(int i = 0; i < clusterCount; ++i)
        {
            SvVarData var = cluster.getSVs().get(i);
            int calcCopyNumber = var.impliedCopyNumber(true);

            if(calcCopyNumber <= minCopyNumber)
                continue;

            int svMultiple = calcCopyNumber / minCopyNumber;

            LOGGER.debug("sample({}) replicating SV({}) {} times, copyNumChg({} vs min={})",
                    sampleId, var.posId(), svMultiple, calcCopyNumber, minCopyNumber);

            var.setReplicatedCount(svMultiple);

            for(int j = 1; j < svMultiple; ++j)
            {
                SvVarData newVar = new SvVarData(var);
                cluster.addVariant(newVar);
            }
        }
    }

    public static double LOW_QUALITY_CN_CHANGE = 0.5;

    public boolean isLowQualityVariant(final SvVarData var)
    {
        return var.copyNumberChange(true) < LOW_QUALITY_CN_CHANGE && var.copyNumberChange(false) < LOW_QUALITY_CN_CHANGE;
    }

    public boolean isEquivSingleBreakend(final SvVarData var)
    {
        if(!var.isNullBreakend())
            return false;

        if(var.isDupBreakend(true) || var.isDupBreakend(false))
            return true;

        return  (var.getAssemblyData(true).contains(ASSEMBLY_TYPE_EQV)
                || var.getAssemblyData(false).contains(ASSEMBLY_TYPE_EQV));
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

            SvCluster lohClusterStart = null;
            SvVarData lohSvStart = null;
            SvCluster lohClusterEnd = null;
            SvVarData lohSvEnd = null;

            for(SvCluster cluster : clusters)
            {
                for(final SvVarData var : cluster.getSVs())
                {
                    if(var.id().equals(lohEvent.StartSV))
                    {
                        lohClusterStart = cluster;
                        lohSvStart = var;
                    }
                    if(var.id().equals(lohEvent.EndSV))
                    {
                        lohClusterEnd = cluster;
                        lohSvEnd = var;
                    }

                    if(lohSvStart != null && lohSvEnd != null)
                        break;
                }
            }

            if(lohClusterEnd == null || lohClusterStart == null)
            {
                // LOGGER.error("sample({}) start varId({}) not found in any cluster", sampleId, lohEvent.StartSV);
                continue;
            }

            if(lohClusterStart == lohClusterEnd)
                continue;

            LOGGER.debug("cluster({} svs={}) merges in other cluster({} svs={}) on LOH event(sv1={} sv2={} len={})",
                    lohClusterStart.getId(), lohClusterStart.getUniqueSvCount(), lohClusterEnd.getId(), lohClusterEnd.getCount(),
                    lohEvent.StartSV, lohEvent.EndSV, lohEvent.Length);

            lohSvStart.addClusterReason(CLUSTER_REASON_LOH, lohSvEnd.id());
            lohSvEnd.addClusterReason(CLUSTER_REASON_LOH, lohSvStart.id());

            lohClusterStart.mergeOtherCluster(lohClusterEnd);
            clusters.remove(lohClusterEnd);
        }
    }

    public static boolean isSmallConsistentCluster(final SvCluster cluster)
    {
        if(cluster.getUniqueSvCount() > SMALL_CLUSTER_SIZE)
            return false;

        return isConsistentCluster(cluster);
    }

    public static boolean isConsistentCluster(final SvCluster cluster)
    {
        if(cluster.isSimpleSVs() || cluster.isResolved())
            return true;

        if(!cluster.isConsistent())
            return false;

        if(cluster.isFullyChained())
            return true;

        // other wise check whether all remaining SVs are either in consistent chains
        // or themselves consistent
        for(final SvChain chain : cluster.getChains())
        {
            if(!chain.isConsistent())
                return false;
        }

        for(final SvVarData var : cluster.getUnlinkedSVs())
        {
            if(!var.isLocal()) // so filters out cross arm, null breakends and translocation
                return false;

            if(calcConsistency(var) != 0)
                return false;
        }

        return true;
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

            if(isSmallConsistentCluster(cluster1))
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

                if(isSmallConsistentCluster(cluster2))
                {
                    // skip balanced inversion & translocations
                    ++index2;
                    continue;
                }

                // look at all SVs which span 2 arms (either same chromosome or not) and match in both clusters, ignoring those in line elements
                boolean commonLinkingSvFound = false;
                String commonArms = "";

                for(final SvVarData var1 : cluster1.getSVs())
                {
                    if (var1.isLocal() || var1.isReplicatedSv() || var1.isNullBreakend())
                        continue;

                    if (var1.isLineElement(true) || var1.isLineElement(false))
                        continue;

                    for (final SvVarData var2 : cluster2.getSVs())
                    {
                        if (var2.isLocal() || var2.isReplicatedSv() || var2.isNullBreakend())
                            continue;

                        if (var2.isLineElement(true) || var2.isLineElement(false))
                            continue;

                        if (var1.chromosome(true).equals(var2.chromosome(true)) && var1.chromosome(false).equals(var2.chromosome(false))
                        && var1.arm(true).equals(var2.arm(true)) && var1.arm(false).equals(var2.arm(false)))
                        {
                            commonLinkingSvFound = true;

                            var1.addClusterReason(CLUSTER_REASON_COMMON_ARMS, var2.id());
                            var2.addClusterReason(CLUSTER_REASON_COMMON_ARMS, var1.id());
                            break;
                        }
                        else if (var1.chromosome(true).equals(var2.chromosome(false)) && var1.chromosome(false).equals(var2.chromosome(true))
                                && var1.arm(true).equals(var2.arm(false)) && var1.arm(false).equals(var2.arm(true)))
                        {
                            commonLinkingSvFound = true;

                            var1.addClusterReason(CLUSTER_REASON_COMMON_ARMS, var2.id());
                            var2.addClusterReason(CLUSTER_REASON_COMMON_ARMS, var1.id());
                            break;
                        }
                    }

                    if(commonLinkingSvFound)
                    {
                        commonArms = makeChrArmStr(var1.chromosome(true), var1.arm(true)) + " - " + makeChrArmStr(var1.chromosome(false), var1.arm(false));
                        break;
                    }
                }

                if(commonLinkingSvFound)
                {
                    LOGGER.debug("cluster({} svs={}) merges in other cluster({} svs={}) on common arms({})",
                            cluster1.getId(), cluster1.getUniqueSvCount(), cluster2.getId(), cluster2.getUniqueSvCount(), commonArms);

                    cluster1.mergeOtherCluster(cluster2);
                    clusters.remove(index2);
                    continue;
                }

                ++index2;
                continue;
            }

            ++index1;
        }

        return clusters.size() < initClusterCount;
    }

    private boolean mergeOnUnresolvedSVs(List<SvCluster> clusters, boolean checkOverlaps, boolean checkSingleLinks)
    {
        // merge any overlapping complex clusters and return true if merges were found
        // ignore simple SVs and simple, consistent cluster-2s
        int initClusterCount = clusters.size();

        int index1 = 0;
        while(index1 < clusters.size())
        {
            SvCluster cluster1 = clusters.get(index1);

            if(isSmallConsistentCluster(cluster1))
            {
                ++index1;
                continue;
            }

            List<SvArmGroup> armGroups1 = cluster1.getArmGroups();

            int index2 = index1 + 1;
            while(index2 < clusters.size())
            {
                SvCluster cluster2 = clusters.get(index2);

                if(isSmallConsistentCluster(cluster2))
                {
                    ++index2;
                    continue;
                }

                List<SvArmGroup> armGroups2 = cluster2.getArmGroups();

                boolean canMergeArms = false;

                for(SvArmGroup armGroup1 : armGroups1)
                {
                    for(SvArmGroup armGroup2 : armGroups2)
                    {
                        if (checkOverlaps)
                        {
                            final SvVarData linkingVar = armGroupsHaveOverlappingSVs(armGroup1, armGroup2, true);

                            if(linkingVar != null)
                            {
                                LOGGER.debug("arm({} svs={}) overlaps with arm({} svs={}) on SV({})",
                                        armGroup1.posId(), armGroup1.getCount(), armGroup2.posId(), armGroup2.getCount(), linkingVar);

                                linkingVar.addClusterReason(CLUSTER_REASON_INV_OVERLAP, "");

                                canMergeArms = true;
                                break;
                            }
                        }

                        if (!canMergeArms && checkSingleLinks && armGroupsHaveLinkingSVs(armGroup1, armGroup2))
                        {
                            LOGGER.debug("cluster({}) arm({} svs={}) has linking SVs with cluster({}) arm({} svs={})",
                                    cluster1.getId(), armGroup1.posId(), armGroup1.getCount(), cluster2.getId(), armGroup2.posId(), armGroup2.getCount());

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

                // merge the second cluster into the first
                cluster1.mergeOtherCluster(cluster2);
                clusters.remove(index2);
            }

            ++index1;
        }

        return clusters.size() < initClusterCount;
    }

    private boolean armGroupsHaveLinkingSVs(final SvArmGroup group1, final SvArmGroup group2)
    {
        if (!group1.chromosome().equals(group2.chromosome()) || !group1.arm().equals(group2.arm()))
            return false;

        // merge if there are INVs in both groups or BNDs in both groups
        if(calcTypeCount(group1.getSVs(), INV) > 0 && calcTypeCount(group2.getSVs(), INV) > 0)
            return true;

        if(calcTypeCount(group1.getSVs(), BND) > 0 && calcTypeCount(group2.getSVs(), BND) > 0)
            return true;

        return false;
    }

    private final SvVarData armGroupsHaveOverlappingSVs(final SvArmGroup group1, final SvArmGroup group2, boolean reqInversion)
    {
        if(!group1.chromosome().equals(group2.chromosome()) || !group1.arm().equals(group2.arm()))
            return null;

        // check for an overlapping INV from one cluster to another
        boolean hasOverlap = false;

        for(int i = 0; i < 2; ++i)
        {
            final SvArmGroup checkGroup = (i == 0) ? group1 : group2;

            if(!checkGroup.hasEndsSet())
                continue;

            final SvArmGroup svGroup = (i == 0) ? group2 : group1;

            if(reqInversion)
            {
                boolean hasInversion = false;

                for (final SvVarData var : svGroup.getSVs())
                {
                    if (var.type() == INV)
                    {
                        hasInversion = true;
                        break;
                    }
                }

                if (!hasInversion)
                    continue;
            }

            for (final SvVarData var : svGroup.getSVs())
            {
                if (var.type() != INV)
                    continue;

                if ((var.position(true) >= checkGroup.posStart() && var.position(true) <= checkGroup.posEnd())
                || (var.position(false) >= checkGroup.posStart() && var.position(false) <= checkGroup.posEnd()))
                {
                    return var;
                }
            }
        }

        return null;
    }

    private boolean mergeLongDupsDels(List<SvCluster> clusters)
    {
        if(mDelDupCutoffLength == 0)
            return false;

        int initClusterCount = clusters.size();

        int index1 = 0;
        while(index1 < clusters.size())
        {
            SvCluster cluster1 = clusters.get(index1);

            if(cluster1.isResolved())
            {
                ++index1;
                continue;
            }

            boolean hasLongDelDup1 = clusterHasLongDelDup(cluster1);

            int index2 = index1 + 1;
            while(index2 < clusters.size())
            {
                SvCluster cluster2 = clusters.get(index2);

                if(cluster2.isResolved())
                {
                    ++index2;
                    continue;
                }

                boolean hasLongDelDup2 = clusterHasLongDelDup(cluster2);

                if((hasLongDelDup1 && canMergeClustersOnLongDelDups(cluster1, cluster2))
                || (hasLongDelDup2 && canMergeClustersOnLongDelDups(cluster2, cluster1)))
                {
                    cluster1.mergeOtherCluster(cluster2);
                    clusters.remove(index2);

                    if(hasLongDelDup1)
                        addClusterReason(cluster2, cluster1, CLUSTER_REASON_LONG_DEL_DUP);

                    if(hasLongDelDup2)
                        addClusterReason(cluster1, cluster2, CLUSTER_REASON_LONG_DEL_DUP);
                }
                else
                {
                    ++index2;
                    continue;
                }
            }

            ++index1;
        }

        return clusters.size() < initClusterCount;
    }

    public static void addClusterReason(SvCluster mergedCluster, SvCluster mergingCluster, final String reason)
    {
        final String otherId = mergingCluster.getCount() == 1 ? mergingCluster.getSVs().get(0).id() : "";

        for(SvVarData var : mergedCluster.getSVs())
        {
            var.addClusterReason(reason, otherId);
        }

    }

    public boolean clusterHasLongDelDup(final SvCluster cluster)
    {
        if(cluster.isSimpleSVs())
        {
            for(final SvVarData var : cluster.getSVs())
            {
                if((var.type() == DUP || var.type() == DEL) && var.length() >= mDelDupCutoffLength)
                {
                    return true;
                }
            }
        }
        else if(cluster.isResolved())
        {
            if(cluster.getResolvedType() == RESOLVED_TYPE_DEL_EXT_TI
            || cluster.getResolvedType() == RESOLVED_TYPE_DUP_INT_TI
            ||cluster.getResolvedType() == RESOLVED_TYPE_DUP_EXT_TI)
            {
                if (cluster.getLengthOverride() >= mDelDupCutoffLength)
                    return true;
            }
        }

        return false;
    }

    public boolean canMergeClustersOnLongDelDups(final SvCluster cluster1, final SvCluster cluster2)
    {
        if(cluster2.isSimpleSVs())
        {
            // merge simple DEL-DUP combinations if both exceed the required length
            if(cluster2.getCount() == 1)
            {
                final SvVarData var1 = cluster1.getSVs().get(0);
                final SvVarData var2 = cluster2.getSVs().get(0);

                if((var2.type() == DEL || var2.type() == DUP) && var2.type() != var1.type()
                        && var2.length() >= mDelDupCutoffLength && isLocalOverlap(var1, var2))
                {
                    LOGGER.debug("long DEL-DUPs({} and {}) overlap with lengths({} and {})",
                            var1.posId(), var2.posId(), var1.length(), var2.length());
                    return true;
                }
            }
        }
        else if(cluster2.getResolvedType() == RESOLVED_TYPE_RECIPROCAL_TRANS
            || cluster2.getResolvedType() == RESOLVED_TYPE_SIMPLE_INS)
        {
            return false;
        }
        else
        {
            List<SvArmGroup> armGroups1 = cluster1.getArmGroups();
            List<SvArmGroup> armGroups2 = cluster2.getArmGroups();

            for (SvArmGroup armGroup1 : armGroups1)
            {
                for (SvArmGroup armGroup2 : armGroups2)
                {
                    final SvVarData linkingVar = armGroupsHaveOverlappingSVs(armGroup1, armGroup2, false);

                    if(linkingVar != null)
                    {
                        LOGGER.debug("long DUP-DEL arm({} svs={}) overlaps with arm({} svs={}) on SV({})",
                                armGroup1.posId(), armGroup1.getCount(), armGroup2.posId(), armGroup2.getCount(), linkingVar.id());

                        linkingVar.addClusterReason(CLUSTER_REASON_LONG_DEL_DUP, "");
                        return true;
                    }
                }
            }
        }

        return false;
    }

    private static int DEL_DUP_LENGTH_TRIM_COUNT = 5;
    private static int MAX_ARM_COUNT = 41; // excluding the 5 short arms

    public void setSimpleVariantLengths(final String sampleId)
    {
        mDelDupCutoffLength = 0;

        List<Long> lengthsList = Lists.newArrayList();

        int simpleArmCount = 0;

        for (final Map.Entry<String, List<SvBreakend>> entry : mChrBreakendMap.entrySet())
        {
            final String chromosome = entry.getKey();
            final List<SvBreakend> breakendList = entry.getValue();

            // first check for complex events on the arm since these will be skipped
            boolean pArmHasInversions = false;
            boolean qArmHasInversions = false;

            for(final SvBreakend breakend : breakendList)
            {
                final SvVarData var = breakend.getSV();

                if(var.type() != INV)
                    continue;

                if(!pArmHasInversions && (var.arm(true) == CHROMOSOME_ARM_P || var.arm(false) == CHROMOSOME_ARM_P))
                    pArmHasInversions = true;

                if(!qArmHasInversions && (var.arm(true) == CHROMOSOME_ARM_Q || var.arm(false) == CHROMOSOME_ARM_Q))
                    qArmHasInversions = true;

                if(pArmHasInversions && qArmHasInversions)
                    break;
            }

            // skip chromosome altogether
            if(pArmHasInversions && qArmHasInversions)
                continue;

            if(!pArmHasInversions)
                ++simpleArmCount;

            if(!qArmHasInversions)
                ++simpleArmCount;

            int armCount = 0;

            for(final SvBreakend breakend : breakendList)
            {
                if(!breakend.usesStart() || !(breakend.getSV().type() == DEL || breakend.getSV().type() == DUP))
                    continue;

                final SvVarData var = breakend.getSV();

                if(pArmHasInversions)
                {
                    if (var.arm(true) == CHROMOSOME_ARM_P || var.arm(false) == CHROMOSOME_ARM_P)
                        continue;
                }

                if(qArmHasInversions)
                {
                    if (var.arm(true) == CHROMOSOME_ARM_Q || var.arm(false) == CHROMOSOME_ARM_Q)
                        continue;
                }

                lengthsList.add(var.length());
                ++armCount;
            }

            // LOGGER.debug("sample({}) chr({}) svCount({} delDups({})", sampleId, chromosome, breakendList.size(), armCount);
        }

        int trimCount = (int)round(simpleArmCount / (double)MAX_ARM_COUNT * DEL_DUP_LENGTH_TRIM_COUNT);

        if(lengthsList.size() > trimCount)
        {
            Collections.sort(lengthsList);
            int lengthIndex = lengthsList.size() - trimCount - 1; // 10 items, index 0 - 9, exclude 5 - 9, select 9
            mDelDupCutoffLength = lengthsList.get(lengthIndex);
        }

        LOGGER.debug("sample({}) simple dels and dups: count({}) cutoff-length({}) simpleArms({}) trimCount({})",
                sampleId, lengthsList.size(), mDelDupCutoffLength, simpleArmCount, trimCount);

        // LOGGER.info("DEL_DUP_CUTOFF: {},{},{},{},{}",
        //        sampleId, lengthsList.size(), mDelDupCutoffLength, simpleArmCount, trimCount);

        mDelDupCutoffLength = min(max(mDelDupCutoffLength, MIN_SIMPLE_DUP_DEL_CUTOFF), MAX_SIMPLE_DUP_DEL_CUTOFF);
    }

    public void populateChromosomeBreakendMap(final List<SvVarData> allVariants)
    {
        mChrBreakendMap.clear();

        // add each SV's breakends to a map keyed by chromosome, with the breakends in order of position lowest to highest
        for(final SvVarData var : allVariants)
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

    private static double MAX_COPY_NUMBER_INTEGER_DIFF = 0.25;

    public void calcCopyNumberData(final String sampleId)
    {
        // look for duplication on each chromosome by examining the copy number heading towards the telomeres
        for (final Map.Entry<String, List<SvBreakend>> entry : mChrBreakendMap.entrySet())
        {
            final String chromosome = entry.getKey();
            final List<SvBreakend> breakendList = entry.getValue();

            int chrCopyNumber = 0;
            boolean isValid = true;

            for(int i = 0; i < 2; ++i)
            {
                final SvBreakend breakend = (i == 0) ? breakendList.get(i) : breakendList.get(breakendList.size() - 1);
                final SvVarData var = breakend.getSV();
                boolean useStart = breakend.usesStart();
                double copyNumber = var.copyNumber(useStart);

                double telomereCopyNumber;
                boolean useCopyNumber;

                if(!var.isTranslocation())
                {
                    if((i == 0 && var.orientation(true) == 1) || (i == 1 && var.orientation(false) == -1))
                        useCopyNumber = true;
                    else
                        useCopyNumber = false;
                }
                else
                {
                    if((i == 0 && var.orientation(useStart) == 1) || (i == 1 && var.orientation(useStart) == -1))
                        useCopyNumber = true;
                    else
                        useCopyNumber = false;
                }

                if(useCopyNumber)
                {
                    // coming from telemore already so just take copy number
                    telomereCopyNumber = copyNumber;
                }
                else
                {
                    // section out to telomere on this chromatid has been lost, so whatever is left needs to be doubled
                    double copyNumberChange = var.copyNumberChange(useStart);
                    telomereCopyNumber = copyNumber - copyNumberChange;
                    telomereCopyNumber *= 2;
                }

                int copyNumRounded = (int)round(telomereCopyNumber);

                /*
                if(abs(telomereCopyNumber - copyNumRounded) > MAX_COPY_NUMBER_INTEGER_DIFF)
                {
                    isValid = false;
                    continue;
                }
                */

                /*
                LOGGER.debug(String.format("chromosome(%s) %s telomere has copyNumber(%d) from SV(%s pos=%d:%d cn=%.1f cnc=%.1f)",
                        chromosome, i == 0 ? "start" : "end", copyNumRounded,
                        var.id(), var.position(useStart), var.orientation(useStart),
                        var.copyNumber(useStart), var.copyNumberChange(useStart)));
                */

                if(i == 0)
                    chrCopyNumber = copyNumRounded;
                else
                    chrCopyNumber = min(chrCopyNumber, copyNumRounded);
            }

            if(chrCopyNumber > 2 && isValid)
            {
                LOGGER.debug("sample({}) chromosome({}) has copyNumber({})", sampleId, chromosome, chrCopyNumber);
                mChrCopyNumberMap.put(chromosome, chrCopyNumber);
            }
            else
            {
                mChrCopyNumberMap.put(chromosome, 2);
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
                final SvVarData var = breakend.getSV();

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
                    final SvVarData nextVar = nextBreakend.getSV();
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
                SvVarData var = breakend.getSV();

                final SvBreakend prevBreakend = (i > 0) ? breakendList.get(i - 1) : null;
                final SvBreakend nextBreakend = (i < breakendCount-1) ? breakendList.get(i + 1) : null;

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

                // checkFoldbackSvs(breakend, prevBreakend);
            }
        }
    }

    public void reportPotentialFoldbacks(final String sampleId, final List<SvCluster> clusters)
    {
        if(!mFoldbackColsLogged)
        {
            mFoldbackColsLogged = true;
            String colNames = "FOLDBACK,Time,SampleId";
            colNames += ",Id1,Position1,Orientation1,SE1,Type1,Ploidy1,CopyNum1,CopyNumChg1,AsmbData1";
            colNames += ",Id2,Position2,Orientation2,SE2,Type2,Ploidy2,CopyNum2,CopyNumChg2,AsmbData2";
            colNames += ",Length,CopyNumDiff,EndsAssembled,IsLine,SameCluster";

            LOGGER.info(colNames);
        }

        for(Map.Entry<String, List<SvBreakend>> entry : mChrBreakendMap.entrySet())
        {
            List<SvBreakend> breakendList = entry.getValue();

            for(int i = 0; i < breakendList.size(); ++i)
            {
                final SvBreakend breakend = breakendList.get(i);
                final SvBreakend prevBreakend = (i > 0) ? breakendList.get(i - 1) : null;

                checkFoldbackSvs(sampleId, breakend, prevBreakend, clusters);
            }
        }
    }

    private void checkFoldbackSvs(final String sampleId, SvBreakend be1, SvBreakend be2, final List<SvCluster> clusters)
    {
        if(be1 == null || be2 == null)
            return;

        if(be1.orientation() != be2.orientation())
            return;

        final SvVarData var1 = be1.getSV();
        final SvVarData var2 = be2.getSV();

        if(var1.type() == INS || var2.type() == INS)
            return;

        // skip unclustered DELs & DUPs, reciprocal INV or reciprocal BNDs
        final SvCluster cluster1 = findCluster(var1, clusters);

        if(cluster1.isSimpleSVs() || (cluster1.getCount() == 2 && cluster1.isConsistent()))
            return;

        final SvCluster cluster2 = findCluster(var2, clusters);

        if(cluster2.isSimpleSVs() || (cluster2.getCount() == 2 && cluster2.isConsistent()))
            return;

        boolean v1Start = be1.usesStart();
        boolean v2Start = be2.usesStart();

        int length = (int)abs(be1.position() - be2.position());

        double copyNumberDiff;
        if((be1.orientation() == 1 && be1.position() < be2.position())
        || (be1.orientation() == -1 && be1.position() > be2.position()))
        {
            double cn1 = var1.copyNumber(v1Start) - var1.copyNumberChange(v1Start);
            double cn2 = var2.copyNumber(v2Start);
            copyNumberDiff = cn2 - cn1;
        }
        else
        {
            double cn2 = var2.copyNumber(v2Start) - var2.copyNumberChange(v2Start);
            double cn1 = var1.copyNumber(v1Start);
            copyNumberDiff = cn1 - cn2;
        }

        boolean otherEndsAssembled = haveLinkedAssemblies(var1, var2, !v1Start, !v2Start);
        boolean isLineElement = var1.isLineElement(v1Start) || var2.isLineElement(v2Start);
        boolean sameCluster = cluster1 == cluster2;

        // Add ploidy of both ends
        // Add assembly field for both ends

        LOGGER.info(String.format("FOLDBACK,%s,%s,%d,%d,%s,%s,%.3f,%.3f,%.3f,%s,%s,%d,%d,%s,%s,%.3f,%.3f,%.3f,%s,%d,%.3f,%s,%s,%s",
                sampleId,
                var1.id(), be1.position(), be1.orientation(), v1Start ? "start" : "end", var1.typeStr(),
                var1.getSvData().ploidy(), var1.copyNumber(v1Start), var1.copyNumberChange(v1Start), var1.getAssemblyData(v1Start),
                var2.id(), be2.position(), be2.orientation(), v2Start ? "start" : "end", var2.typeStr(),
                var2.getSvData().ploidy(), var2.copyNumber(v2Start), var2.copyNumberChange(v2Start), var2.getAssemblyData(v2Start),
                length, copyNumberDiff, otherEndsAssembled, isLineElement, sameCluster));
    }

    public void logInversionPairData(final String sampleId, final List<SvCluster> clusters)
    {
        if (!mInversionPairColsLogged)
        {
            mInversionPairColsLogged = true;
            String colNames = "INV_PAIRS,Time,SampleId,Id1,Chromsome,Orient1,PosStart1,PosEnd1";
            colNames += ",Id2,Orient2,PosStart2,PosEnd2,PairType,PairDesc,DelDupLen,TIAssembly,TILen";

            LOGGER.info(colNames);
        }

        for (final SvCluster cluster : clusters)
        {
            if (cluster.getUniqueSvCount() != 2 || cluster.getTypeCount(INV) != 2)
                continue;

            if (!cluster.isConsistent())
                continue;

            markInversionPairTypes(cluster, true, sampleId);
        }

    }

    public void markInversionPairTypes(SvCluster cluster, boolean logData, final String sampleId)
    {
        if(cluster.getLinkedPairs().isEmpty())
        {
            LOGGER.warn("cluster({}) missing creating linked pairs");
            return;
        }

        final SvVarData var1 = cluster.getSVs().get(0);
        final SvVarData var2 = cluster.getSVs().get(1);

        // determine overlap configurations

        /* 4 types of inversion pairs
        1. DEL with enclosed inverted TI (also know as 'Reciprocal INV') - Have 2 DSB and a ’TI' from the middle which is inverted.
            - outer breakends face out (the DEL)
            - TI enclosed
        2. DEL with external inverted TI
            - resultant type = DEL
            - length = other 2 breakends
        3. DUP with external inverted TI
            - 2 x TIs, but TI breakends don't overlap
            - type = DUP
            - TI and DUP length are interchangable, but choose shorter for TI
        4. DUP with enclosed inverted TI
            - no overlapping breakends
            - TI from the innermost 2
            - outer breakends face in
            - resultant type = DUP
            - length is outside 2 breakends (ie the other 2)
         */

        // first work out if there are 1 or 2 templated insertions

        final SvLinkedPair linkedPair1 = cluster.getLinkedPairs().get(0);

        boolean v1OpenOnStart = linkedPair1.first().equals(var1) ? linkedPair1.firstUnlinkedOnStart() : linkedPair1.secondUnlinkedOnStart();
        boolean v2OpenOnStart = linkedPair1.first().equals(var2) ? linkedPair1.firstUnlinkedOnStart() : linkedPair1.secondUnlinkedOnStart();

        SvLinkedPair linkedPair2 = cluster.getLinkedPairs().size() == 2 ? cluster.getLinkedPairs().get(1) : null;

        if(linkedPair2 != null && !linkedPair1.hasLinkClash(linkedPair2))
        {
            // existing second link is fine to consider
        }
        else if(areLinkedSection(var1, var2, v1OpenOnStart, v2OpenOnStart, false))
        {
            // could one be formed from the other 2 breakends
            linkedPair2 = new SvLinkedPair(var1, var2, LINK_TYPE_TI, v1OpenOnStart, v2OpenOnStart);
        }
        else if(areSectionBreak(var1, var2, v1OpenOnStart, v2OpenOnStart))
        {
            linkedPair2 = new SvLinkedPair(var1, var2, LINK_TYPE_DB, v1OpenOnStart, v2OpenOnStart);
        }
        else
        {
            LOGGER.error("sample({}) ids({} & {}) must be one or the other", sampleId, var1.id(), var2.id());
            return;
        }

        // set the first link to be the longer of the 2 for classification purposes
        SvLinkedPair lp1;
        SvLinkedPair lp2;

        if(linkedPair1.linkType() == LINK_TYPE_TI && linkedPair2.linkType() == LINK_TYPE_TI)
        {
            lp1 = linkedPair1.length() > linkedPair2.length() ? linkedPair1 : linkedPair2;
            lp2 = linkedPair1 == lp1 ? linkedPair2 : linkedPair1;
        }
        else if(linkedPair1.linkType() == LINK_TYPE_TI)
        {
            // take the DB for the main DEL if one exists
            lp2 = linkedPair1;
            lp1 = linkedPair2;
        }
        else if(linkedPair2.linkType() == LINK_TYPE_TI)
        {
            lp1 = linkedPair1;
            lp2 = linkedPair2;
        }
        else
        {
            // 2 deletion bridges
            lp1 = linkedPair1.length() > linkedPair2.length() ? linkedPair1 : linkedPair2;
            lp2 = linkedPair1 == lp1 ? linkedPair2 : linkedPair1;
        }

        long mainPos1 = lp1.first().position(lp1.firstLinkOnStart());
        long mainPos2 = lp1.second().position(lp1.secondLinkOnStart());

        long otherPos1 = lp2.first().position(lp2.firstLinkOnStart());
        long otherPos2 = lp2.second().position(lp2.secondLinkOnStart());

        boolean isTIEnclosed = false;

        long lowerMainPosLimit = mainPos1 < mainPos2 ? mainPos1 : mainPos2;
        lowerMainPosLimit -= MIN_TEMPLATED_INSERTION_LENGTH;

        long upperMainPosLimit = mainPos1 > mainPos2 ? mainPos1 : mainPos2;
        upperMainPosLimit += MIN_TEMPLATED_INSERTION_LENGTH;

        if(otherPos1 >= lowerMainPosLimit && otherPos1 <= upperMainPosLimit
        && otherPos2 >= lowerMainPosLimit && otherPos2 <= upperMainPosLimit)
        {
            isTIEnclosed = true;
        }

        StructuralVariantType invPairType;
        String invPairDesc = "";

        if(lp1.linkType() == LINK_TYPE_DB && lp2 != null && lp2.linkType() == LINK_TYPE_DB)
        {
            invPairType = DEL;
            invPairDesc = RESOLVED_TYPE_DEL_INT_TI;
        }
        else
        {
            if (lp1.linkType() == LINK_TYPE_DB || lp2.linkType() == LINK_TYPE_DB)
            {
                invPairType = DEL;

                if (!isTIEnclosed)
                {
                    invPairDesc = RESOLVED_TYPE_DEL_EXT_TI;
                }
                else
                {
                    invPairDesc = RESOLVED_TYPE_DEL_INT_TI;
                }
            }
            else
            {
                invPairType = DUP;

                if (!isTIEnclosed)
                {
                    invPairDesc = RESOLVED_TYPE_DUP_EXT_TI;
                }
                else
                {
                    invPairDesc = RESOLVED_TYPE_DUP_INT_TI;
                }
            }
        }

        long otherProximity = min(min(abs(mainPos1 - otherPos1), abs(mainPos1 - otherPos2)), min(abs(mainPos2 - otherPos1), abs(mainPos2 - otherPos2)));

        // check if there are variants falling within the outer breakends, indicating that this pair
        // are possibly not actually an isolated cluster 2
        boolean overlapsOtherVariants = hasEnclosedBreakends(
                var1.chromosome(true),
                min(var1.position(true), var2.position(true)),
                max(var1.position(false), var2.position(false)),
                cluster.getSVs());

        if(logData)
        {
            LOGGER.info(String.format("CL2_PAIR,%s,%s,%s,%s,%d,%d,%d,%s,%d,%d,%d,%s,%s,%d,%s,%d,%s,%d,%.2f,%.2f",
                    sampleId, cluster.getDesc(), var1.id(), var1.chromosome(true), var1.orientation(true), var1.position(true), var1.position(false),
                    var2.id(), var2.orientation(true), var2.position(true), var2.position(false),
                    invPairType, invPairDesc, lp1.length(), !lp1.isInferred(), lp2.length(), overlapsOtherVariants, otherProximity,
                    var1.getSvData().ploidy(), var2.getSvData().ploidy()));
        }

        cluster.setResolved(true, invPairDesc);
        cluster.setLengthOverride(lp1.length());
    }

    private boolean hasEnclosedBreakends(final String chromosome, long startPosition, long endPosition, List<SvVarData> skipList)
    {
        final List<SvBreakend> breakendList = mChrBreakendMap.get(chromosome);

        if(breakendList == null)
            return true;

        for(final SvBreakend breakend : breakendList)
        {
            if(skipList.contains(breakend.getSV()))
                continue;

            if(breakend.position() > startPosition && breakend.position() < endPosition)
                return true;
        }

        return false;
    }

    public void markBndPairTypes(SvCluster cluster, boolean logData, final String sampleId)
    {
        if(cluster.getLinkedPairs().isEmpty())
        {
            cluster.setResolved(true, RESOLVED_TYPE_RECIPROCAL_TRANS);
            return;
        }

        if(cluster.getLinkedPairs().size() == 2 && cluster.getArmGroups().size() == 2
        && cluster.getLinkedPairs().get(0).linkType() == LINK_TYPE_DB
        && cluster.getLinkedPairs().get(1).linkType() == LINK_TYPE_DB)
        {
            cluster.setResolved(true, RESOLVED_TYPE_RECIPROCAL_TRANS);
            return;
        }

        final SvVarData var1 = cluster.getSVs().get(0);
        final SvVarData var2 = cluster.getSVs().get(1);

        // first work out if there are 1 or 2 templated insertions
        final SvLinkedPair lp1;
        if(cluster.getLinkedPairs().size() > 0 && cluster.getLinkedPairs().get(0).linkType() == LINK_TYPE_TI)
            lp1 = cluster.getLinkedPairs().get(0);
        else
            lp1 = null;

        final SvLinkedPair lp2;
        if(cluster.getLinkedPairs().size() > 1 && cluster.getLinkedPairs().get(1).linkType() == LINK_TYPE_TI)
            lp2 = cluster.getLinkedPairs().get(1);
        else
            lp2 = null;

        if(lp1 == null && lp2 == null)
        {
            LOGGER.warn("cluster({} {}) has no linked pairs", cluster.getId(), cluster.getDesc());
            return;
        }

        final SvLinkedPair tiLinkedPair;

        if(lp1 != null && lp2 != null)
        {
            // take assembly pair if exists over inferred pair
            // otherwise take the shorter of the 2
            if(!lp1.isInferred())
                tiLinkedPair = lp1;
            else if(!lp2.isInferred())
                tiLinkedPair = lp2;
            else if(lp1.length() < lp2.length())
                tiLinkedPair = lp1;
            else
                tiLinkedPair = lp2;
        }
        else
        {
            tiLinkedPair = lp1 != null ? lp1 : lp2;
        }

        // must start and finish on the same chromosome
        boolean v1OpenOnStart = tiLinkedPair.first().equals(var1) ? tiLinkedPair.firstUnlinkedOnStart() : tiLinkedPair.secondUnlinkedOnStart();
        boolean v2OpenOnStart = tiLinkedPair.first().equals(var2) ? tiLinkedPair.firstUnlinkedOnStart() : tiLinkedPair.secondUnlinkedOnStart();

        if(!var1.chromosome(v1OpenOnStart).equals(var2.chromosome(v2OpenOnStart)))
            return;

        StructuralVariantType pairType;
        String pairDesc = "";

        if(areSectionBreak(var1, var2, v1OpenOnStart, v2OpenOnStart))
        {
            pairType = DEL;
            pairDesc = RESOLVED_TYPE_DEL_EXT_TI;
        }
        else
        {
            pairType = DUP;
            pairDesc = RESOLVED_TYPE_DUP_EXT_TI;
        }

        long pairTypeLength = abs(var1.position(v1OpenOnStart) - var2.position(v2OpenOnStart));

        boolean overlapsOtherVariants = hasEnclosedBreakends(
                var1.chromosome(v1OpenOnStart),
                min(var1.position(v1OpenOnStart), var2.position(v2OpenOnStart)),
                max(var1.position(v1OpenOnStart), var2.position(v2OpenOnStart)),
                cluster.getSVs());

        if(logData)
        {
            LOGGER.info(String.format("CL2_PAIR,%s,%s,%s,%s,%d,%d,%d,%s,%d,%d,%d,%s,%s,%d,%s,%d,%s,0,%.2f,%.2f",
                    sampleId, cluster.getDesc(), var1.id(), var1.chromosome(v1OpenOnStart), var1.orientation(true), var1.position(true), var1.position(false),
                    var2.id(), var2.orientation(true), var2.position(true), var2.position(false),
                    pairType, pairDesc, pairTypeLength, !tiLinkedPair.isInferred(), tiLinkedPair.length(), overlapsOtherVariants,
                    var1.getSvData().ploidy(), var2.getSvData().ploidy()));
        }

        cluster.setResolved(true, pairDesc);
        cluster.setLengthOverride(tiLinkedPair.length());

    }

    public boolean markDelDupPairTypes(SvCluster cluster, boolean logData, final String sampleId)
    {
        final SvVarData var1 = cluster.getSVs().get(0);
        final SvVarData var2 = cluster.getSVs().get(1);

        // determine overlap configurations

        /* 2 types of DEL and/or DUP pairs
        1. DEL with external TI
            - 2 DELs must be non-overlapping, forming a TI between them from innermost breakends OR
            - DEL-DUP overlapping with TI formed from right-most breakends
        2. DUP with external TI
            - DUP-DEL with DEL enclosed by DUP, TI formed from right-most breakends
            - length = other 2 breakends
        3. DUP with internal TI
            - either 2 DUPs overlapping or enclosed, TI formed from inner-most breakends OR
            - DUP-DEL with DEL enclosed by DUP, TI formed from right-most breakends
            - length = other 2 breakends
            - TI and DUP length are interchangable, but choose shorter for TI
        4. DEL with internal TI
            - from non-overlapping DUPs
            - TI formed from 1st and 3rd breakends and is longer than the DEL length
         */

        // first work out if there are 1 or 2 templated insertions
        final SvLinkedPair lp1;
        if(cluster.getLinkedPairs().size() > 0 && cluster.getLinkedPairs().get(0).linkType() == LINK_TYPE_TI)
            lp1 = cluster.getLinkedPairs().get(0);
        else
            lp1 = null;

        final SvLinkedPair lp2;
        if(cluster.getLinkedPairs().size() > 1 && cluster.getLinkedPairs().get(1).linkType() == LINK_TYPE_TI)
            lp2 = cluster.getLinkedPairs().get(1);
        else
            lp2 = null;

        if(lp1 == null && lp2 == null)
        {
            LOGGER.warn("cluster({} {}) has no linked pairs", cluster.getId(), cluster.getDesc());
            return false;
        }

        final SvLinkedPair tiLinkedPair;

        if(lp1 != null && lp2 != null)
        {
            // take assembly pair if exists over inferred pair
            // otherwise take the shorter of the 2
            if(!lp1.isInferred())
                tiLinkedPair = lp1;
            else if(!lp2.isInferred())
                tiLinkedPair = lp2;
            else if(lp1.length() < lp2.length())
                tiLinkedPair = lp1;
            else
                tiLinkedPair = lp2;
        }
        else
        {
            tiLinkedPair = lp1 != null ? lp1 : lp2;
        }

        StructuralVariantType pairType;
        String pairDesc = "";

        long tiPos1 = tiLinkedPair.first().position(tiLinkedPair.firstLinkOnStart());
        long tiPos2 = tiLinkedPair.second().position(tiLinkedPair.secondLinkOnStart());

        long otherPos1 = tiLinkedPair.first().position(!tiLinkedPair.firstLinkOnStart());
        long otherPos2 = tiLinkedPair.second().position(!tiLinkedPair.secondLinkOnStart());

        long tiUpperPos = tiPos1 < tiPos2 ? tiPos2 : tiPos1;
        long tiLowerPos = tiPos2 < tiPos1 ? tiPos2 : tiPos1;

        long nonTiUpperPos = otherPos1 < otherPos2 ? otherPos2 : otherPos1;
        long nonTiLowerPos = otherPos2 < otherPos1 ? otherPos2 : otherPos1;

        // boolean isTIEnclosed = tiUpperPos < nonTiUpperPos && tiLowerPos > nonTiLowerPos;
        boolean isTIEnclosing = tiUpperPos > nonTiUpperPos && tiLowerPos < nonTiLowerPos;

        boolean noOverlap = var1.position(false) < var2.position(true) || var2.position(false) < var1.position(true);

        long pairTypeLength = abs(otherPos1 - otherPos2);

        boolean isValid = true;

        if(var1.type() == DEL && var2.type() == DEL)
        {
            // must be non-overlapping otherwise wouldn't be clustered
            pairType = DEL;
            pairDesc = RESOLVED_TYPE_DEL_INT_TI;

            isValid = noOverlap;
        }
        else if(var1.type() == DUP && var2.type() == DUP)
        {
            if(noOverlap)
            {
                pairType = DEL;
                pairDesc = RESOLVED_TYPE_DEL_INT_TI;

                isValid = isTIEnclosing;
            }
            else
            {
                pairType = DUP;
                pairDesc = RESOLVED_TYPE_DUP_INT_TI;

                // TI can straddle other breakends
                // isValid = isTIEnclosed;
            }
        }
        else
        {
            if((var1.type() == DEL && var1.position(true) < var2.position(true))
            || (var2.type() == DEL && var2.position(true) < var1.position(true)))
            {
                // DEL first
                pairType = DEL;
                pairDesc = RESOLVED_TYPE_DEL_EXT_TI;

                // isValid = tiLowerPos > nonTiUpperPos;
            }
            else
            {
                pairType = DUP;
                pairDesc = RESOLVED_TYPE_DUP_EXT_TI;

                // isValid = tiLowerPos > nonTiUpperPos;
            }
        }

        boolean overlapsOtherVariants = hasEnclosedBreakends(
                var1.chromosome(true),
                min(var1.position(true), var2.position(true)),
                max(var1.position(false), var2.position(false)),
                cluster.getSVs());

        long otherProximity = min(min(abs(tiPos1 - otherPos1), abs(tiPos1 - otherPos2)), min(abs(tiPos2 - otherPos1), abs(tiPos2 - otherPos2)));

        if(!isValid)
        {
            LOGGER.debug("cluster({}) invalid classification", cluster.getId());
            return false;
        }

        if(logData)
        {
            LOGGER.info(String.format("CL2_PAIR,%s,%s,%s,%s,%d,%d,%d,%s,%d,%d,%d,%s,%s,%d,%s,%d,%s,%d,%.2f,%.2f",
                    sampleId, cluster.getDesc(), var1.id(), var1.chromosome(true), var1.orientation(true), var1.position(true), var1.position(false),
                    var2.id(), var2.orientation(true), var2.position(true), var2.position(false),
                    pairType, pairDesc, pairTypeLength, !tiLinkedPair.isInferred(), tiLinkedPair.length(), overlapsOtherVariants, otherProximity,
                    var1.getSvData().ploidy(), var2.getSvData().ploidy()));
        }

        cluster.setResolved(true, pairDesc);
        cluster.setLengthOverride(tiLinkedPair.length());
        return true;
    }

    public void setChromosomalArmStats(final List<SvVarData> allVariants)
    {
        mChrArmSvCount.clear();
        mChrArmSvExpected.clear();
        mChrArmSvRate.clear();
        mMedianChrArmRate = 0;

        // form a map of unique arm to SV count
        for(final SvVarData var : allVariants)
        {
            String chrArmStart = mUtils.getVariantChrArm(var,true);
            String chrArmEnd = mUtils.getVariantChrArm(var,false);

            // ensure an entry exists
            if (!mChrArmSvCount.containsKey(chrArmStart))
            {
                mChrArmSvCount.put(chrArmStart, 0);
            }

            // exclude LINE elements from back-ground rates
            if(var.isLineElement(true))
            {
                mChrArmSvCount.replace(chrArmStart, mChrArmSvCount.get(chrArmStart) + 1);
            }

            if(!var.isNullBreakend())
            {
                if (!chrArmStart.equals(chrArmEnd) && !mChrArmSvCount.containsKey(chrArmEnd))
                {
                    mChrArmSvCount.put(chrArmEnd, 0);
                }

                if (var.isLineElement(false))
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

    public String getChrArmData(final SvVarData var)
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
