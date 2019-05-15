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
import static com.hartwig.hmftools.svanalysis.analysis.LinkFinder.MIN_TEMPLATED_INSERTION_LENGTH;
import static com.hartwig.hmftools.svanalysis.analysis.LinkFinder.areLinkedSection;
import static com.hartwig.hmftools.svanalysis.analysis.LinkFinder.arePairedDeletionBridges;
import static com.hartwig.hmftools.svanalysis.analysis.LinkFinder.areSectionBreak;
import static com.hartwig.hmftools.svanalysis.analysis.LinkFinder.getMinTemplatedInsertionLength;
import static com.hartwig.hmftools.svanalysis.analysis.SvUtilities.CHROMOSOME_ARM_P;
import static com.hartwig.hmftools.svanalysis.analysis.SvUtilities.CHROMOSOME_ARM_Q;
import static com.hartwig.hmftools.svanalysis.analysis.SvUtilities.PERMITED_DUP_BE_DISTANCE;
import static com.hartwig.hmftools.svanalysis.analysis.SvUtilities.addSvToChrBreakendMap;
import static com.hartwig.hmftools.svanalysis.analysis.SvUtilities.copyNumbersEqual;
import static com.hartwig.hmftools.svanalysis.types.SvCluster.RESOLVED_TYPE_DUP_BE;
import static com.hartwig.hmftools.svanalysis.types.SvCluster.RESOLVED_TYPE_LOW_CNC;
import static com.hartwig.hmftools.svanalysis.types.SvCluster.RESOLVED_TYPE_DEL_EXT_TI;
import static com.hartwig.hmftools.svanalysis.types.SvCluster.RESOLVED_TYPE_DEL_INT_TI;
import static com.hartwig.hmftools.svanalysis.types.SvCluster.RESOLVED_TYPE_DUP_EXT_TI;
import static com.hartwig.hmftools.svanalysis.types.SvCluster.RESOLVED_TYPE_DUP_INT_TI;
import static com.hartwig.hmftools.svanalysis.types.SvCluster.RESOLVED_TYPE_LINE;
import static com.hartwig.hmftools.svanalysis.types.SvCluster.RESOLVED_TYPE_NONE;
import static com.hartwig.hmftools.svanalysis.types.SvCluster.RESOLVED_TYPE_POLY_G_C;
import static com.hartwig.hmftools.svanalysis.types.SvCluster.RESOLVED_TYPE_RECIPROCAL_TRANS;
import static com.hartwig.hmftools.svanalysis.types.SvCluster.RESOLVED_TYPE_SGL_PAIR_DEL;
import static com.hartwig.hmftools.svanalysis.types.SvCluster.RESOLVED_TYPE_SGL_PAIR_DUP;
import static com.hartwig.hmftools.svanalysis.types.SvCluster.RESOLVED_TYPE_SGL_PAIR_INS;
import static com.hartwig.hmftools.svanalysis.types.SvCluster.RESOLVED_TYPE_SGL_PLUS_INCONSISTENT;
import static com.hartwig.hmftools.svanalysis.types.SvCluster.RESOLVED_TYPE_SIMPLE_SV;
import static com.hartwig.hmftools.svanalysis.types.SvLOH.LOH_NO_SV;
import static com.hartwig.hmftools.svanalysis.types.SvLinkedPair.LINK_TYPE_DB;
import static com.hartwig.hmftools.svanalysis.types.SvLinkedPair.LINK_TYPE_TI;
import static com.hartwig.hmftools.svanalysis.types.SvVarData.ASSEMBLY_TYPE_EQV;
import static com.hartwig.hmftools.svanalysis.types.SvVarData.RELATION_TYPE_NEIGHBOUR;
import static com.hartwig.hmftools.svanalysis.types.SvVarData.RELATION_TYPE_OVERLAP;

import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantType;
import com.hartwig.hmftools.svanalysis.types.SvBreakend;
import com.hartwig.hmftools.svanalysis.types.SvCluster;
import com.hartwig.hmftools.svanalysis.types.SvLinkedPair;
import com.hartwig.hmftools.svanalysis.types.SvVarData;
import com.hartwig.hmftools.svanalysis.types.SvLOH;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class SvClusteringMethods {

    private static final Logger LOGGER = LogManager.getLogger(SvClusteringMethods.class);

    private int mNextClusterId;

    private Map<String, List<SvBreakend>> mChrBreakendMap; // every breakend on a chromosome, ordered by ascending position
    private Map<String, List<SvLOH>> mSampleLohData;

    private Map<String, double[]> mChromosomeCopyNumberMap; // p-arm telomere, centromere and q-arm telemore CN data

    private long mDelDupCutoffLength;
    private int mProximityDistance;

    public static int MAX_SIMPLE_DUP_DEL_CUTOFF = 5000000;
    public static int MIN_SIMPLE_DUP_DEL_CUTOFF = 100000;
    public static int DEFAULT_PROXIMITY_DISTANCE = 5000;

    public static String CLUSTER_REASON_PROXIMITY = "Prox";
    public static String CLUSTER_REASON_LOH = "LOH";
    public static String CLUSTER_REASON_COMMON_ARMS = "ComArm";
    public static String CLUSTER_REASON_FOLDBACKS = "Foldback";
    public static String CLUSTER_REASON_SOLO_SINGLE = "Single";
    public static String CLUSTER_REASON_LOOSE_OVERLAP = "LooseOverlap";
    public static String CLUSTER_REASON_LOH_CHAIN = "LOHChain";
    public static String CLUSTER_REASON_NET_ARM_END_PLOIDY = "ArmEndPloidy";
    public static String CLUSTER_REASON_BE_PLOIDY_DROP = "BEPloidy";
    public static String CLUSTER_REASON_LONG_DEL_DUP_OR_INV = "DelDupInv";

    public SvClusteringMethods(int proximityLength)
    {
        mNextClusterId = 0;

        mDelDupCutoffLength = 0;
        mProximityDistance = proximityLength;

        mChrBreakendMap = new HashMap();
        mChromosomeCopyNumberMap = new HashMap();
        mSampleLohData = null;
    }

    public final Map<String, List<SvBreakend>> getChrBreakendMap() { return mChrBreakendMap; }
    public final Map<String, List<SvLOH>> getSampleLohData() { return mSampleLohData; }
    public int getNextClusterId() { return mNextClusterId++; }
    public void setSampleLohData(final Map<String, List<SvLOH>> data) { mSampleLohData = data; }
    public void setChrCopyNumberMap(final Map<String, double[]> data) { mChromosomeCopyNumberMap = data; }
    public final Map<String, double[]> getChrCopyNumberMap() { return mChromosomeCopyNumberMap; }
    public long getDelDupCutoffLength() { return mDelDupCutoffLength; }
    public int getProximityDistance() { return mProximityDistance; }

    public static boolean isFilteredResolvedType(final String resolvedType)
    {
        return resolvedType.equals(RESOLVED_TYPE_LOW_CNC) || resolvedType.equals(RESOLVED_TYPE_POLY_G_C)
                || resolvedType.equals(RESOLVED_TYPE_DUP_BE);
    }

    public void clusterByProximity(List<SvVarData> allVariants, List<SvCluster> clusters)
    {
        clusterByProximity(clusters);
    }

    private void clusterByProximity(List<SvCluster> clusters)
    {
        mNextClusterId = 0;

        // walk through each chromosome and breakend list
        for (final Map.Entry<String, List<SvBreakend>> entry : mChrBreakendMap.entrySet())
        {
            final List<SvBreakend> breakendList = entry.getValue();

            int currentIndex = 0;
            while (currentIndex < breakendList.size())
            {
                final SvBreakend breakend = breakendList.get(currentIndex);
                SvVarData var = breakend.getSV();

                // isSpecificSV(var);

                String sglExcludedReason = getSingleBreakendUnclusteredType(var);
                if (sglExcludedReason != RESOLVED_TYPE_NONE)
                {
                    if(var.getCluster() == null)
                    {
                        SvCluster newCluster = new SvCluster(getNextClusterId());
                        newCluster.addVariant(var);
                        newCluster.setResolved(true, sglExcludedReason);
                        clusters.add(newCluster);
                    }
                    ++currentIndex;
                    continue;
                }

                SvBreakend nextBreakend = null;
                int nextIndex = currentIndex + 1;

                for (; nextIndex < breakendList.size(); ++nextIndex)
                {
                    final SvBreakend nextBe = breakendList.get(nextIndex);
                    final SvVarData nextVar = nextBe.getSV();

                    // isSpecificSV(nextVar);
                    sglExcludedReason = getSingleBreakendUnclusteredType(nextVar);
                    if (sglExcludedReason != RESOLVED_TYPE_NONE)
                    {
                        if(nextVar.getCluster() == null)
                        {
                            SvCluster newCluster = new SvCluster(getNextClusterId());
                            newCluster.addVariant(nextVar);
                            newCluster.setResolved(true, sglExcludedReason);
                            clusters.add(newCluster);
                        }
                        continue;
                    }

                    nextBreakend = nextBe;
                    break;
                }

                if (nextBreakend == null)
                {
                    // no more breakends on this chromosome
                    if (var.getCluster() == null)
                    {
                        SvCluster cluster = new SvCluster(getNextClusterId());
                        cluster.addVariant(var);
                        clusters.add(cluster);
                    }

                    break;
                }

                SvCluster cluster = var.getCluster();
                SvVarData nextVar = nextBreakend.getSV();
                SvCluster nextCluster = nextVar.getCluster();

                if (cluster != null && cluster == nextCluster)
                {
                    // already clustered
                }
                else if (abs(nextBreakend.position() - breakend.position()) > mProximityDistance)
                {
                    // too far between the breakends
                    if (cluster == null)
                    {
                        cluster = new SvCluster(getNextClusterId());
                        cluster.addVariant(var);
                        clusters.add(cluster);
                    }
                    else
                    {
                        // nothing more to do for this variant - already clustered
                    }
                }
                else
                {
                    // 2 breakends are close enough to cluster
                    if (var == nextVar)
                    {
                        if (cluster == null)
                        {
                            cluster = new SvCluster(getNextClusterId());
                            cluster.addVariant(var);
                            clusters.add(cluster);
                        }
                    }
                    else
                    {
                        // one or both SVs could already be a part of clusters, or neither may be
                        if (cluster == null && nextCluster == null)
                        {
                            cluster = new SvCluster(getNextClusterId());
                            cluster.addVariant(var);
                            cluster.addVariant(nextVar);
                            cluster.addClusterReason(CLUSTER_REASON_PROXIMITY);
                            clusters.add(cluster);
                        }
                        else if (cluster != null && nextCluster != null)
                        {
                            // keep one and remove the other
                            cluster.mergeOtherCluster(nextCluster, false);
                            clusters.remove(nextCluster);
                        }
                        else
                        {
                            if (cluster == null)
                            {
                                nextCluster.addVariant(var);
                            }
                            else
                            {
                                cluster.addVariant(nextVar);
                            }
                        }

                        if (var.getClusterReason().isEmpty())
                            var.addClusterReason(CLUSTER_REASON_PROXIMITY, nextVar.id());

                        if (nextVar.getClusterReason().isEmpty())
                            nextVar.addClusterReason(CLUSTER_REASON_PROXIMITY, var.id());
                    }

                }

                // move the index to the SV which was just proximity cluster so the next comparison is with the closest candidate
                currentIndex = nextIndex;
            }
        }

        splitDelClusters(clusters);
    }

    private void splitDelClusters(List<SvCluster> clusters)
    {
        // no 2 DELs can overlap, so if the cluster is entirely comprised of DELs, split out any overlapping SVs
        List<SvCluster> clustersToRemove = Lists.newArrayList();

        int clusterCount = clusters.size();
        for(int i = 0; i < clusterCount; ++i)
        {
            SvCluster cluster = clusters.get(i);

            if(cluster.getSvCount() == 1 || cluster.getSvCount() != cluster.getTypeCount(DEL))
                continue;

            List<SvVarData> delSVs = Lists.newArrayList(cluster.getSVs());
            List<SvCluster> newClusters = Lists.newArrayList();

            for(SvVarData var : delSVs)
            {
                boolean addedToCluster = false;
                for(final SvCluster newCluster : newClusters)
                {
                    boolean hasOverlap = false;

                    for(final SvVarData existingVar : newCluster.getSVs())
                    {
                        if(var.position(true) >= existingVar.position(false)
                        || var.position(false) <= existingVar.position(true))
                        {
                            continue;
                        }

                        // has an overlap
                        hasOverlap = true;
                        break;
                    }

                    if(!hasOverlap)
                    {
                        newCluster.addVariant(var);
                        addedToCluster = true;
                        break;
                    }
                }

                if(!addedToCluster)
                {
                    SvCluster newCluster = new SvCluster(getNextClusterId());
                    newClusters.add(newCluster);
                    newCluster.addVariant(var);
                    var.clearClusterReason();
                }
            }

            clustersToRemove.add(cluster);
            clusters.addAll(newClusters);
        }

        for(SvCluster cluster : clustersToRemove)
        {
            clusters.remove(cluster);
        }
    }

    public boolean validateClustering(final List<SvCluster> clusters)
    {
        // validation that every SV was put into a cluster
        for (final Map.Entry<String, List<SvBreakend>> entry : mChrBreakendMap.entrySet())
        {
            final List<SvBreakend> breakendList = entry.getValue();

            for (int i = 0; i < breakendList.size(); ++i)
            {
                final SvBreakend breakend = breakendList.get(i);
                SvVarData var = breakend.getSV();
                if(var.getCluster() == null)
                {
                    LOGGER.error("var({}) not clustered", var.posId());
                    return false;
                }
            }
        }

        // check that no 2 clusters contain the same SV
        for(int i = 0; i < clusters.size(); ++i)
        {
            SvCluster cluster1 = clusters.get(i);
            // isSpecificCluster(cluster1);

            // check all SVs in this cluster reference it
            for(SvVarData var : cluster1.getSVs())
            {
                if(var.getCluster() != cluster1)
                {
                    LOGGER.error("var({}) in cluster({}) has incorrect ref", var.posId(), cluster1.id());
                    return false;
                }
            }

            for(int j = i+1; j < clusters.size(); ++j)
            {
                SvCluster cluster2 = clusters.get(j);

                for(SvVarData var : cluster1.getSVs())
                {
                    if(cluster2.getSVs().contains(var))
                    {
                        LOGGER.error("var({}) in 2 clusters({} and {})", var.posId(), cluster1.id(), cluster2.id());
                        return false;
                    }
                }
            }
        }

        return true;
    }

    public static boolean checkClusterDuplicates(List<SvCluster> clusters)
    {
        for(int i = 0; i < clusters.size(); ++i)
        {
            final SvCluster cluster1 = clusters.get(i);
            for(int j = i + 1; j < clusters.size(); ++j)
            {
                final SvCluster cluster2 = clusters.get(j);

                if(cluster1 == cluster2 || cluster1.id() == cluster2.id())
                {
                    LOGGER.error("cluster({}) exists twice in list", cluster1.id());
                    return false;
                }
            }
        }

        return true;
    }

    public void mergeClusters(final String sampleId, List<SvCluster> clusters)
    {
        // first apply replication rules since this can affect consistency
        for(SvCluster cluster : clusters)
        {
            // check for sub-clonal / low-copy-number supported variants
            if(cluster.getSvCount() == 1 && hasLowCNChangeSupport(cluster.getSV(0)))
            {
                cluster.setResolved(true, RESOLVED_TYPE_LOW_CNC);
                continue;
            }

            if(cluster.hasLinkingLineElements())
                continue;

            markClusterLongDelDups(cluster);
            markClusterInversions(cluster);
        }

        int initClusterCount = clusters.size();

        int iterations = 0;

        // the merge must be run a few times since as clusters grow, more single SVs and other clusters
        // will then fall within the bounds of the new larger clusters
        boolean foundMerges = true;

        while(foundMerges && iterations < 5)
        {
            if(mergeOnOverlappingInvDupDels(clusters))
                foundMerges = true;

            if(mergeOnUnresolvedSingles(clusters))
                foundMerges = true;

            ++iterations;
        }

        mergeOnLOHEvents(sampleId, clusters);

        if(clusters.size() < initClusterCount)
        {
            LOGGER.debug("reduced cluster count({} -> {}) iterations({})", initClusterCount, clusters.size(), iterations);
        }
    }

    public static int MAX_SV_REPLICATION_MULTIPLE = 32;
    public static int MAX_CLUSTER_COUNT_REPLICATION = 500;

    public static void applyCopyNumberReplication(SvCluster cluster)
    {
        // isSpecificCluster(cluster);

        // use the relative copy number change to replicate some SVs within a cluster
        if(!cluster.hasVariedCopyNumber())
            return;

        // avoid what will likely be balanced translocations
        if(cluster.getSvCount() == 2 && cluster.getTypeCount(BND) == 2)
            return;

        int maxReplication = cluster.getSvCount() > MAX_CLUSTER_COUNT_REPLICATION ? 8 : MAX_SV_REPLICATION_MULTIPLE;

        // first establish the lowest copy number change
        double minCopyNumber = cluster.getMinCNChange();
        double maxCopyNumber = cluster.getMaxCNChange();

        if(minCopyNumber <= 0)
        {
            LOGGER.debug("cluster({}) warning: invalid CN variation(min={} max={})",
                    cluster.id(), minCopyNumber, maxCopyNumber);
            return;
        }

        double replicationFactor = 1;
        if(maxCopyNumber > MAX_SV_REPLICATION_MULTIPLE * minCopyNumber)
        {
            LOGGER.debug("cluster({}) warning: large CN variation(min={} max={})",
                    cluster.id(), minCopyNumber, maxCopyNumber);

            // scale the replication down
            replicationFactor = maxReplication / (maxCopyNumber / minCopyNumber);
        }

        // replicate the SVs which have a higher copy number than their peers
        int clusterCount = cluster.getSvCount();

        for(int i = 0; i < clusterCount; ++i)
        {
            SvVarData var = cluster.getSV(i);
            double calcCopyNumber = var.getRoundedCNChange();

            int svMultiple = (int)round(calcCopyNumber / minCopyNumber);

            svMultiple = max((int)round(svMultiple * replicationFactor), 1);
            // svMultiple = min(svMultiple, MAX_SV_REPLICATION_MULTIPLE);

            if(svMultiple <= 1)
                continue;

            LOGGER.debug("cluster({}) replicating SV({}) {} times, copyNumChg({} vs min={})",
                    cluster.id(), var.posId(), svMultiple, calcCopyNumber, minCopyNumber);

            var.setReplicatedCount(svMultiple);

            for(int j = 1; j < svMultiple; ++j)
            {
                SvVarData newVar = new SvVarData(var);
                cluster.addVariant(newVar);
            }
        }
    }

    public static double LOW_QUALITY_CN_CHANGE = 0.5;

    public boolean hasLowCNChangeSupport(final SvVarData var)
    {
        if(var.type() == INS)
            return false;

        if(var.isNullBreakend())
            return var.copyNumberChange(true) < LOW_QUALITY_CN_CHANGE;
        else
            return var.copyNumberChange(true) < LOW_QUALITY_CN_CHANGE && var.copyNumberChange(false) < LOW_QUALITY_CN_CHANGE;
    }

    private static String POLY_C_MOTIF = "CCCCCCCCCCCCCCCC";
    private static String POLY_G_MOTIF = "GGGGGGGGGGGGGGGG";

    private String getSingleBreakendUnclusteredType(final SvVarData var)
    {
        if(var.type() != SGL)
            return RESOLVED_TYPE_NONE;

        if(isEquivSingleBreakend(var))
            return RESOLVED_TYPE_DUP_BE;

        if(var.getSvData().insertSequence().contains(POLY_C_MOTIF) || var.getSvData().insertSequence().contains(POLY_G_MOTIF))
            return RESOLVED_TYPE_POLY_G_C;

        return RESOLVED_TYPE_NONE;
    }

    private boolean isEquivSingleBreakend(final SvVarData var)
    {
        if(!var.isNullBreakend())
            return false;

        if(var.isDupBreakend(true) || var.isDupBreakend(false))
            return true;

        return  (var.getAssemblyData(true).contains(ASSEMBLY_TYPE_EQV)
                || var.getAssemblyData(false).contains(ASSEMBLY_TYPE_EQV));
    }

    public void clearLOHBreakendData(final String sampleId)
    {
        if(mSampleLohData == null || sampleId.isEmpty())
            return;

        List<SvLOH> lohList = mSampleLohData.get(sampleId);

        if(lohList == null)
            return;

        lohList = lohList.stream().filter(x -> !x.Skipped).collect(Collectors.toList());

        for(final SvLOH lohEvent : lohList)
        {
            lohEvent.setBreakend(null, true);
            lohEvent.setBreakend(null, false);
        }
    }

    private void mergeOnLOHEvents(final String sampleId, List<SvCluster> clusters)
    {
        if(mSampleLohData == null)
            return;

        // first extract all the SVs from the LOH events
        List<SvLOH> lohList = mSampleLohData.get(sampleId);

        if(lohList == null)
            return;

        // note that LOH-breakend links are established here and then must be

        lohList = lohList.stream().filter(x -> !x.Skipped).collect(Collectors.toList());

        String currentChromosome = "";
        List<SvBreakend> breakendList = null;

        for(final SvLOH lohEvent : lohList)
        {
            if((lohEvent.StartSV.isEmpty() && lohEvent.EndSV.isEmpty())
            || (lohEvent.StartSV.equals(LOH_NO_SV) && lohEvent.EndSV.equals(LOH_NO_SV)))
                continue;

            SvCluster lohClusterStart = null;
            SvVarData lohSvStart = null;
            SvCluster lohClusterEnd = null;
            SvVarData lohSvEnd = null;

            // use the breakend table to find matching SVs
            if(breakendList == null || !currentChromosome.equals(lohEvent.Chromosome))
            {
                breakendList = mChrBreakendMap.get(lohEvent.Chromosome);
                currentChromosome = lohEvent.Chromosome;
            }

            if(breakendList == null)
                continue;

            for(final SvBreakend breakend : breakendList)
            {
                if(breakend.orientation() == 1 && breakend.getSV().id().equals(lohEvent.StartSV))
                {
                    lohClusterStart = breakend.getSV().getCluster();
                    lohSvStart = breakend.getSV();
                    lohClusterStart.addLohEvent(lohEvent);
                    lohEvent.setBreakend(breakend, true);
                }

                if(breakend.orientation() == -1 && breakend.getSV().id().equals(lohEvent.EndSV))
                {
                    lohClusterEnd = breakend.getSV().getCluster();

                    if(!lohClusterEnd.getLohEvents().contains(lohEvent))
                        lohClusterEnd.addLohEvent(lohEvent);

                    lohSvEnd = breakend.getSV();
                    lohEvent.setBreakend(breakend, false);
                }

                if(lohEvent.matchedBothSVs())
                    break;
            }

            if(!lohEvent.IsValid)
                continue; // cannot be used for clustering

            if(lohClusterEnd == null || lohClusterStart == null)
            {
                // LOGGER.error("sample({}) start varId({}) not found in any cluster", sampleId, lohEvent.StartSV);
                continue;
            }

            if(lohClusterStart == lohClusterEnd)
                continue;

            if(lohClusterStart.hasLinkingLineElements() || lohClusterEnd.hasLinkingLineElements())
                continue;

            LOGGER.debug("cluster({} svs={}) merges in other cluster({} svs={}) on LOH event: SVs({} and {}) length({})",
                    lohClusterStart.id(), lohClusterStart.getSvCount(), lohClusterEnd.id(), lohClusterEnd.getSvCount(),
                    lohEvent.StartSV, lohEvent.EndSV, lohEvent.Length);

            lohSvStart.addClusterReason(CLUSTER_REASON_LOH, lohSvEnd.id());
            lohSvEnd.addClusterReason(CLUSTER_REASON_LOH, lohSvStart.id());
            lohClusterStart.addClusterReason(CLUSTER_REASON_LOH);

            lohClusterStart.mergeOtherCluster(lohClusterEnd);
            clusters.remove(lohClusterEnd);
        }
    }

    private void markClusterInversions(final SvCluster cluster)
    {
        if(cluster.getTypeCount(INV) == 0 || cluster.hasLinkingLineElements())
            return;

        // skip cluster-2s which resolved to a simple type
        if(cluster.isSyntheticSimpleType(true))
            return;

        for (final SvVarData var : cluster.getSVs())
        {
            if(var.type() == INV && !var.isCrossArm())
            {
                cluster.registerInversion(var);
            }
        }
    }

    private void markClusterLongDelDups(final SvCluster cluster)
    {
        if(cluster.hasLinkingLineElements())
            return;

        // isSpecificCluster(cluster);

        if(cluster.isSyntheticSimpleType(false))
        {
            if(cluster.getSynDelDupLength() >= mDelDupCutoffLength)
            {
                for (final SvVarData var : cluster.getSVs())
                {
                    cluster.registerLongDelDup(var);
                }
            }

            return;
        }

        if(cluster.getTypeCount(DEL) > 0 || cluster.getTypeCount(DUP) > 0)
        {
            for(final SvVarData var : cluster.getSVs())
            {
                if((var.type() == DUP || var.type() == DEL) && var.length() >= mDelDupCutoffLength && !var.isCrossArm())
                {
                    cluster.registerLongDelDup(var);
                }
            }
        }
    }

    private boolean mergeOnOverlappingInvDupDels(List<SvCluster> clusters)
    {
        // merge any clusters with overlapping inversions, long dels or long dups on the same arm
        int initClusterCount = clusters.size();

        final List<StructuralVariantType> requiredTypes = Lists.newArrayList();
        requiredTypes.add(INV);

        int index1 = 0;
        while(index1 < clusters.size())
        {
            SvCluster cluster1 = clusters.get(index1);

            if(cluster1.getInversions().isEmpty() && cluster1.getLongDelDups().isEmpty() || cluster1.hasLinkingLineElements())
            {
                ++index1;
                continue;
            }

            List<SvVarData> cluster1Svs = Lists.newArrayList();
            cluster1Svs.addAll(cluster1.getLongDelDups());
            cluster1Svs.addAll(cluster1.getInversions());

            int index2 = index1 + 1;
            while(index2 < clusters.size())
            {
                SvCluster cluster2 = clusters.get(index2);

                if(cluster2.getInversions().isEmpty() && cluster2.getLongDelDups().isEmpty() || cluster2.hasLinkingLineElements())
                {
                    ++index2;
                    continue;
                }

                List<SvVarData> cluster2Svs = Lists.newArrayList();
                cluster2Svs.addAll(cluster2.getLongDelDups());
                cluster2Svs.addAll(cluster2.getInversions());

                boolean canMergeClusters = false;
                String clusterReason = "";

                for (final SvVarData var1 : cluster1Svs)
                {
                    for (final SvVarData var2 : cluster2Svs)
                    {
                        if(!var1.chromosome(true).equals(var2.chromosome(true)))
                            continue;

                        if(var1.position(false) < var2.position(true) || var1.position(true) > var2.position(false))
                            continue;

                        LOGGER.debug("cluster({}) SV({} {}) and cluster({}) SV({} {}) have inversion or longDelDup overlap",
                                cluster1.id(), var1.posId(), var1.type(), cluster2.id(), var2.posId(), var2.type());

                        var1.addClusterReason(CLUSTER_REASON_LONG_DEL_DUP_OR_INV, var2.id());
                        var2.addClusterReason(CLUSTER_REASON_LONG_DEL_DUP_OR_INV, var1.id());

                        canMergeClusters = true;
                        break;
                    }


                    if(canMergeClusters)
                        break;
                }

                if(canMergeClusters)
                {
                    cluster1.mergeOtherCluster(cluster2);
                    cluster1.addClusterReason(clusterReason);
                    clusters.remove(index2);
                }
                else
                {
                    ++index2;
                }
            }

            ++index1;
        }

        return clusters.size() < initClusterCount;
    }

    private boolean mergeOnUnresolvedSingles(List<SvCluster> clusters)
    {
        // merge clusters with 1 unresolved single with the following rules:
        // 2 x cluster-1s with SGLs that are each other's nearest neighbours
        //
        // use the chr-breakend map to walk through and find the closest links
        // only apply a rule between the 2 closest breakends at the exclusions of the cluster on their other end
        // unless the other breakend is a short, simple SV

        boolean foundMerges = false;

        for (final Map.Entry<String, List<SvBreakend>> entry : mChrBreakendMap.entrySet())
        {
            final List<SvBreakend> breakendList = entry.getValue();
            int breakendCount = breakendList.size();

            for(int i = 0; i < breakendList.size(); ++i)
            {
                final SvBreakend breakend = breakendList.get(i);
                SvVarData var = breakend.getSV();
                final SvCluster cluster = var.getCluster();

                // take the point of view of the cluster with the solo single
                if(cluster.getSvCount() != 1 || cluster.getSV(0).type() != SGL || cluster.isResolved())
                    continue;

                // now look for a proximate cluster with either another solo single or needing one to be resolved
                // check previous and next breakend's cluster
                SvBreakend prevBreakend = (i > 0) ? breakendList.get(i - 1) : null;
                SvBreakend nextBreakend = (i < breakendCount - 1) ? breakendList.get(i + 1) : null;

                if(prevBreakend != null && prevBreakend.getSV().isSimpleType() && prevBreakend.getSV().length() < mDelDupCutoffLength)
                    prevBreakend = null;

                if(nextBreakend != null && nextBreakend.getSV().isSimpleType() && nextBreakend.getSV().length() < mDelDupCutoffLength)
                    nextBreakend = null;

                // additionally check that breakend after the next one isn't a closer SGL to the next breakend,
                // which would invalidate this one being the nearest neighbour
                long prevProximity = prevBreakend != null ? abs(breakend.position() - prevBreakend.position()) : -1;
                long nextProximity = nextBreakend != null ? abs(breakend.position() - nextBreakend.position()) : -1;

                if(nextBreakend != null && i < breakendCount - 2)
                {
                    final SvBreakend followingBreakend = breakendList.get(i + 2);
                    final SvCluster followingCluster = followingBreakend.getSV().getCluster();

                    if(followingCluster.getSvCount() == 1 && followingBreakend.getSV().type() == SGL && !followingCluster.isResolved())
                    {
                        long followingProximity = abs(nextBreakend.position() - followingBreakend.position());

                        if (followingProximity < nextProximity)
                            nextBreakend = null;
                    }
                }

                if(nextBreakend == null && prevBreakend == null)
                    continue;

                SvBreakend otherBreakend = null;

                if(nextBreakend != null && prevBreakend != null)
                {
                    otherBreakend = nextProximity < prevProximity ? nextBreakend : prevBreakend;
                }
                else if(nextBreakend != null)
                {
                    otherBreakend = nextBreakend;
                }
                else
                {
                    otherBreakend = prevBreakend;
                }

                SvVarData otherVar = otherBreakend.getSV();
                final SvCluster otherCluster = otherVar.getCluster();

                final String resolvedType = canResolveWithSoloSingle(otherCluster, cluster);

                if(resolvedType != RESOLVED_TYPE_NONE)
                {
                    otherCluster.mergeOtherCluster(cluster);
                    otherCluster.addClusterReason(CLUSTER_REASON_SOLO_SINGLE);
                    otherCluster.setResolved(true, resolvedType);
                    otherCluster.setSynDelDupData(abs(otherVar.position(true) - var.position(true)), 0);

                    clusters.remove(cluster);
                    foundMerges = true;
                    break;
                }
            }
        }

        return foundMerges;
    }

    private final String canResolveWithSoloSingle(SvCluster otherCluster, SvCluster soloSingleCluster)
    {
        // 3 cases:
        // - 2 x SGLs could form a simple DEL or DUP
        // - 2 x SGLs + another SV could form a simple cluster-2 resolved type
        // - a SGL + another SV could form a simple cluster-2 resolved type

        final SvVarData soloSingle = soloSingleCluster.getSV(0);

        if(otherCluster.getSvCount() == 1)
        {
            final SvVarData otherVar = otherCluster.getSV(0);

            if(otherVar.type() == SGL)
            {
                // either both must be NONEs or one be a SGL but without centromeric or telomeric support
                if(otherCluster.isResolved())
                    return RESOLVED_TYPE_NONE;

                String resolvedType = markSinglePairResolvedType(otherVar, soloSingle);

                if(resolvedType == RESOLVED_TYPE_NONE)
                    return RESOLVED_TYPE_NONE;

                LOGGER.debug("cluster({}) SV({}) and cluster({}) SV({}) syntheticType({})",
                        soloSingleCluster.id(), soloSingle.posId(), otherCluster.id(), otherVar.posId(), resolvedType);

                String syntheticType = resolvedType == RESOLVED_TYPE_SGL_PAIR_INS ? "INS" :
                        (resolvedType == RESOLVED_TYPE_SGL_PAIR_DUP ? "DUP" : "DEL");

                soloSingle.addClusterReason(CLUSTER_REASON_SOLO_SINGLE, syntheticType.toString() + "_" + otherVar.id());
                otherVar.addClusterReason(CLUSTER_REASON_SOLO_SINGLE, syntheticType.toString() + "_" + soloSingle.id());

                return resolvedType;
            }
            else
            {
                boolean inconsistentOnStart;
                if(otherVar.hasInconsistentCopyNumberChange(true) && otherVar.chromosome(false).equals(soloSingle.chromosome(true)))
                {
                    inconsistentOnStart = true;
                }
                else if(otherVar.hasInconsistentCopyNumberChange(false) && otherVar.chromosome(true).equals(soloSingle.chromosome(true)))
                {
                    inconsistentOnStart = false;
                }
                else
                {
                    return RESOLVED_TYPE_NONE;
                }

                double cnInconsistency = otherVar.ploidy() - otherVar.copyNumberChange(inconsistentOnStart);

                if(round(cnInconsistency) != round(soloSingle.copyNumberChange(true)))
                    return RESOLVED_TYPE_NONE;

                LOGGER.debug(String.format("cluster(%s) SV(%s) and cluster(%s) SV(%s) potentially resolve CN inconsistency(%.2f vs %.2f)",
                        soloSingleCluster.id(), soloSingle.posId(), otherCluster.id(), otherVar.posId(),
                        cnInconsistency, soloSingle.copyNumberChange(true)));

                soloSingle.addClusterReason(CLUSTER_REASON_SOLO_SINGLE, "CnInc_" + otherVar.id());
                otherVar.addClusterReason(CLUSTER_REASON_SOLO_SINGLE, "CnInc_" + soloSingle.id());

                return RESOLVED_TYPE_SGL_PLUS_INCONSISTENT;
            }
        }
        else
        {
            return RESOLVED_TYPE_NONE;
        }
    }

    public static String markSinglePairResolvedType(final SvVarData sgl1, final SvVarData sgl2)
    {
        if(!sgl1.isNoneSegment() && !sgl2.isNoneSegment())
            return RESOLVED_TYPE_NONE;

        if(sgl1.sglToCentromereOrTelomere() || sgl2.sglToCentromereOrTelomere())
            return RESOLVED_TYPE_NONE;

        final SvBreakend breakend1 = sgl1.getBreakend(true);
        final SvBreakend breakend2 = sgl2.getBreakend(true);

        // to form a simple del or dup, they need to have different orientations
        if(breakend1.orientation() == breakend2.orientation())
            return RESOLVED_TYPE_NONE;

        // check copy number consistency
        double cn1 = sgl2.copyNumberChange(true);
        double cn2 = sgl1.copyNumberChange(true);

        if(!copyNumbersEqual(cn1, cn2))
            return RESOLVED_TYPE_NONE;

        long length = abs(breakend1.position() - breakend2.position());
        int minTiLength = getMinTemplatedInsertionLength(breakend1, breakend2);

        if(length < minTiLength)
        {
            return RESOLVED_TYPE_SGL_PAIR_INS;
        }
        else
        {
            boolean v1First = breakend1.position() < breakend2.position();
            boolean v1PosOrientation = (breakend1.orientation() == 1);

            if(v1First == v1PosOrientation)
                return RESOLVED_TYPE_SGL_PAIR_DEL;
            else
                return RESOLVED_TYPE_SGL_PAIR_DUP;
        }
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

        mDelDupCutoffLength = min(max(mDelDupCutoffLength, MIN_SIMPLE_DUP_DEL_CUTOFF), MAX_SIMPLE_DUP_DEL_CUTOFF);

        LOGGER.debug("sample({}) simple dels and dups: count({}) cutoff-length({}) simpleArms({}) trimCount({})",
                sampleId, lengthsList.size(), mDelDupCutoffLength, simpleArmCount, trimCount);
    }

    public void populateChromosomeBreakendMap(final List<SvVarData> allVariants)
    {
        mChrBreakendMap.clear();

        // add each SV's breakends to a map keyed by chromosome, with the breakends in order of position lowest to highest
        for (final SvVarData var : allVariants)
        {
            addSvToChrBreakendMap(var, mChrBreakendMap);
        }

        // cache indicies for faster look-up
        for (Map.Entry<String, List<SvBreakend>> entry : mChrBreakendMap.entrySet())
        {
            List<SvBreakend> breakendList = entry.getValue();

            for (int i = 0; i < breakendList.size(); ++i)
            {
                final SvBreakend breakend = breakendList.get(i);
                breakend.setChrPosIndex(i);
            }
        }
    }

    public void annotateNearestSvData()
    {
        // mark each SV's nearest other SV and its relationship - neighbouring or overlapping
        // and any duplicate breakends as well
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

                    if(distance <= PERMITED_DUP_BE_DISTANCE && breakend.orientation() == nextBreakend.orientation())
                    {
                        var.setIsDupBreakend(true, breakend.usesStart());
                        nextBreakend.getSV().setIsDupBreakend(true, nextBreakend.usesStart());
                    }
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

    public void markInversionPairTypes(SvCluster cluster)
    {
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

        final SvVarData var1 = cluster.getSV(0);
        final SvVarData var2 = cluster.getSV(1);

        // isSpecificCluster(cluster);

        if(cluster.getLinkedPairs().isEmpty())
        {
            if(var1.orientation(true) == var2.orientation(true))
                return;

            // 2 INVs without any overlap can form 2 deletion bridges
            long syntheticLength;

            if(var1.orientation(true) == 1
            && var1.getDBLink(false) != null && var1.getDBLink(false) == var2.getDBLink(true))
            {
                // take the outermost positions for the DEL length
                syntheticLength = abs(var2.position(false) - var1.position(true));
            }
            else if(var2.orientation(true) == 1
            && var2.getDBLink(false) != null && var2.getDBLink(false) == var1.getDBLink(true))
            {
                syntheticLength = abs(var1.position(false) - var2.position(true));
            }
            else
            {
                LOGGER.debug("cluster({}) inversion-pair no DBs or TIs", cluster.id());
                return;
            }

            cluster.setSynDelDupData(syntheticLength, 0);
            boolean isResolved = syntheticLength < mDelDupCutoffLength;
            cluster.setResolved(isResolved, RESOLVED_TYPE_DEL_INT_TI);
            return;
        }

        resolveSyntheticDelDupCluster(cluster);
    }

    private boolean resolveSyntheticDelDupCluster(SvCluster cluster)
    {
        if(cluster.getLinkedPairs().isEmpty() || cluster.hasReplicatedSVs())
            return false;

        // first work out if there are 1 or 2 templated insertions
        SvLinkedPair linkedPair1 = cluster.getLinkedPairs().get(0);

        SvBreakend lowerBreakend = linkedPair1.getBreakend(true);
        SvBreakend upperBreakend = linkedPair1.getBreakend(false);

        SvVarData lowerSV = lowerBreakend.getSV();
        SvVarData upperSV = upperBreakend.getSV();

        SvBreakend lowerOtherBreakend = lowerBreakend.getOtherBreakend();
        SvBreakend upperOtherBreakend = upperBreakend.getOtherBreakend();

        SvLinkedPair linkedPair2 = cluster.getLinkedPairs().size() == 2 ? cluster.getLinkedPairs().get(1) : null;

        if(linkedPair2 != null && !linkedPair1.hasLinkClash(linkedPair2))
        {
            // existing second link is fine to consider
        }
        else if(lowerOtherBreakend.orientation() != upperOtherBreakend.orientation())
        {
            // convert to either a DB or a TI
            if((lowerOtherBreakend.orientation() == 1 && lowerOtherBreakend.position() <= upperOtherBreakend.position())
            || (upperOtherBreakend.orientation() == 1 && upperOtherBreakend.position() <= lowerOtherBreakend.position()))
            {
                linkedPair2 = new SvLinkedPair(lowerSV, upperSV, LINK_TYPE_DB, lowerOtherBreakend.usesStart(), upperOtherBreakend.usesStart());
            }
            else
            {
                // if this doesn't overlap enough to be a TI, it will be converted into a DB
                linkedPair2 = new SvLinkedPair(lowerSV, upperSV, LINK_TYPE_TI, lowerOtherBreakend.usesStart(), upperOtherBreakend.usesStart());
            }
        }
        else
        {
            LOGGER.debug("cluster({}) ids({} & {}) neither TI nor DB", cluster.id(), lowerSV.id(), upperSV.id());
            return false;
        }

        // set the TI link to be the shorter of the 2
        SvLinkedPair tiPair;
        SvLinkedPair otherPair;

        if(linkedPair1.linkType() == LINK_TYPE_TI && linkedPair2.linkType() == LINK_TYPE_TI)
        {
            tiPair = linkedPair1.length() < linkedPair2.length() ? linkedPair1 : linkedPair2;
            otherPair = linkedPair1 == tiPair ? linkedPair2 : linkedPair1;
        }
        else if(linkedPair1.linkType() == LINK_TYPE_TI)
        {
            // take the DB for the main DEL if one exists
            tiPair = linkedPair1;
            otherPair = linkedPair2;
        }
        else if(linkedPair2.linkType() == LINK_TYPE_TI)
        {
            tiPair = linkedPair2;
            otherPair = linkedPair1;
        }
        else
        {
            // 2 deletion bridges
            tiPair = linkedPair1.length() < linkedPair2.length() ? linkedPair1 : linkedPair2;
            otherPair = linkedPair1 == tiPair ? linkedPair2 : linkedPair1;
        }

        // the other breakends cannot be cross arm
        if(lowerOtherBreakend.arm() != upperOtherBreakend.arm())
            return false;

        if(tiPair.isInferred() && tiTraversesComplexSVs(cluster, tiPair))
            return false;

        long tiPos1 = tiPair.first().position(tiPair.firstLinkOnStart());
        long tiPos2 = tiPair.second().position(tiPair.secondLinkOnStart());

        long otherPos1 = otherPair.first().position(otherPair.firstLinkOnStart());
        long otherPos2 = otherPair.second().position(otherPair.secondLinkOnStart());

        boolean isTIEnclosed = false;

        long lowerOtherPosLimit = otherPos1 < otherPos2 ? otherPos1 : otherPos2;
        lowerOtherPosLimit -= MIN_TEMPLATED_INSERTION_LENGTH;

        long upperOtherPosLimit = otherPos1 > otherPos2 ? otherPos1 : otherPos2;
        upperOtherPosLimit += MIN_TEMPLATED_INSERTION_LENGTH;

        if(tiPos1 >= lowerOtherPosLimit && tiPos1 <= upperOtherPosLimit && tiPos2 >= lowerOtherPosLimit && tiPos2 <= upperOtherPosLimit)
        {
            isTIEnclosed = true;
        }

        String resolvedType = "";

        if(tiPair.linkType() == LINK_TYPE_DB && otherPair != null && otherPair.linkType() == LINK_TYPE_DB)
        {
            resolvedType = RESOLVED_TYPE_DEL_INT_TI;
        }
        else
        {
            if (tiPair.linkType() == LINK_TYPE_DB || otherPair.linkType() == LINK_TYPE_DB)
            {
                if (!isTIEnclosed)
                {
                    resolvedType = RESOLVED_TYPE_DEL_EXT_TI;
                }
                else
                {
                    resolvedType = RESOLVED_TYPE_DEL_INT_TI;
                }
            }
            else
            {
                if (!isTIEnclosed)
                {
                    resolvedType = RESOLVED_TYPE_DUP_EXT_TI;
                }
                else
                {
                    resolvedType = RESOLVED_TYPE_DUP_INT_TI;
                }
            }
        }

        boolean isResolved = otherPair.length() < mDelDupCutoffLength;
        cluster.setResolved(isResolved, resolvedType);

        // correct for DB subtracting 1
        int delDupLength = otherPair.linkType() == LINK_TYPE_DB ? otherPair.length() + 1 : otherPair.length();

        cluster.setSynDelDupData(delDupLength, tiPair.length());

        return true;
    }

    public void markBndPairTypes(SvCluster cluster)
    {
        final SvVarData var1 = cluster.getSV(0);
        final SvVarData var2 = cluster.getSV(1);

        /* possible configurations:
            1. Reciprocal Translocation
            - 2 DBs, no overlappying breakends OR
            - TIs converted to DBs since too short

            2. One set of breakends facing (the TI) the other facing away (the DEL)
            - DEL with TI

            3. Two sets of facing breakends so 2 TIs
            - but rather than a closed loop, one set remain unlinked (the overlap being the DUP)

            Other configurations are nothing
         */

        if(cluster.getLinkedPairs().isEmpty() && arePairedDeletionBridges(var1, var2))
        {
            cluster.setResolved(true, RESOLVED_TYPE_RECIPROCAL_TRANS);
            return;
        }

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
            LOGGER.debug("cluster({} {}) has no linked pairs", cluster.id(), cluster.getDesc());
            return;
        }

        final SvLinkedPair tiLinkedPair;

        if(lp1 != null && lp2 != null)
        {
            // take assembly pair if exists over inferred pair, otherwise take the shorter of the 2
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

        if(var1.orientation(v1OpenOnStart) == var2.orientation(v2OpenOnStart))
            return;

        // the other breakends cannot be cross arm
        if(var1.arm(v1OpenOnStart) != var2.arm(v2OpenOnStart))
            return;

        if(tiTraversesComplexSVs(cluster, tiLinkedPair))
            return;

        String pairDesc = "";

        if(areSectionBreak(var1, var2, v1OpenOnStart, v2OpenOnStart))
        {
            pairDesc = RESOLVED_TYPE_DEL_EXT_TI;
        }
        else
        {
            pairDesc = RESOLVED_TYPE_DUP_EXT_TI;
        }

        long syntheticLength = abs(var1.position(v1OpenOnStart) - var2.position(v2OpenOnStart));

        boolean isResolved = syntheticLength < mDelDupCutoffLength;
        cluster.setResolved(isResolved, pairDesc);
        cluster.setSynDelDupData(syntheticLength, tiLinkedPair.length());
    }

    public boolean markDelDupPairTypes(SvCluster cluster)
    {
        if(cluster.getTypeCount(DUP) == 2)
        {
            // to prevent misclassification of otherwise randomly clustered DUPs,
            // require an assembled TI
            if(cluster.getAssemblyLinkedPairs().size() != 1)
                return false;
        }

        // isSpecificCluster(cluster);

        return resolveSyntheticDelDupCluster(cluster);
    }

    private boolean tiTraversesComplexSVs(final SvCluster cluster, final SvLinkedPair pair)
    {
        // count any non-trivial cluster's SVs crossed by this pair
        final SvBreakend firstBreakend = pair.first().getBreakend(pair.firstLinkOnStart());
        final SvBreakend secondBreakend = pair.second().getBreakend(pair.secondLinkOnStart());
        int lowerIndex = min(firstBreakend.getChrPosIndex(), secondBreakend.getChrPosIndex());
        int upperIndex = max(firstBreakend.getChrPosIndex(), secondBreakend.getChrPosIndex());

        if(lowerIndex >= upperIndex - 1)
            return false;

        final List<SvBreakend> breakendList = mChrBreakendMap.get(firstBreakend.chromosome());

        for (int i = lowerIndex + 1; i <= upperIndex - 1; ++i)
        {
            final SvCluster otherCluster = breakendList.get(i).getSV().getCluster();

            if(otherCluster == cluster)
                continue;

            if (otherCluster.getResolvedType() == RESOLVED_TYPE_SIMPLE_SV
            || otherCluster.getResolvedType() == RESOLVED_TYPE_LINE
            || otherCluster.getResolvedType() == RESOLVED_TYPE_RECIPROCAL_TRANS
            || isFilteredResolvedType(otherCluster.getResolvedType()))
            {
                continue;
            }

            return true;
        }

        return false;
    }

}
