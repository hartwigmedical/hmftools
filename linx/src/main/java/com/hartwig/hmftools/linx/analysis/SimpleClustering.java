package com.hartwig.hmftools.linx.analysis;

import static java.lang.Math.abs;
import static java.lang.Math.max;

import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.isStart;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.NEG_ORIENT;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.POS_ORIENT;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.DEL;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.DUP;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.INS;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.INV;
import static com.hartwig.hmftools.linx.LinxConfig.LNX_LOGGER;
import static com.hartwig.hmftools.linx.analysis.ClusterClassification.getSyntheticLength;
import static com.hartwig.hmftools.linx.analysis.ClusterClassification.isSimpleSingleSV;
import static com.hartwig.hmftools.linx.analysis.ClusteringReason.HIGH_JCN;
import static com.hartwig.hmftools.linx.analysis.ClusteringReason.HOM_LOSS;
import static com.hartwig.hmftools.linx.analysis.ClusteringReason.LOH;
import static com.hartwig.hmftools.linx.analysis.ClusteringReason.LOH_CHAIN;
import static com.hartwig.hmftools.linx.analysis.ClusteringReason.LONG_DEL_DUP_INV;
import static com.hartwig.hmftools.linx.analysis.ClusteringReason.MAJOR_ALLELE_JCN;
import static com.hartwig.hmftools.linx.analysis.ClusteringReason.PROXIMITY;
import static com.hartwig.hmftools.linx.analysis.SvUtilities.MAX_COPY_NUM_DIFF;
import static com.hartwig.hmftools.linx.analysis.SvUtilities.copyNumbersEqual;
import static com.hartwig.hmftools.linx.analysis.SvUtilities.formatJcn;
import static com.hartwig.hmftools.linx.analysis.SvUtilities.getProximity;
import static com.hartwig.hmftools.linx.types.LinxConstants.ADJACENT_JCN_RATIO;
import static com.hartwig.hmftools.linx.types.LinxConstants.HIGH_JCN_THRESHOLD;
import static com.hartwig.hmftools.linx.types.LinxConstants.LOW_JCN_THRESHOLD;
import static com.hartwig.hmftools.linx.types.LinxConstants.MAX_MERGE_DISTANCE;
import static com.hartwig.hmftools.linx.types.ResolvedType.LINE;
import static com.hartwig.hmftools.linx.types.SvVarData.RELATION_TYPE_NEIGHBOUR;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.sv.StructuralVariantType;
import com.hartwig.hmftools.linx.CohortDataWriter;
import com.hartwig.hmftools.linx.CohortFileInterface;
import com.hartwig.hmftools.linx.LinxConfig;
import com.hartwig.hmftools.linx.cn.HomLossEvent;
import com.hartwig.hmftools.linx.cn.LohEvent;
import com.hartwig.hmftools.linx.types.ResolvedType;
import com.hartwig.hmftools.linx.types.SvBreakend;
import com.hartwig.hmftools.linx.types.SvCluster;
import com.hartwig.hmftools.linx.types.SvVarData;

public class SimpleClustering implements CohortFileInterface
{
    private ClusteringState mState;
    private final LinxConfig mConfig;
    private String mSampleId;

    private int mClusteringIndex;
    private final CohortDataWriter mCohortDataWriter;

    public SimpleClustering(final LinxConfig config, final ClusteringState state, final CohortDataWriter cohortDataWriter)
    {
        mState = state;
        mConfig = config;
        mCohortDataWriter = cohortDataWriter;
    }

    private int getNextClusterId() { return mState.getNextClusterId(); }

    public void initialise(final String sampleId)
    {
        mSampleId = sampleId;
        mClusteringIndex = 0;
    }

    public void clusterByProximity(final List<SvCluster> clusters)
    {
        int proximityDistance = mConfig.ProximityDistance;

        // walk through each chromosome and breakend list
        for (final Map.Entry<String, List<SvBreakend>> entry : mState.getChrBreakendMap().entrySet())
        {
            final List<SvBreakend> breakendList = entry.getValue();

            int currentIndex = 0;
            while (currentIndex < breakendList.size())
            {
                final SvBreakend breakend = breakendList.get(currentIndex);
                SvVarData var = breakend.getSV();

                SvBreakend nextBreakend = null;
                int nextIndex = currentIndex + 1;

                for (; nextIndex < breakendList.size(); ++nextIndex)
                {
                    final SvBreakend nextBe = breakendList.get(nextIndex);
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
                else if (abs(nextBreakend.position() - breakend.position()) > proximityDistance)
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
                            cluster.addClusterReason(PROXIMITY);
                            clusters.add(cluster);
                        }
                        else if (cluster != null && nextCluster != null)
                        {
                            // keep one and remove the other
                            cluster.mergeOtherCluster(nextCluster, false);
                            cluster.addClusterReason(PROXIMITY);
                            clusters.remove(nextCluster);
                        }
                        else
                        {
                            if (cluster == null)
                            {
                                nextCluster.addVariant(var);
                                nextCluster.addClusterReason(PROXIMITY);
                            }
                            else
                            {
                                cluster.addVariant(nextVar);
                                cluster.addClusterReason(PROXIMITY);
                            }
                        }

                        if (!var.hasClusterReason(PROXIMITY))
                            var.addClusterReason(PROXIMITY, nextVar.id());

                        if (!nextVar.hasClusterReason(PROXIMITY))
                            nextVar.addClusterReason(PROXIMITY, var.id());

                        // checkClusteringClonalDiscrepancy(var, nextVar, PROXIMITY);
                    }
                }

                // move the index to the SV which was just proximity cluster so the next comparison is with the closest candidate
                currentIndex = nextIndex;
            }
        }
    }

    public void addClusterReasons(final SvVarData var1, final SvVarData var2, final ClusteringReason clusterReason)
    {
        var1.addClusterReason(clusterReason, var2.id());
        var2.addClusterReason(clusterReason, var1.id());

        if(mConfig.hasMultipleSamples())
        {
            logClusteringDetails(var1, var2, clusterReason);
        }
    }

    private static final String COHORT_WRITER_CLUSTER_HISTORY = "ClusterHistory";

    @Override
    public String fileType() { return COHORT_WRITER_CLUSTER_HISTORY; }

    @Override
    public BufferedWriter createWriter(final String outputDir)
    {
        try
        {
            String outputFileName = outputDir + "LNX_CLUSTERING_HISTORY.csv";

            BufferedWriter writer = createBufferedWriter(outputFileName, false);

            // definitional fields
            writer.write("SampleId,MergeIndex,ClusterId1,SvId1,ClusterCount1,ClusterId2,SvId2,ClusterCount2");
            writer.write(",Reason,MinDistance");
            writer.newLine();
            return writer;
        }
        catch(IOException e)
        {
            LNX_LOGGER.error("failed to open and write output file headers");
            return null;
        }
    }

    protected void logClusteringDetails(final SvVarData var1, final SvVarData var2, final ClusteringReason reason)
    {
        if(!mConfig.Output.WriteClusterHistory)
            return;

        StringBuilder sb = new StringBuilder();
        sb.append(String.format("%s,%d", mSampleId, mClusteringIndex));

        sb.append(String.format(",%d,%d,%d,%d,%d,%d",
                var1.getCluster().id(), var1.id(), var1.getCluster().getSvCount(),
                var2.getCluster().id(), var2.id(), var2.getCluster().getSvCount()));

        int breakendDistance = getProximity(var1, var2);
        sb.append(String.format(",%s,%d", reason, breakendDistance));

        mCohortDataWriter.write(this, Lists.newArrayList(sb.toString()));

        ++mClusteringIndex;
    }

    public void mergeClusters(List<SvCluster> clusters)
    {
        // first apply replication rules since this can affect consistency
        for(SvCluster cluster : clusters)
        {
            if(cluster.isResolved())
                continue;

            markClusterLongDelDups(cluster);
            markClusterInversions(cluster);
        }

        int initClusterCount = clusters.size();

        mergeOnLOHEvents(clusters);

        mergeOnHighFacingJcn(clusters);

        mergeOnMajorAlleleJcnBounds(clusters);

        mergeLOHResolvingClusters(clusters);

        // the merge must be run a few times since as clusters grow, more single SVs and other clusters
        // will then fall within the bounds of the new larger clusters
        boolean foundMerges = true;
        int iterations = 0;

        while(foundMerges)
        {
            foundMerges = false;

            if(mergeOnOverlappingInvDupDels(clusters, false))
                foundMerges = true;

            ++iterations;

            if(iterations >= 10)
            {
                if(foundMerges)
                    LNX_LOGGER.warn("sample({}) exiting simple merge loop after {} iterations with merge just found", mSampleId, iterations);
                break;
            }
        }

        if(clusters.size() < initClusterCount)
        {
            LNX_LOGGER.debug("reduced cluster count({} -> {}) iterations({})", initClusterCount, clusters.size(), iterations);
        }
    }

    public static boolean hasLowJcn(final SvVarData var)
    {
        if(var.type() == INS)
            return false;

        return var.jcnMax() < LOW_JCN_THRESHOLD;
    }

    private void mergeOnLOHEvents(List<SvCluster> clusters)
    {
        if (mState.getLohEventList().isEmpty() && mState.getHomLossList().isEmpty())
            return;

        // first link up breakends joined by an LOH with no multi-SV hom-loss events within
        for (final LohEvent lohEvent : mState.getLohEventList())
        {
            if (!lohEvent.matchedBothSVs() || lohEvent.hasIncompleteHomLossEvents())
                continue; // cannot be used for clustering

            SvBreakend lohBeStart = lohEvent.getBreakend(true);
            SvCluster lohClusterStart = lohBeStart.getCluster();
            SvVarData lohSvStart = lohBeStart.getSV();
            SvBreakend lohBeEnd = lohEvent.getBreakend(false);
            SvCluster lohClusterEnd = lohBeEnd.getCluster();
            SvVarData lohSvEnd = lohBeEnd.getSV();

            if (lohClusterStart != lohClusterEnd)
            {
                if (lohClusterStart.hasLinkingLineElements() || lohClusterEnd.hasLinkingLineElements())
                    continue;

                if(variantsHaveDifferentJcn(lohBeStart, lohBeEnd))
                    continue;

                LNX_LOGGER.debug("cluster({} svs={}) merges in other cluster({} svs={}) on LOH event({})",
                        lohClusterStart.id(), lohClusterStart.getSvCount(), lohClusterEnd.id(), lohClusterEnd.getSvCount(), lohEvent);

                addClusterReasons(lohSvStart, lohSvEnd, LOH);
                lohClusterStart.addClusterReason(LOH);

                lohClusterStart.mergeOtherCluster(lohClusterEnd);
                clusters.remove(lohClusterEnd);
            }
        }

        // search for LOH events which are clustered but which contain hom-loss events which aren't clustered
        // (other than simple DELs) which already will have been handled
        for (final LohEvent lohEvent : mState.getLohEventList())
        {
            if (!lohEvent.hasIncompleteHomLossEvents())
                continue;

            // already clustered - search for hom-loss events contained within this LOH that isn't clustered
            if(lohEvent.clustered() || lohEvent.wholeArmLoss())
            {
                List<HomLossEvent> unclusteredHomLossEvents = lohEvent.getHomLossEvents().stream()
                        .filter(HomLossEvent::matchedBothSVs)
                        .filter(x -> !x.sameSV())
                        .filter(x -> x.PosStart > lohEvent.PosStart)
                        .filter(x -> x.PosEnd < lohEvent.PosEnd)
                        .filter(x -> !x.clustered())
                        .collect(Collectors.toList());

                for (final HomLossEvent homLoss : unclusteredHomLossEvents)
                {
                    SvBreakend homLossBeStart = homLoss.getBreakend(true);
                    SvBreakend homLossBeEnd = homLoss.getBreakend(false);

                    SvCluster cluster = homLossBeStart.getCluster();
                    SvCluster otherCluster = homLossBeEnd.getCluster();

                    if(cluster == otherCluster) // protect against clusters already merged or removed
                        continue;

                    if(variantsHaveDifferentJcn(homLossBeStart, homLossBeEnd))
                        continue;

                    LNX_LOGGER.debug("cluster({} svs={}) merges in other cluster({} svs={}) on hom-loss({}: {} -> {}) inside LOH event({} -> {})",
                            cluster.id(), cluster.getSvCount(), otherCluster.id(), otherCluster.getSvCount(),
                            homLoss.Chromosome, homLoss.PosStart, homLoss.PosEnd, lohEvent.PosStart, lohEvent.PosEnd);

                    addClusterReasons(homLossBeStart.getSV(), homLossBeEnd.getSV(), HOM_LOSS);
                    cluster.addClusterReason(HOM_LOSS);

                    cluster.mergeOtherCluster(otherCluster);
                    clusters.remove(otherCluster);
               }

                continue;
            }

            if (!lohEvent.matchedBothSVs() || lohEvent.clustered())
                continue;

            if (lohEvent.getBreakend(true).getCluster().hasLinkingLineElements()
            || lohEvent.getBreakend(false).getCluster().hasLinkingLineElements())
            {
                continue;
            }

            // now look for an LOH with unclustered breakends but which contains only clustered hom-loss events
            boolean hasIncompleteHomLossEvent = false;

            for(final HomLossEvent homLossEvent : lohEvent.getHomLossEvents())
            {
                if(!(homLossEvent.PosStart > lohEvent.PosStart && homLossEvent.PosEnd < lohEvent.PosEnd))
                {
                    // handle overlapping separately
                    hasIncompleteHomLossEvent = true;
                    break;
                }

                if(!homLossEvent.matchedBothSVs())
                {
                    hasIncompleteHomLossEvent = true;
                    break;
                }

                if(!homLossEvent.clustered())
                {
                    hasIncompleteHomLossEvent = true;
                    break;
                }
            }

            if(!hasIncompleteHomLossEvent)
            {
                // all hom-loss events involving more than 1 SV were clustered, so clustered the LOH SVs
                SvBreakend lohBeStart = lohEvent.getBreakend(true);
                SvBreakend lohBeEnd = lohEvent.getBreakend(false);

                if(variantsHaveDifferentJcn(lohBeStart, lohBeEnd))
                    continue;

                SvCluster cluster = lohBeStart.getCluster();
                SvCluster otherCluster = lohBeEnd.getCluster();

                lohEvent.setIsValid(true);

                LNX_LOGGER.debug("cluster({} svs={}) merges in other cluster({} svs={}) on no unclustered hom-loss events",
                        cluster.id(), cluster.getSvCount(), otherCluster.id(), otherCluster.getSvCount());

                addClusterReasons(lohBeStart.getSV(), lohBeEnd.getSV(), HOM_LOSS);
                cluster.addClusterReason(HOM_LOSS);

                cluster.mergeOtherCluster(otherCluster);
                clusters.remove(otherCluster);
            }

            // finally look for overlapping LOH and hom-loss events where all but 2 of the breakends are clustered
            // resulting in a clustering of a breakend from the LOH with one of the hom-loss breakends

            List<SvBreakend> unclusteredBreakends = Lists.newArrayList();
            boolean allBreakendsValid = true;

            for(int se = SE_START; se <= SE_END; ++se)
            {
                boolean isStart = isStart(se);
                unclusteredBreakends.add(lohEvent.getBreakend(isStart));

                for(final HomLossEvent homLossEvent : lohEvent.getHomLossEvents())
                {
                    if(!homLossEvent.matchedBothSVs())
                    {
                        allBreakendsValid = false;
                        break;
                    }

                    unclusteredBreakends.add(homLossEvent.getBreakend(isStart));
                }

                if(!allBreakendsValid)
                    break;
            }

            if(allBreakendsValid)
            {
                int i = 0;
                while(i < unclusteredBreakends.size())
                {
                    final SvBreakend be1 = unclusteredBreakends.get(i);

                    boolean found = false;
                    for(int j = i+1; j < unclusteredBreakends.size(); ++j)
                    {
                        final SvBreakend be2 = unclusteredBreakends.get(j);

                        if(be1.getCluster() == be2.getCluster())
                        {
                            unclusteredBreakends.remove(j);
                            unclusteredBreakends.remove(i);
                            found = true;
                            break;
                        }
                    }

                    if(!found)
                        ++i;
                }
            }

            if(unclusteredBreakends.size() == 2)
            {
                SvBreakend breakend1 = unclusteredBreakends.get(0);
                SvBreakend breakend2 = unclusteredBreakends.get(1);
                SvCluster cluster = breakend1.getCluster();
                SvCluster otherCluster = breakend2.getCluster();

                lohEvent.setIsValid(true);

                if(cluster == otherCluster)
                    continue;

                if(variantsHaveDifferentJcn(breakend1, breakend2))
                    continue;

                LNX_LOGGER.debug("cluster({} svs={}) merges in other cluster({} svs={}) on unclustered LOH and hom-loss breakends",
                        cluster.id(), cluster.getSvCount(), otherCluster.id(), otherCluster.getSvCount());

                addClusterReasons(breakend1.getSV(), breakend2.getSV(), HOM_LOSS);
                cluster.addClusterReason(HOM_LOSS);

                cluster.mergeOtherCluster(otherCluster);
                clusters.remove(otherCluster);
            }
        }
    }

    private void markClusterInversions(final SvCluster cluster)
    {
        if(cluster.getTypeCount(INV) == 0)
            return;

        // skip cluster-2s which resolved to a simple type and not long
        if(cluster.isResolved() && cluster.getResolvedType().isSimple())
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
        // find and record any long DEL or DUP for merging, including long synthetic ones
        if(cluster.isSyntheticType() && cluster.getResolvedType().isSimple() && !cluster.isResolved())
        {
            int syntheticLength = getSyntheticLength(cluster);

            if ((cluster.getResolvedType() == ResolvedType.DEL && syntheticLength >= mState.getDelCutoffLength())
            || (cluster.getResolvedType() == ResolvedType.DUP && syntheticLength >= mState.getDupCutoffLength()))
            {
                for (final SvVarData var : cluster.getSVs())
                {
                    cluster.registerLongDelDup(var);
                }
            }
        }

        if(cluster.getTypeCount(DEL) > 0 || cluster.getTypeCount(DUP) > 0)
        {
            for(final SvVarData var : cluster.getSVs())
            {
                if(var.isCrossArm())
                    continue;

                if(exceedsDupDelCutoffLength(var.type(), var.length()))
                {
                    cluster.registerLongDelDup(var);
                }
            }
        }
    }

    public void mergeLongDelDupClusters(List<SvCluster> clusters)
    {
        boolean foundMerges = true;
        int iterations = 0;

        while(foundMerges)
        {
            foundMerges = false;

            if(mergeOnOverlappingInvDupDels(clusters, true))
                foundMerges = true;

            ++iterations;

            if(iterations >= 5)
            {
                if(foundMerges)
                    LNX_LOGGER.warn("sample({}) exiting simple merge loop after {} iterations with merge just found", mSampleId, iterations);
                break;
            }
        }
    }

    private boolean mergeOnOverlappingInvDupDels(List<SvCluster> clusters, boolean allowDelDupOverlaps)
    {
        // merge any clusters with overlapping inversions, long dels or long dups on the same arm
        List<SvCluster> longDDIClusters = clusters.stream()
                .filter(x -> !x.getInversions().isEmpty() || !x.getLongDelDups().isEmpty())
                .filter(x -> !x.hasLinkingLineElements())
                .collect(Collectors.toList());

        if(longDDIClusters.size() <= 1)
            return false;

        LNX_LOGGER.debug("checking long {}} overlaps for {} clusters",
                !allowDelDupOverlaps ? "DEL_DUP-requiring-INV" : "multiple DDI overlaps", longDDIClusters.size());

        List<SvCluster> mergedClusters = Lists.newArrayList();

        int index1 = 0;
        while(index1 < longDDIClusters.size())
        {
            SvCluster cluster1 = longDDIClusters.get(index1);

            if(mergedClusters.contains(cluster1))
            {
                ++index1;
                continue;
            }

            boolean mergedOtherClusters = false;
            ClusteringReason mergeReason = ClusteringReason.NONE;

            List<SvVarData> cluster1Svs = Lists.newArrayList(cluster1.getLongDelDups());
            cluster1Svs.addAll(cluster1.getInversions());

            int index2 = index1 + 1;
            while(index2 < longDDIClusters.size())
            {
                SvCluster cluster2 = longDDIClusters.get(index2);

                if(mergedClusters.contains(cluster2))
                {
                    ++index2;
                    continue;
                }

                List<SvVarData> cluster2Svs = Lists.newArrayList(cluster2.getLongDelDups());
                cluster2Svs.addAll(cluster2.getInversions());

                int delDupOverlapCount = 0;
                int closeLinkPairs = 0;

                for (final SvVarData var1 : cluster1Svs)
                {
                    for (final SvVarData var2 : cluster2Svs)
                    {
                        boolean pairContainsInv = var1.type() == INV || var2.type() == INV;

                        if(!allowDelDupOverlaps && !pairContainsInv)
                            continue;

                        if(!var1.chromosome(true).equals(var2.chromosome(true)))
                            continue;

                        if(var1.position(false) < var2.position(true) || var1.position(true) > var2.position(false))
                            continue;

                        boolean enclosed = (var1.position(true) < var2.position(true) && var1.position(false) > var2.position(false))
                                || (var2.position(true) < var1.position(true) && var2.position(false) > var1.position(false));

                        if(allowDelDupOverlaps && enclosed)
                            continue;

                        // check for conflicting LOH / hom-loss events
                        if(variantsViolateLohHomLoss(var1, var2))
                        {
                            LNX_LOGGER.trace("cluster({}) SV({}) and cluster({}) var({}) have conflicting LOH & hom-loss events",
                                    cluster1.id(), var1.id(), cluster2.id(), var2.id());
                            continue;
                        }

                        if(variantsHaveDifferentJcn(var1, var2))
                            continue;

                        if(allowDelDupOverlaps && (!(copyNumbersEqual(var1.copyNumber(true), var2.copyNumber(true)))
                                || !copyNumbersEqual(var1.copyNumber(false), var2.copyNumber(false))))
                        {
                            continue;
                        }

                        if(!pairContainsInv)
                            ++delDupOverlapCount;

                        boolean[] closeLinkData = breakendsInCloseLink(var1, var2);

                        if(closeLinkData[CLOSE_BREAKS_TI_DB])
                            ++closeLinkPairs;

                        boolean closeLink = closeLinkData[CLOSE_BREAKS_PROXIMATE];

                        // either require an INV to be a part of the overlap, or at least 3 DELs or DUPS
                        // and either 1 closer TI or DB link or at least 3 outside the range

                        if((closeLink && pairContainsInv) || (delDupOverlapCount >= 3 && closeLinkPairs >= 3))
                        {
                            if(closeLink && pairContainsInv)
                            {
                                LNX_LOGGER.debug("cluster({}) SV({} {}) and cluster({}) SV({} {}) have INV-DEL-DUP overlap",
                                        cluster1.id(), var1.posId(), var1.type(), cluster2.id(), var2.posId(), var2.type());
                            }
                            else
                            {
                                LNX_LOGGER.debug("cluster({}) and cluster({}) have {} DEL-DUP overlaps and {} close pairs",
                                        cluster1.id(), cluster2.id(), delDupOverlapCount, closeLinkPairs);
                            }

                            mergeReason = LONG_DEL_DUP_INV;
                            addClusterReasons(var1, var2, mergeReason);
                            mergedOtherClusters = true;
                            break;
                        }
                    }

                    if(mergedOtherClusters)
                        break;
                }

                if(mergedOtherClusters)
                {
                    cluster1.mergeOtherCluster(cluster2);
                    cluster1.addClusterReason(mergeReason);
                    mergedClusters.add(cluster2);
                    break;
                }
                else
                {
                    ++index2;
                }
            }

            if(mergedOtherClusters)
                continue; // repeat this cluster after merging in another's SVs
            else
                ++index1;
        }

        if(mergedClusters.isEmpty())
            return false;

        mergedClusters.forEach(x -> clusters.remove(x));
        return true;
    }

    protected static boolean variantsViolateLohHomLoss(final SvVarData var1, final SvVarData var2)
    {
        for(int se1 = SE_START; se1 <= SE_END; ++se1)
        {
            if(se1 == SE_END && var1.isSglBreakend())
                continue;

            for(int se2 = SE_START; se2 <= SE_END; ++se2)
            {
                if(se2 == SE_END && var2.isSglBreakend())
                    continue;

                if(breakendsViolateLohHomLoss(var1.getBreakend(se1), var2.getBreakend(se2)))
                    return true;

                if(breakendsViolateLohHomLoss(var2.getBreakend(se2), var1.getBreakend(se1)))
                    return true;
            }
        }

        return false;
    }

    protected static boolean breakendsViolateLohHomLoss(final SvBreakend breakend, final SvBreakend otherBreakend)
    {
        // cannot merge to clusters if the reason for merging is 2 of their breakends forming an LOH and the other inside its bounds
        if(!breakend.chromosome().equals(otherBreakend.chromosome()))
            return false;

        List<LohEvent> lohEvents = breakend.getCluster().getLohEvents().stream()
                .filter(x -> x.getBreakend(true) == breakend || x.getBreakend(false) == breakend)
                .collect(Collectors.toList());

        if(lohEvents.isEmpty())
            return false;

        return lohEvents.stream().anyMatch(x -> otherBreakend.position() > x.PosStart && otherBreakend.position() < x.PosEnd);
    }

    protected static boolean variantsHaveDifferentJcn(final SvBreakend breakend1, final SvBreakend breakend2)
    {
        return variantsHaveDifferentJcn(breakend1.getSV(), breakend2.getSV());
    }

    protected static boolean variantsHaveDifferentJcn(final SvVarData var1, final SvVarData var2)
    {
        // We identify high confidence subclonal variants using the uncertainty bounds, using the threshold of maximum ploidy < 0.75.
        // Clonal and subclonal variants are unlikely to have occurred at the same time and hence all subclonal variants are excluded from
        // clustering with any variant that does not overlap in ploidy uncertainty and does not have a ploidy within 0.5 of the subclonal
        // variants. Proximity clustering is still allowed, since the ploidy estimates for proximate variants are more uncertain.
        boolean var1IsLowJcn = var1.jcn() < LOW_JCN_THRESHOLD;
        boolean var2IsLowJcn = var2.jcn() < LOW_JCN_THRESHOLD;

        if(var1IsLowJcn == var2IsLowJcn)
            return false;

        double ploidyDiff = abs(var1.jcn() - var2.jcn());

        if(ploidyDiff < MAX_COPY_NUM_DIFF)
            return false;

        if(var1IsLowJcn && var1.jcnMax() < var2.jcnMin())
            return true;
        else if(var2IsLowJcn && var2.jcnMax() < var1.jcnMin())
            return true;
        else
            return false;
    }

    protected boolean exceedsDupDelCutoffLength(StructuralVariantType type, int length)
    {
        if(type == DEL)
            return length > mState.getDelCutoffLength();
        else if(type == DUP)
            return length > mState.getDelCutoffLength();
        else
            return false;
    }

    private static int CLOSE_BREAKS_TI_DB = 0;
    private static int CLOSE_BREAKS_PROXIMATE = 1;

    private static boolean[] breakendsInCloseLink(final SvVarData var1, final SvVarData var2)
    {
        // test whether the breakends form a DB or TI within the long distance threshold
        boolean[] data = {false, false};

        for(int se1 = SE_START; se1 <= SE_END; ++se1)
        {
            if(se1 == SE_END && var1.isSglBreakend())
                continue;

            final SvBreakend breakend1 = var1.getBreakend(se1);

            for(int se2 = SE_START; se2 <= SE_END; ++se2)
            {
                final SvBreakend breakend2 = var2.getBreakend(se2);

                if(!breakend1.getChrArm().equals(breakend2.getChrArm()))
                    continue;

                if(breakend1.orientation() != breakend2.orientation())
                {
                    data[CLOSE_BREAKS_TI_DB] = true;

                    if(abs(breakend1.position() - breakend2.position()) < MAX_MERGE_DISTANCE)
                    {
                        data[CLOSE_BREAKS_PROXIMATE] = true;
                        return data;
                    }
                }
           }
        }

        return data;
    }

    protected static boolean skipClusterType(final SvCluster cluster)
    {
        if(cluster.getResolvedType() == LINE)
            return true;

        if(isSimpleSingleSV(cluster) && cluster.getSV(0).getNearestSvRelation() == RELATION_TYPE_NEIGHBOUR)
            return true;

        if(cluster.isResolved())
            return true;

        return false;
    }

    private boolean mergeOnMajorAlleleJcnBounds(List<SvCluster> clusters)
    {
        /* The major allele ploidy of a segment is the maximum ploidy any derivative chromosome which includes that segment can have.
        Hence a breakend cannot chain completely across a region with major allele ploidy < ploidy of the breakend, or
        partially across the region with a chain of ploidy more than the major allele.

    	Therefore cluster any breakend with the next 1 or more facing breakends (excluding LINE & assembled & simple non overlapping DEL/DUP)
    	if the major allele ploidy in the segment immediately after the facing breakend is lower than the breakend ploidy - sum(facing breakend ploidies).
    	We limit this clustering to a proximity of 5 million bases and bounded by the centromere, since although more distal events are on
    	the same chromosome may be definitely on the same derivative chromosome, this does necessarily imply they occurred concurrently.
        */

        // method:
        // walk through each cluster's chromosome breakend lists in turn in both directions
        // from each breakend walk forward and subtract any facing breakend's ploidy in the same cluster
        // if an opposing unclustered breakend is encountered and the major AP in the segment after the unclustered breakend is less than
        // the clustered net breakend JCN, then merge in the unclustered breakend, subtract its ploidy and continue

        List<SvCluster> mergedClusters = Lists.newArrayList();

        int clusterIndex = 0;
        while(clusterIndex < clusters.size())
        {
            SvCluster cluster = clusters.get(clusterIndex);

            if(mergedClusters.contains(cluster) || skipClusterType(cluster))
            {
                ++clusterIndex;
                continue;
            }

            boolean mergedOtherClusters = false;

            for (final Map.Entry<String, List<SvBreakend>> entry : cluster.getChrBreakendMap().entrySet())
            {
                List<SvBreakend> breakendList = entry.getValue();

                List<SvBreakend> fullBreakendList = mState.getChrBreakendMap().get(entry.getKey());

                // walk through this list from each direction
                for(int i = 0; i <= 1; ++i)
                {
                    boolean traverseUp = (i == 0);
                    int index = traverseUp ? -1 : breakendList.size(); // will be incremented on first pass

                    while(true)
                    {
                        index += traverseUp ? 1 : -1;

                        if(index < 0 || index >= breakendList.size())
                            break;

                        SvBreakend breakend = breakendList.get(index);

                        if((breakend.orientation() == -1) != traverseUp)
                            continue;

                        double breakendJcn = breakend.jcn();

                        // now walk from this location onwards using the full arm breakend list
                        int chrIndex = breakend.getChrPosIndex();

                        List<SvBreakend> opposingBreakends = Lists.newArrayList();

                        while(true)
                        {
                            chrIndex += traverseUp ? 1 : -1;

                            if(chrIndex < 0 || chrIndex >= fullBreakendList.size())
                                break;

                            SvBreakend nextBreakend = fullBreakendList.get(chrIndex);

                            if(abs(nextBreakend.position() - breakend.position()) > MAX_MERGE_DISTANCE)
                                break;

                            if(nextBreakend.arm() != breakend.arm())
                                break;

                            if(nextBreakend.orientation() == breakend.orientation())
                                continue;

                            if(skipClusterType(nextBreakend.getCluster()))
                                continue;

                            if(nextBreakend.getCluster() == cluster)
                            {
                                breakendJcn -= nextBreakend.jcn();

                                if(breakendJcn <= 0)
                                    break;

                                continue;
                            }

                            if(nextBreakend.isAssembledLink())
                                continue;

                            if(variantsHaveDifferentJcn(breakend, nextBreakend))
                                continue;

                            opposingBreakends.add(nextBreakend);

                            // should this next breakend be merged in?
                            double followingMajorAP = nextBreakend.majorAlleleJcn(!traverseUp);

                            if(!copyNumbersEqual(breakendJcn, followingMajorAP) && breakendJcn > followingMajorAP)
                            {
                                // take the highest of the opposing breakends which were encountered, and if the highest match, then the first
                                double maxOpposingJcn = opposingBreakends.stream().mapToDouble(x -> x.jcn()).max().getAsDouble();

                                SvBreakend opposingBreakend = opposingBreakends.stream()
                                        .filter(x -> copyNumbersEqual(x.jcn(), maxOpposingJcn)).findFirst().get();

                                SvCluster otherCluster = opposingBreakend.getCluster();

                                LNX_LOGGER.debug("cluster({}) breakend({} netJCN={}) merges cluster({}) breakend({} ploidy={}) prior to MAP drop({})",
                                        cluster.id(), breakend, formatJcn(breakendJcn), otherCluster.id(), opposingBreakend.toString(),
                                        formatJcn(opposingBreakend.jcn()), formatJcn(followingMajorAP));

                                addClusterReasons(breakend.getSV(), opposingBreakend.getSV(), MAJOR_ALLELE_JCN);
                                otherCluster.addClusterReason(MAJOR_ALLELE_JCN);
                                cluster.addClusterReason(MAJOR_ALLELE_JCN);

                                cluster.mergeOtherCluster(otherCluster);

                                mergedClusters.add(otherCluster);

                                mergedOtherClusters = true;
                                break;
                            }
                        }

                        if(mergedOtherClusters)
                            break;
                    }

                    if(mergedOtherClusters)
                        break;
                }

                if(mergedOtherClusters)
                    break;
            }

            if(mergedOtherClusters)
            {
                // repeat this cluster
            }
            else
            {
                ++clusterIndex;
            }
        }

        if(mergedClusters.isEmpty())
            return false;

        mergedClusters.forEach(x -> clusters.remove(x));

        return true;
    }

    private boolean mergeOnHighFacingJcn(List<SvCluster> clusters)
    {
        // merge any facing breakends whose JCNs exceed the threshold, regardless of distance, as long as the region in between
        // has continuous major allele copy number at or above this same threshold
        final List<SvCluster> mergedClusters = Lists.newArrayList();

        for (Map.Entry<String, List<SvBreakend>> entry : mState.getChrBreakendMap().entrySet())
        {
            final List<SvBreakend> breakendList = entry.getValue();

            for(final SvBreakend breakend : breakendList)
            {
                if(breakend.orientation() == POS_ORIENT)
                    continue;

                if(breakend.getCluster().getResolvedType() == LINE)
                    continue;

                if(breakend.jcn() < HIGH_JCN_THRESHOLD)
                    continue;

                double adjacentMaJcn = breakend.majorAlleleJcn(true);

                if(breakend.jcn() / max(adjacentMaJcn, 0.01) < ADJACENT_JCN_RATIO)
                    continue;

                // now look for from here through a region of sustained high major allele CN until a facing breakend is reached
                for(int index = breakend.getChrPosIndex() + 1; index < breakendList.size(); ++index)
                {
                    final SvBreakend nextBreakend = breakendList.get(index);

                    if(nextBreakend.orientation() == NEG_ORIENT)
                        continue;

                    double nextAdjacentMaJcn = nextBreakend.majorAlleleJcn(false);

                    boolean isHighFacingBreakend = nextBreakend.jcn() >= HIGH_JCN_THRESHOLD &&
                            nextBreakend.jcn() / max(nextAdjacentMaJcn, 0.01) >= ADJACENT_JCN_RATIO;

                    if(isHighFacingBreakend)
                    {
                        if(nextBreakend.getCluster() == breakend.getCluster() || nextBreakend.getCluster().getResolvedType() == LINE)
                            continue;

                        // found a cluster to merge
                        SvCluster cluster = breakend.getCluster();
                        SvCluster otherCluster = nextBreakend.getCluster();

                        LNX_LOGGER.debug("cluster({}) breakend({} netJCN={}) merges cluster({}) breakend({} ploidy={}) high facing JCN",
                                cluster.id(), breakend, formatJcn(breakend.jcn()), otherCluster.id(), nextBreakend.toString(),
                                formatJcn(nextBreakend.jcn()));

                        addClusterReasons(breakend.getSV(), nextBreakend.getSV(), HIGH_JCN);
                        otherCluster.addClusterReason(HIGH_JCN);
                        cluster.addClusterReason(HIGH_JCN);

                        cluster.mergeOtherCluster(otherCluster);
                        mergedClusters.add(otherCluster);
                    }

                    if(nextBreakend.majorAlleleJcn(false) < HIGH_JCN_THRESHOLD)
                        break;

                    /*
                    // find where the region drops below the required threshold
                    if(nextBreakend.majorAlleleJcn(false) < HIGH_JCN_THRESHOLD)
                    {
                        // check if it's due to another high JCN breakend
                        if(nextBreakend.jcn() < HIGH_JCN_THRESHOLD)
                            break;

                        double nextAdjacentMaJcn = nextBreakend.majorAlleleJcn(false);

                        if(nextBreakend.jcn() / max(nextAdjacentMaJcn, 0.01) < ADJACENT_JCN_RATIO)
                            break;

                        if(nextBreakend.getCluster() == breakend.getCluster() || nextBreakend.getCluster().getResolvedType() == LINE)
                            break;

                        // found a cluster to merge
                        SvCluster cluster = breakend.getCluster();
                        SvCluster otherCluster = nextBreakend.getCluster();

                        LNX_LOGGER.debug("cluster({}) breakend({} netJCN={}) merges cluster({}) breakend({} ploidy={}) high facing JCN",
                                cluster.id(), breakend, formatJcn(breakend.jcn()), otherCluster.id(), nextBreakend.toString(),
                                formatJcn(nextBreakend.jcn()));

                        addClusterReasons(breakend.getSV(), nextBreakend.getSV(), HIGH_JCN);
                        otherCluster.addClusterReason(HIGH_JCN);
                        cluster.addClusterReason(HIGH_JCN);

                        cluster.mergeOtherCluster(otherCluster);
                        mergedClusters.add(otherCluster);
                    }
                    */
                }
            }
        }

        if(mergedClusters.isEmpty())
            return false;

        mergedClusters.forEach(x -> clusters.remove(x));

        return true;
    }

    private boolean mergeLOHResolvingClusters(List<SvCluster> clusters)
    {
        // No breakend in a cluster can chain across an LOH which has been caused by a breakend in the same cluster.
        // Hence if the other breakend of a DUP type variant bounding an LOH can only chain to only one available (not assembled, not LINE)
        // breakend prior to the LOH, then we cluster the DUP and the other breakend.

        List<SvCluster> clustersWithLohEvents = clusters.stream()
                .filter(x -> !x.getLohEvents().isEmpty())
                .filter(x -> !x.hasLinkingLineElements())
                .collect(Collectors.toList());

        List<SvCluster> mergedClusters = Lists.newArrayList();

        for(SvCluster lohCluster : clustersWithLohEvents)
        {
            if(mergedClusters.contains(lohCluster)) // if this has LOH events they will have been added to the parent cluster
                continue;

            List<LohEvent> lohEvents = lohCluster.getLohEvents();

            int lohIndex = 0; // used since merging another cluster can add more LOH events
            while(lohIndex < lohEvents.size())
            {
                LohEvent lohEvent = lohEvents.get(lohIndex);

                if(!lohEvent.isValid())
                {
                    ++lohIndex;
                    continue;
                }

                for(int be = SE_START; be <= SE_END; ++be)
                {
                    SvBreakend lohBreakend = lohEvent.getBreakend(isStart(be));

                    if(lohBreakend == null || lohBreakend.type() != DUP)
                        continue;

                    // it's possible that the breakends for this LOH are not clustered, eg if one is LINE
                    if(lohBreakend.getCluster() != lohCluster)
                        continue;

                    // walk towards the LOH from the other end of this DUP to see if it can find a resolving event within the cluster
                    List<SvBreakend> fullBreakendList = mState.getChrBreakendMap().get(lohBreakend.chromosome());

                    SvBreakend otherBreakend = lohBreakend.getOtherBreakend();
                    int index = otherBreakend.getChrPosIndex();
                    boolean traverseUp = otherBreakend.orientation() == -1;
                    SvCluster resolvingCluster = null;
                    SvBreakend resolvingBreakend = null;

                    while(true)
                    {
                        index += traverseUp ? 1 : -1;

                        if(index < 0 || index >= fullBreakendList.size())
                            break;

                        SvBreakend nextBreakend = fullBreakendList.get(index);

                        if(nextBreakend == lohBreakend)
                        {
                            // the LOH was reached without finding an offsetting SV
                            if(resolvingCluster == null)
                                break;

                            LNX_LOGGER.debug("cluster({}) SV({}) resolved prior to LOH by other cluster({}) breakend({})",
                                    lohCluster.id(), lohBreakend.getSV().posId(), resolvingCluster.id(), resolvingBreakend.toString());

                            addClusterReasons(lohBreakend.getSV(), resolvingBreakend.getSV(), LOH_CHAIN);

                            resolvingCluster.addClusterReason(LOH_CHAIN);
                            lohCluster.addClusterReason(LOH_CHAIN);

                            lohCluster.mergeOtherCluster(resolvingCluster);

                            mergedClusters.add(resolvingCluster);
                            break;
                        }

                        if(nextBreakend.orientation() == otherBreakend.orientation())
                            continue;

                        if(nextBreakend.isAssembledLink())
                            continue;

                        if(variantsHaveDifferentJcn(lohBreakend, nextBreakend))
                            continue;

                        SvCluster otherCluster = nextBreakend.getCluster();

                        if(otherCluster == lohCluster)
                            break; // own cluster resolves this LOH breakend

                        if(mergedClusters.contains(otherCluster) || otherCluster.isResolved())
                            continue;

                        if(resolvingCluster != null)
                            break; // cannot apply this rule if more than 1 cluster meet the conditions

                        // found an option, but continue on to see if any other clusters also satisfy the same conditions
                        resolvingBreakend = nextBreakend;
                        resolvingCluster = otherCluster;
                    }
                }

                ++lohIndex;
            }
        }

        if(mergedClusters.isEmpty())
            return false;

        mergedClusters.forEach(x -> clusters.remove(x));
        return true;
    }


    public boolean validateClustering(final List<SvCluster> clusters)
    {
        // validation that every SV was put into a cluster
        for (final Map.Entry<String, List<SvBreakend>> entry : mState.getChrBreakendMap().entrySet())
        {
            final List<SvBreakend> breakendList = entry.getValue();

            for (int i = 0; i < breakendList.size(); ++i)
            {
                final SvBreakend breakend = breakendList.get(i);
                SvVarData var = breakend.getSV();
                if(var.getCluster() == null)
                {
                    LNX_LOGGER.error("var({}) not clustered", var.posId());
                    return false;
                }
            }
        }

        // check that no 2 clusters contain the same SV
        for(int i = 0; i < clusters.size(); ++i)
        {
            SvCluster cluster1 = clusters.get(i);

            // check all SVs in this cluster reference it
            for(SvVarData var : cluster1.getSVs())
            {
                if(var.getCluster() != cluster1)
                {
                    LNX_LOGGER.error("var({}) in cluster({}) has incorrect ref", var.posId(), cluster1.id());
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
                        LNX_LOGGER.error("var({}) in 2 clusters({} and {})", var.posId(), cluster1.id(), cluster2.id());
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
                    LNX_LOGGER.error("cluster({}) exists twice in list", cluster1.id());
                    return false;
                }
            }
        }

        return true;
    }

}
