package com.hartwig.hmftools.linx.analysis;

import static java.lang.Math.abs;

import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.isStart;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.SGL;
import static com.hartwig.hmftools.linx.LinxConfig.LNX_LOGGER;
import static com.hartwig.hmftools.linx.analysis.ClusteringReason.COMMON_ARMS;
import static com.hartwig.hmftools.linx.analysis.ClusteringReason.CONSEC_BREAKS;
import static com.hartwig.hmftools.linx.analysis.ClusteringReason.FOLDBACKS;
import static com.hartwig.hmftools.linx.analysis.ClusteringReason.OVERLAP_FOLDBACKS;
import static com.hartwig.hmftools.linx.analysis.ClusteringReason.SATELLITE_SGL;
import static com.hartwig.hmftools.linx.analysis.ClusteringReason.TI_JCN_MATCH;
import static com.hartwig.hmftools.linx.analysis.SimpleClustering.skipClusterType;
import static com.hartwig.hmftools.linx.analysis.SimpleClustering.variantsHaveDifferentJcn;
import static com.hartwig.hmftools.linx.analysis.SimpleClustering.variantsViolateLohHomLoss;
import static com.hartwig.hmftools.linx.analysis.SvUtilities.copyNumbersEqual;
import static com.hartwig.hmftools.linx.analysis.SvUtilities.makeChrArmStr;
import static com.hartwig.hmftools.linx.types.LinxConstants.MAX_MERGE_DISTANCE;
import static com.hartwig.hmftools.linx.types.SvVarData.haveSameChrArms;

import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.linx.chaining.ChainJcnLimits;
import com.hartwig.hmftools.linx.types.SvBreakend;
import com.hartwig.hmftools.linx.types.SvCluster;
import com.hartwig.hmftools.linx.types.SvVarData;

public class ComplexClustering
{
    // references only
    private final List<SvCluster> mClusters;
    private final ClusteringState mState;
    private final SimpleClustering mSimpleClustering;
    private String mSampleId;

    public ComplexClustering(final ClusteringState state, final List<SvCluster> clusters, final SimpleClustering simpleClustering)
    {
        mClusters = clusters;
        mState = state;
        mSimpleClustering = simpleClustering;
        mSampleId = "";
    }

    public void applyRules(final String sampleId, boolean foldbacksOnly)
    {
        mSampleId = sampleId;

        // second round of cluster merging on more complex criteria and inconsistencies:
        // merge on foldbacks on the same arm
        // merge on links between common arms
        // merge if one cluster has footprints which overlap unresolved complex SVs
        // merge clusters which resolve another's LOH DUP

        // first collect the clusters for which these complex rules apply
        List<SvCluster> complexClusters = mClusters.stream()
                .filter(x -> !x.isResolved())
                .collect(Collectors.toList());

        int iterations = 1;
        boolean foundMerges = true;

        while(foundMerges)
        {
            foundMerges = false;

            int index1 = 0;
            while(index1 < complexClusters.size())
            {
                SvCluster cluster1 = complexClusters.get(index1);

                if(cluster1.isResolved())
                {
                    ++index1;
                    continue;
                }

                int index2 = index1 + 1;
                while(index2 < complexClusters.size())
                {
                    SvCluster cluster2 = complexClusters.get(index2);

                    if(cluster2.isResolved())
                    {
                        ++index2;
                        continue;
                    }

                    // try each merge reason in turn
                    boolean canMergeClusters = false;

                    canMergeClusters = canMergeClustersOnFoldbacks(cluster1, cluster2);

                    if(!canMergeClusters && !foldbacksOnly)
                        canMergeClusters = canMergeClustersOnCommonArms(cluster1, cluster2);

                    if(!canMergeClusters)
                    {
                        ++index2;
                        continue;
                    }

                    foundMerges = true;
                    cluster1.mergeOtherCluster(cluster2, false);

                    complexClusters.remove(index2);
                    mClusters.remove(cluster2);
                }

                ++index1;
            }

            if(!foldbacksOnly)
            {
                if(mergeBreakendStraddledClusters(complexClusters))
                    foundMerges = true;

                if(mergeFacingJcnLinkedClusters(complexClusters))
                    foundMerges = true;

                if(mergeSingleSatelliteRepeats(complexClusters))
                    foundMerges = true;
            }

            ++iterations;

            if(iterations > 20)
            {
                LNX_LOGGER.warn("sample({}) reached {} iterations of clustering merging", mSampleId, iterations);
                break;
            }
        }
    }

    private boolean canMergeClustersOnFoldbacks(final SvCluster cluster1, final SvCluster cluster2)
    {
        // merge any clusters with foldbacks on the same arm
        final List<SvVarData> cluster1Foldbacks = cluster1.getFoldbacks();
        final List<SvVarData> cluster2Foldbacks = cluster2.getFoldbacks();

        if (cluster1Foldbacks.isEmpty() && cluster2Foldbacks.isEmpty())
            return false;

        if (!cluster1Foldbacks.isEmpty() && !cluster2Foldbacks.isEmpty())
        {
            for (final SvVarData var1 : cluster1Foldbacks)
            {
                for (int be1 = SE_START; be1 <= SE_END; ++be1)
                {
                    boolean v1Start = isStart(be1);

                    if (be1 == SE_END && var1.isSglBreakend())
                        continue;

                    if (var1.getFoldbackBreakend(v1Start) == null)
                        continue;

                    for (final SvVarData var2 : cluster2Foldbacks)
                    {
                        for (int be2 = SE_START; be2 <= SE_END; ++be2)
                        {
                            boolean v2Start = isStart(be2);

                            if (be2 == SE_END && var2.isSglBreakend())
                                continue;

                            if (var2.getFoldbackBreakend(v2Start) == null)
                                continue;

                            if (!var1.chromosome(v1Start).equals(var2.chromosome(v2Start)) || !var1.arm(v1Start).equals(var2.arm(v2Start)))
                                continue;

                            if(variantsHaveDifferentJcn(var1, var2))
                                continue;

                            LNX_LOGGER.debug("cluster({}) SV({}) and cluster({}) SV({}) have foldbacks on same arm",
                                    cluster1.id(), var1.posId(), cluster2.id(), var2.posId());

                            mSimpleClustering.addClusterReasons(var1, var2, FOLDBACKS);

                            cluster1.addClusterReason(FOLDBACKS);
                            cluster2.addClusterReason(FOLDBACKS);
                            return true;
                        }
                    }
                }
            }
        }

        return false;
    }

    private boolean canMergeClustersOnCommonArms(final SvCluster cluster1, final SvCluster cluster2)
    {
        // merge if the 2 clusters have BNDs linking the same 2 inconsistent or long arms
        // re-check which BNDs may link arms
        cluster1.setArmLinks();
        cluster2.setArmLinks();

        final List<SvVarData> crossArmList1 = cluster1.getUnlinkedRemoteSVs();
        final List<SvVarData> crossArmList2 = cluster2.getUnlinkedRemoteSVs();

        // now that the candidate arm groups have been established, just need to find a single BND
        // from each cluster that falls into the same par of arm groups

        for (final SvVarData var1 : crossArmList1)
        {
            for (final SvVarData var2 : crossArmList2)
            {
                if (!haveSameChrArms(var1, var2))
                    continue;

                if (variantsViolateLohHomLoss(var1, var2))
                    continue;

                if(variantsHaveDifferentJcn(var1, var2))
                    continue;

                LNX_LOGGER.debug("cluster({}) and cluster({}) have common links with SV({}) and SV({})",
                        cluster1.id(), cluster2.id(), var1.posId(), var2.posId());

                mSimpleClustering.addClusterReasons(var1, var2, COMMON_ARMS);

                cluster1.addClusterReason(COMMON_ARMS);
                cluster2.addClusterReason(COMMON_ARMS);
                return true;
            }
        }

        return false;
    }

    private static final int MAX_COMPLEX_CLUSTER_MERGE_SIZE = 5;

    private boolean mergeBreakendStraddledClusters(List<SvCluster> clusters)
    {
        // Breakends straddled by consecutive same orientation breakends
        // Merge any non resolved breakend to a cluster which straddles it immediately on both sides with 2 breakends
        // facing the same direction, and where the facing breakends have matching ploidy

        // OR

        // Breakends straddled by foldbacks
        // Merge any non resolved breakend into a cluster which has 2 different foldbacks straddling it immediately on
        // both sides and at least one of the foldbacks faces the breakend

        List<SvCluster> mergedClusters = Lists.newArrayList();

        int clusterIndex = 0;
        while(clusterIndex < clusters.size())
        {
            SvCluster cluster = clusters.get(clusterIndex);

            if(mergedClusters.contains(cluster) || cluster.isResolved())
            {
                ++clusterIndex;
                continue;
            }

            boolean mergedOtherClusters = false;

            for (final Map.Entry<String, List<SvBreakend>> entry : cluster.getChrBreakendMap().entrySet())
            {
                List<SvBreakend> breakendList = entry.getValue();

                List<SvBreakend> fullBreakendList = mState.getChrBreakendMap().get(entry.getKey());

                for (int i = 0; i < breakendList.size() - 1; ++i)
                {
                    final SvBreakend lowerBreakend = breakendList.get(i);
                    SvBreakend upperBreakend = breakendList.get(i + 1);

                    if (lowerBreakend.arm() != upperBreakend.arm())
                        continue;

                    // check for 2 different foldbacks with at least one facing in
                    boolean isFoldbackPair = lowerBreakend.isFoldback() && upperBreakend.isFoldback()
                            && lowerBreakend.getFoldbackBreakend() != upperBreakend
                            && !(lowerBreakend.orientation() == 1 && upperBreakend.orientation() == -1);

                    if (!isFoldbackPair)
                    {
                        if(lowerBreakend.orientation() == upperBreakend.orientation())
                        {
                            // exclude this pair if the front breakend is a DB
                            final SvBreakend frontBE = lowerBreakend.orientation() == 1 ? lowerBreakend : upperBreakend;

                            if (frontBE.getDBLink() != null && frontBE.getDBLink().length() < 0)
                                continue;

                            // all ok
                        }
                        else
                        {
                            if (i >= breakendList.size() - 2 || breakendList.get(i + 2).orientation() != lowerBreakend.orientation())
                                continue;

                            // include this pair if the back breakend is consecutive but masked by a DB
                            final SvBreakend backBE = lowerBreakend.orientation() == -1 ? lowerBreakend : upperBreakend;

                            if (backBE.getDBLink() != null && backBE.getDBLink().length() < 0)
                            {
                                upperBreakend = breakendList.get(i + 2);
                            }
                            else
                            {
                                continue;
                            }
                        }

                        if(!ChainJcnLimits.jcnMatch(lowerBreakend, upperBreakend))
                            continue;
                    }

                    // now look for a lone breakend which falls between these 2 consecutive breakends
                    // if there is a DB in between the FB or consecutive breakends then there will be a gap of 2 breakends in between
                    // ie index 0 (lower), 1 (other), 2 (DB) and 3 (upper)
                    int chrIndexLower = lowerBreakend.getChrPosIndex();
                    int chrIndexUpper = upperBreakend.getChrPosIndex();
                    SvBreakend otherBreakend = null;

                    int unclusteredStraddledBreakends = 0;

                    for (int j = chrIndexLower + 1; j <= chrIndexUpper - 1; ++j)
                    {
                        final SvBreakend breakend = fullBreakendList.get(j);
                        final SvCluster otherCluster = breakend.getCluster();

                        if(otherCluster.isResolved() || otherCluster == cluster)
                            continue;

                        ++unclusteredStraddledBreakends;
                        otherBreakend = breakend;
                    }

                    if(unclusteredStraddledBreakends != 1)
                        continue;

                    final SvCluster otherCluster = otherBreakend.getCluster();

                    if (otherCluster == cluster || otherCluster.isResolved() || mergedClusters.contains(otherCluster))
                        continue;

                    if(cluster.getSvCount() > MAX_COMPLEX_CLUSTER_MERGE_SIZE && otherCluster.getSvCount() > MAX_COMPLEX_CLUSTER_MERGE_SIZE)
                        continue;

                    // if not straddled by a foldback pair, then the breakend must be facing the consecutive straddling breakends
                    if(!isFoldbackPair && otherBreakend.orientation() == lowerBreakend.orientation())
                        continue;

                    LNX_LOGGER.debug("cluster({}) {} breakends({} & {}) overlap cluster({}) breakend({})",
                            cluster.id(), isFoldbackPair ? "foldback" : "consecutive",
                            lowerBreakend.toString(), upperBreakend.toString(), otherCluster.id(), otherBreakend.toString());

                    final ClusteringReason reason = isFoldbackPair ? OVERLAP_FOLDBACKS : CONSEC_BREAKS;
                    mSimpleClustering.addClusterReasons(otherBreakend.getSV(), lowerBreakend.getSV(), reason);

                    otherCluster.addClusterReason(reason);
                    cluster.addClusterReason(reason);

                    cluster.mergeOtherCluster(otherCluster);

                    mergedClusters.add(otherCluster);

                    mergedOtherClusters = true;
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
        mergedClusters.forEach(x -> mClusters.remove(x));
        return true;
    }

    private boolean mergeFacingJcnLinkedClusters(List<SvCluster> clusters)
    {
        // Extended proximity for complex and incomplete events
        // Merge any neighbouring (excluding resolved events) remaining incomplete or complex clusters that are within 5M bases and which
        // have facing flanking breakends on each cluster which could form a TI with matching ploidy.
        // In the case of a foldback the ploidy of the facing breakend is also permitted to match 2x the ploidy.
        List<SvCluster> mergedClusters = Lists.newArrayList();

        int clusterIndex = 0;
        while(clusterIndex < clusters.size())
        {
            SvCluster cluster = clusters.get(clusterIndex);

            if(mergedClusters.contains(cluster) || cluster.isResolved())
            {
                ++clusterIndex;
                continue;
            }

            boolean mergedOtherClusters = false;

            for (final Map.Entry<String, List<SvBreakend>> entry : cluster.getChrBreakendMap().entrySet())
            {
                List<SvBreakend> breakendList = entry.getValue();

                List<SvBreakend> fullBreakendList = mState.getChrBreakendMap().get(entry.getKey());

                for(int i = 0; i <= 1; ++i)
                {
                    boolean traverseUp = (i == 0);

                    SvBreakend boundaryBreakend = traverseUp ? breakendList.get(breakendList.size() - 1) : breakendList.get(0);

                    // breakend needs to face out
                    if(traverseUp == (boundaryBreakend.orientation() == 1))
                        continue;

                    // walk up from here to the next unresolved breakend
                    int index = boundaryBreakend.getChrPosIndex();

                    while(true)
                    {
                        index += traverseUp ? 1 : -1;

                        if(index < 0 || index >= fullBreakendList.size())
                            break;

                        final SvBreakend nextBreakend = fullBreakendList.get(index);

                        if(SimpleClustering.variantsHaveDifferentJcn(boundaryBreakend, nextBreakend))
                            continue;

                        if(abs(nextBreakend.position() - boundaryBreakend.position()) > MAX_MERGE_DISTANCE)
                            break;

                        final SvCluster otherCluster = nextBreakend.getCluster();

                        if(otherCluster == cluster || skipClusterType(otherCluster))
                            continue;

                        if(cluster.getSvCount() > MAX_COMPLEX_CLUSTER_MERGE_SIZE && otherCluster.getSvCount() > MAX_COMPLEX_CLUSTER_MERGE_SIZE)
                            continue;

                        if(nextBreakend.orientation() == boundaryBreakend.orientation())
                            break;

                        if(copyNumbersEqual(boundaryBreakend.jcn(), nextBreakend.jcn()))
                        {
                            LNX_LOGGER.debug("cluster({}) boundary breakend({}) ploidy TI match with cluster({}) breakend({})",
                                    cluster.id(), boundaryBreakend, otherCluster.id(), nextBreakend);

                            mSimpleClustering.addClusterReasons(boundaryBreakend.getSV(), nextBreakend.getSV(), TI_JCN_MATCH);

                            otherCluster.addClusterReason(TI_JCN_MATCH);
                            cluster.addClusterReason(TI_JCN_MATCH);

                            cluster.mergeOtherCluster(otherCluster);

                            mergedClusters.add(otherCluster);
                            mergedOtherClusters = true;
                        }

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
        mergedClusters.forEach(x -> mClusters.remove(x));
        return true;
    }

    private boolean mergeSingleSatelliteRepeats(List<SvCluster> clusters)
    {
        // merge any cluster with less than or 1 non SGL and non INF breakend with any other cluster which contains a SGL on the same
        // chromosome with matching repeat class or type marked as satellite

        // To protect against false positives and joining complex clusters which both touch repeats, but otherwise don’t appear to overlap,
        // we avoid clustering 2 clusters which already have multiple non SGL breakends
        List<SvCluster> clustersWithSatelliteRepeats = clusters.stream()
                .filter(x -> x.getSVs().stream().anyMatch(y -> y.sglToSatelliteRepeats()))
                .collect(Collectors.toList());

        if(clustersWithSatelliteRepeats.isEmpty())
            return false;

        List<SvCluster> clustersWithOneOrNoSgls = clustersWithSatelliteRepeats.stream()
                .filter(x -> x.getTypeCount(SGL) <= 1)
                .collect(Collectors.toList());

        List<SvCluster> mergedClusters = Lists.newArrayList();

        for(SvCluster srCluster : clustersWithSatelliteRepeats)
        {
            if(mergedClusters.contains(srCluster))
                continue;

            List<String> satelliteChrArms1 = srCluster.getSVs().stream()
                    .filter(x -> x.sglToSatelliteRepeats())
                    .map(x -> makeChrArmStr(x.chromosome(true), x.arm(true)))
                    .collect(Collectors.toList());

            int index = 0;
            while(index < clustersWithOneOrNoSgls.size())
            {
                SvCluster sglCluster = clustersWithOneOrNoSgls.get(index);

                if(sglCluster == srCluster)
                {
                    ++index;
                    continue;
                }


                boolean merged = false;

                for(SvVarData var : sglCluster.getSVs())
                {
                    if(var.sglToSatelliteRepeats() && satelliteChrArms1.contains(makeChrArmStr(var.chromosome(true), var.arm(true))))
                    {
                        final String chromosome = var.chromosome(true);

                        // find the other linking SGL
                        final SvVarData otherSV = srCluster.getSVs().stream()
                                .filter(x -> x.sglToSatelliteRepeats() && x.chromosome(true).equals(chromosome))
                                .findFirst().get();

                        if(variantsHaveDifferentJcn(var, otherSV))
                            continue;

                        LNX_LOGGER.debug("cluster({}) has same chromosome({}) link with satellite cluster({}) SV({})",
                                srCluster.id(), chromosome, sglCluster.id(), var.id());

                        mSimpleClustering.addClusterReasons(otherSV, var, SATELLITE_SGL);
                        srCluster.addClusterReason(SATELLITE_SGL);

                        srCluster.mergeOtherCluster(sglCluster);
                        mergedClusters.add(sglCluster);
                        merged = true;
                        break;
                    }
                }

                if(merged)
                    clustersWithOneOrNoSgls.remove(index);
                else
                    ++index;
            }
        }

        if(mergedClusters.isEmpty())
            return false;

        mergedClusters.forEach(x -> clusters.remove(x));
        mergedClusters.forEach(x -> mClusters.remove(x));
        return true;
    }

}
