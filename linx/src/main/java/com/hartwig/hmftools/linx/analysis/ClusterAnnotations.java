package com.hartwig.hmftools.linx.analysis;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.ITEM_DELIM;
import static com.hartwig.hmftools.common.utils.Strings.appendStr;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.isStart;
import static com.hartwig.hmftools.linx.LinxConfig.LNX_LOGGER;
import static com.hartwig.hmftools.linx.analysis.ClusterClassification.isFilteredResolvedType;
import static com.hartwig.hmftools.linx.analysis.SvUtilities.calcConsistency;
import static com.hartwig.hmftools.linx.analysis.SvUtilities.copyNumbersEqual;
import static com.hartwig.hmftools.linx.analysis.SvUtilities.formatJcn;
import static com.hartwig.hmftools.linx.chaining.ChainMetrics.extractChainMetrics;
import static com.hartwig.hmftools.linx.chaining.LinkFinder.getMinTemplatedInsertionLength;
import static com.hartwig.hmftools.linx.types.LinkedPair.LOCATION_TYPE_EXTERNAL;
import static com.hartwig.hmftools.linx.types.LinkedPair.LOCATION_TYPE_INTERNAL;
import static com.hartwig.hmftools.linx.types.LinkedPair.LOCATION_TYPE_REMOTE;
import static com.hartwig.hmftools.linx.types.LinxConstants.MAX_MERGE_DISTANCE;
import static com.hartwig.hmftools.linx.types.LinxConstants.SHORT_TI_LENGTH;
import static com.hartwig.hmftools.linx.types.ResolvedType.COMPLEX;
import static com.hartwig.hmftools.linx.types.ResolvedType.LINE;
import static com.hartwig.hmftools.linx.types.ResolvedType.SIMPLE_GRP;
import static com.hartwig.hmftools.linx.types.SvCluster.CLUSTER_ANNOT_REP_REPAIR;
import static com.hartwig.hmftools.linx.types.SvCluster.CLUSTER_ANNOT_SHATTERING;

import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.linx.chaining.ChainMetrics;
import com.hartwig.hmftools.linx.chaining.SvChain;
import com.hartwig.hmftools.linx.types.ArmCluster;
import com.hartwig.hmftools.linx.types.ArmClusterType;
import com.hartwig.hmftools.linx.types.ArmGroup;
import com.hartwig.hmftools.linx.types.LinkedPair;
import com.hartwig.hmftools.linx.types.SvBreakend;
import com.hartwig.hmftools.linx.types.SvCluster;
import com.hartwig.hmftools.linx.types.SvVarData;

import org.apache.commons.math3.distribution.PoissonDistribution;

// post-clustering and chaining routines for annotating clusters, chains and links

public class ClusterAnnotations
{
    public static void annotateTemplatedInsertions(final SvCluster cluster, final Map<String, List<SvBreakend>> chrBreakendMap)
    {
        /* for each TI cache if:
            - it has a DB on either or both sides in the same cluster
            - local, the distance to the next SV in the cluster
            - is enclosed by the chain ends or remote if clear
            - how many other TIs it overlaps
         */

        if(cluster.getChains().isEmpty())
            return;

        // gather up start and end arms from each chain, to determine origin arms for the cluster
        final List<String> startEndArms = Lists.newArrayList();

        for(final SvChain chain : cluster.getChains())
        {
            for(int be1 = SE_START; be1 <= SE_END; ++be1)
            {
                final SvBreakend chainEnd = chain.getOpenBreakend(isStart(be1));
                if(chainEnd == null)
                    continue;

                if(!startEndArms.contains(chainEnd.getChrArm()))
                    startEndArms.add(chainEnd.getChrArm());
            }
        }

        for(final SvChain chain : cluster.getChains())
        {
            final SvBreakend chainStart = chain.getOpenBreakend(true);
            final SvBreakend chainEnd = chain.getOpenBreakend(false);

            SvBreakend chainLowerBe = null;
            SvBreakend chainUpperBe = null;
            boolean startEndSameArm = false;
            boolean chainEndsCNMatch = false;
            double chainEndsCN = 0;

            if(chainStart != null && chainEnd != null)
            {
                chainLowerBe = chainStart.position() < chainEnd.position() ? chainStart : chainEnd;
                chainUpperBe = chainStart == chainLowerBe ? chainEnd : chainStart;

                startEndSameArm = chainEnd.getChrArm().equals(chainStart.getChrArm());

                double chainStartCN = chainStart.copyNumber();
                double chainEndCN = chainEnd.copyNumber();
                chainEndsCNMatch = copyNumbersEqual(chainStartCN, chainEndCN);

                if(chainEndsCNMatch)
                    chainEndsCN = max(chainStartCN, chainEndCN);
            }

            for(final LinkedPair pair : chain.getLinkedPairs())
            {
                if(pair.first().isSglBreakend() || pair.second().isSglBreakend())
                    continue;

                SvBreakend pairLowerBe = pair.getBreakend(true);
                SvBreakend pairUpperBe = pair.getBreakend(false);
                double pairCN = getLinkPairCopyNumber(pair);
                boolean pairCNExceedsChainEnds = chainEndsCNMatch && pairCN > chainEndsCN && !copyNumbersEqual(pairCN, chainEndsCN);

                if(pairCNExceedsChainEnds)
                {
                    pair.setHasCopyNumberGain(true);
                }

                if(startEndSameArm)
                {
                    if(pairLowerBe.chromosome().equals(chainLowerBe.chromosome()))
                    {
                        // need to account for DBs
                        int lowerBuffer = getMinTemplatedInsertionLength(pairLowerBe, chainLowerBe);
                        int upperBuffer = getMinTemplatedInsertionLength(pairUpperBe, chainUpperBe);

                        if(pairLowerBe.position() >= chainLowerBe.position() - lowerBuffer
                        && pairUpperBe.position() <= chainUpperBe.position() + upperBuffer)
                        {
                            pair.setLocationType(LOCATION_TYPE_INTERNAL);
                        }
                        else
                        {
                            pair.setLocationType(LOCATION_TYPE_EXTERNAL);
                        }
                    }
                    else
                    {
                        pair.setLocationType(LOCATION_TYPE_REMOTE);
                    }
                }
                else
                {
                    if((chainStart != null && pairLowerBe.getChrArm().equals(chainStart.getChrArm()))
                    || (chainEnd != null && pairLowerBe.getChrArm().equals(chainEnd.getChrArm())))
                    {
                        pair.setLocationType(LOCATION_TYPE_EXTERNAL);
                    }
                    else
                    {
                        pair.setLocationType(LOCATION_TYPE_REMOTE);
                    }
                }

                List<LinkedPair> uniqueOverlaps = Lists.newArrayList();
                int overlapCount = 0;

                for(final LinkedPair otherPair : chain.getLinkedPairs())
                {
                    if(pair == otherPair)
                        continue;

                    if(uniqueOverlaps.stream().anyMatch(x -> x.matches(otherPair)))
                        continue;

                    if(!pair.chromosome().equals(otherPair.chromosome()) || pair.positionDistance() <= otherPair.positionDistance())
                        continue;

                    uniqueOverlaps.add(otherPair);

                    // make note of this overlap if the overlap distance exceeds the short TI length

                    if((pair.getBreakend(true).position() > otherPair.getBreakend(false).position())
                    || (pair.getBreakend(false).position() < otherPair.getBreakend(true).position()))
                    {
                        continue;
                    }

                    int overlapStart = max(pair.getBreakend(true).position(), otherPair.getBreakend(true).position());
                    int overlapEnd = min(pair.getBreakend(false).position(), otherPair.getBreakend(false).position());

                    if(overlapEnd - overlapStart >= SHORT_TI_LENGTH) // longer than a max DB length
                    {
                        ++overlapCount;
                    }
                }

                pair.setOverlapCount(overlapCount);

                // find closest SV in this cluster
                final List<SvBreakend> breakendList = chrBreakendMap.get(pair.chromosome());

                int[] nextSVData = getNextClusterSVData(cluster, breakendList, pair);
                pair.setNextSVData(nextSVData[NEXT_SV_DISTANCE], nextSVData[NEXT_CLUSTERED_SV_DISTANCE]);

                // how many SVs are traversed by this link
                pair.setTraversedSVCount(getTraversedSvCount(cluster, breakendList,
                        pair.getBreakend(true).getChrPosIndex(), pair.getBreakend(false).getChrPosIndex()));
            }
        }
    }

    private static double getLinkPairCopyNumber(final LinkedPair pair)
    {
        // skip past any DBs
        double[] breakendCN = {0, 0};

        for(int se = SE_START; se <= SE_END; ++se)
        {
            SvBreakend breakend = pair.getBreakend(se);

            if(breakend.getDBLink() != null && breakend.getDBLink().length() < 1)
                breakendCN[se] = max(breakend.copyNumber() - breakend.getDBLink().getOtherBreakend(breakend).jcn(), 0);
            else
                breakendCN[se] = breakend.copyNumber();
        }

        return (breakendCN[SE_START] + breakendCN[SE_END]) * 0.5;
    }

    private static int getTraversedSvCount(final SvCluster cluster, final List<SvBreakend> breakendList, int lowerIndex, int upperIndex)
    {
        String traversedInfo = getTraversedSvData(cluster, breakendList, lowerIndex, upperIndex);

        if(traversedInfo.isEmpty())
            return 0;

        String[] items = traversedInfo.split(ITEM_DELIM);
        return items.length;
    }

    private static String getTraversedSvData(final SvCluster cluster, final List<SvBreakend> breakendList, int lowerIndex, int upperIndex)
    {
        if(lowerIndex >= upperIndex - 1)
        {
            return "";
        }
        else if(lowerIndex < 0 || upperIndex >= breakendList.size())
        {
            LNX_LOGGER.error("invalid indices({} & {}) vs breakend list size({})", lowerIndex, upperIndex, breakendList.size());
            return "";
        }

        String traversedInfo = "";

        for(int i = lowerIndex + 1; i <= upperIndex - 1; ++i)
        {
            final SvBreakend breakend = breakendList.get(i);
            final SvCluster otherCluster = breakend.getCluster();

            if(otherCluster == cluster || otherCluster.isResolved())
                continue;

            traversedInfo = appendStr(traversedInfo,
                    String.format("%d %.2f", breakend.orientation(), breakend.getSV().copyNumberChange(breakend.usesStart())),
                    ';');
        }

        return traversedInfo;
    }

    private static final int NEXT_SV_DISTANCE = 0;
    private static final int NEXT_CLUSTERED_SV_DISTANCE = 1;

    private static int[] getNextClusterSVData(final SvCluster cluster, final List<SvBreakend> breakendList, final LinkedPair pair)
    {
        // walk forward and backwards from this pair to the closest SV in the same cluster
        // counting the number of non-trivial SVs traversed in the process
        int[] nextSvData = {-1, -1, 0};

        SvBreakend lowerBreakend = pair.getBreakend(true);
        SvBreakend upperBreakend = pair.getBreakend(false);
        int lowerIndex = lowerBreakend.getChrPosIndex();
        int upperIndex = upperBreakend.getChrPosIndex();

        if(lowerIndex > 0)
        {
            for(int i = lowerIndex - 1; i >= 0; --i)
            {
                final SvBreakend breakend = breakendList.get(i);
                int distance = lowerBreakend.position() - breakend.position();

                if(breakend.getCluster() == cluster)
                {
                    nextSvData[NEXT_CLUSTERED_SV_DISTANCE] = distance;
                    break;
                }
                else if(!isFilteredResolvedType(breakend.getCluster().getResolvedType()))
                {
                    if(nextSvData[NEXT_SV_DISTANCE] == -1)
                        nextSvData[NEXT_SV_DISTANCE] = distance;
                }
            }
        }

        if(upperIndex < breakendList.size() - 1)
        {
            for(int i = upperIndex + 1; i < breakendList.size(); ++i)
            {
                final SvBreakend breakend = breakendList.get(i);

                int distance = breakend.position() - upperBreakend.position();

                if(breakend.getCluster() == cluster)
                {
                    if(nextSvData[NEXT_CLUSTERED_SV_DISTANCE] == -1 || distance < nextSvData[NEXT_CLUSTERED_SV_DISTANCE])
                    {
                        nextSvData[NEXT_CLUSTERED_SV_DISTANCE] = distance;
                    }

                    break;
                }
                else if(!isFilteredResolvedType(breakend.getCluster().getResolvedType()))
                {
                    if(nextSvData[NEXT_SV_DISTANCE] == -1 || distance < nextSvData[NEXT_SV_DISTANCE])
                        nextSvData[NEXT_SV_DISTANCE] = distance;
                }
            }
        }

        if(nextSvData[NEXT_SV_DISTANCE] == -1)
            nextSvData[NEXT_SV_DISTANCE] = nextSvData[NEXT_CLUSTERED_SV_DISTANCE];

        return nextSvData;
    }

    public static void reportUnderclustering(final String sampleId, final List<SvCluster> clusters, Map<String, List<SvBreakend>> chrBreakendMap)
    {
        List<SvCluster> complexClusters = clusters.stream()
                .filter(x -> x.getSvCount() > 3)
                .filter(x -> !x.isResolved())
                .collect(Collectors.toList());

        List<SvCluster[]> reportedClusters = Lists.newArrayList();

        for(SvCluster cluster : complexClusters)
        {
            // calculate a total genomic span for this cluster by extending out 5M from each breakend in the cluster
            long genomicSpan = 0;

            Map<SvCluster, Integer> otherClusterMatches = Maps.newHashMap();

            for(final Map.Entry<String, List<SvBreakend>> entry : cluster.getChrBreakendMap().entrySet())
            {
                final List<SvBreakend> breakendList = entry.getValue();
                final List<SvBreakend> fullBreakendList = chrBreakendMap.get(entry.getKey());

                genomicSpan += MAX_MERGE_DISTANCE * 2; // account for outer breakends in the cluster

                for(int i = 0; i < breakendList.size() - 1; ++i)
                {
                    final SvBreakend breakend = breakendList.get(i);
                    final SvBreakend nextBreakend = breakendList.get(i + 1);

                    int breakendDistance = nextBreakend.position() - breakend.position();

                    genomicSpan += min(MAX_MERGE_DISTANCE * 2, breakendDistance);

                    // look for gaps in this cluster's breakend list
                    if(nextBreakend.getChrPosIndex() == breakend.getChrPosIndex() + 1)
                        continue;

                    List<SvCluster> encounteredClusters = Lists.newArrayList();

                    for(int j = breakend.getChrPosIndex() + 1; j < nextBreakend.getChrPosIndex(); ++j)
                    {
                        final SvBreakend otherBreakend = fullBreakendList.get(j);
                        final SvCluster otherCluster = otherBreakend.getCluster();

                        if(!complexClusters.contains(otherCluster) || encounteredClusters.contains(otherCluster))
                            continue;

                        int distance = min(
                                abs(breakend.position() - otherBreakend.position()),
                                abs(nextBreakend.position() - otherBreakend.position()));

                        if(distance > MAX_MERGE_DISTANCE)
                            continue;

                        encounteredClusters.add(otherBreakend.getCluster());
                    }

                    for(SvCluster otherCluster : encounteredClusters)
                    {
                        Integer matchCount = otherClusterMatches.get(otherCluster);
                        if(matchCount == null)
                            otherClusterMatches.put(otherCluster, 1);
                        else
                            otherClusterMatches.put(otherCluster, matchCount + 1);
                    }
                }

                // also walk from the ends in each direction
                for(int i = 0; i <= 1; ++i)
                {
                    boolean traverseUp = (i == 0);

                    final SvBreakend refBreakend = !traverseUp ? breakendList.get(0) : breakendList.get(breakendList.size() - 1);
                    int index = refBreakend.getChrPosIndex();

                    List<SvCluster> encounteredClusters = Lists.newArrayList();

                    while(true)
                    {
                        index += traverseUp ? 1 : -1;

                        if(index < 0 || index >= fullBreakendList.size())
                            break;

                        final SvBreakend otherBreakend = fullBreakendList.get(index);
                        final SvCluster otherCluster = otherBreakend.getCluster();

                        if(!complexClusters.contains(otherCluster) || encounteredClusters.contains(otherCluster))
                            continue;

                        if(abs(otherBreakend.position() - refBreakend.position()) > MAX_MERGE_DISTANCE)
                            break;

                        encounteredClusters.add(otherBreakend.getCluster());
                    }

                    for(SvCluster otherCluster : encounteredClusters)
                    {
                        Integer matchCount = otherClusterMatches.get(otherCluster);
                        if(matchCount == null)
                            otherClusterMatches.put(otherCluster, 1);
                        else
                            otherClusterMatches.put(otherCluster, matchCount + 1);
                    }
                }
            }

            for(final Map.Entry<SvCluster, Integer> entry : otherClusterMatches.entrySet())
            {
                if(entry.getValue() <= 1)
                    continue;

                final SvCluster otherCluster = entry.getKey();

                if(reportedClusters.stream().anyMatch(x -> x[0] == otherCluster && x[1] == cluster))
                    continue;

                reportedClusters.add(new SvCluster[]{cluster, otherCluster});

                // ignore hom-loss and LOH-related clusters
                boolean lossRelated = cluster.getLohEvents().stream()
                        .anyMatch(x -> x.getHomLossEvents().stream()
                        .filter(y -> y.matchedBothSVs() && y.getBreakend(true).getCluster() == otherCluster).count() > 0);

                if(!lossRelated)
                {
                    lossRelated = otherCluster.getLohEvents().stream()
                            .anyMatch(x -> x.getHomLossEvents().stream()
                                    .filter(y -> y.matchedBothSVs() && y.getBreakend(true).getCluster() == cluster).count() > 0);
                }

                int matchCount = entry.getValue();

                double hitProbability = genomicSpan / 3e9;
                int otherClusterBreakendCount = otherCluster.getSvCount() * 2 - otherCluster.getSglBreakendCount();
                double expectedHits = hitProbability * otherClusterBreakendCount;

                if(expectedHits >= matchCount)
                    continue;

                PoissonDistribution poissonDist = new PoissonDistribution(expectedHits);
                double poissonProb = 1 - poissonDist.cumulativeProbability(matchCount - 1);

                if(poissonProb < 0.001)
                {
                    String overlappingChrStr = "";

                    for(final ArmGroup group : cluster.getArmGroups())
                    {
                        if(otherCluster.getArmGroups().stream().anyMatch(x -> x.id().equals(group.id())))
                        {
                            overlappingChrStr = appendStr(overlappingChrStr, group.id(), ';');
                        }
                    }

                    LNX_LOGGER.info("sampleId({}) cluster({} count={}) proxCount({}) with cluster({} count={}) prob({}) chromosomes({}) lossRelated({})",
                            sampleId, cluster.id(), cluster.getSvCount(), matchCount, otherCluster.id(), otherCluster.getSvCount(),
                            String.format("%.6f exp=%.2f span=%.1fM", poissonProb, expectedHits, genomicSpan/1e6),
                            overlappingChrStr, lossRelated);
                }
            }
        }
    }

    public static void annotateClusterDeletions(final SvCluster cluster, final Map<String, List<SvBreakend>> chrBreakendMap)
    {
        if(cluster.getSvCount() == 1 || cluster.getResolvedType() == LINE)
            return;
        cluster.getMetrics().populate(cluster, chrBreakendMap);
    }

    private static final int INCONSISTENT_ARM = -2;
    private static final int CONSISTENT_ARM = 0;

    public static void annotateClusterChains(final SvCluster cluster)
    {
        // skip simple SVs and single arm clusters without chains
        if(cluster.getSvCount() == 1 || (cluster.getArmCount() == 1 && cluster.getChains().isEmpty()))
            return;

        /* data to gather for each arm in the chain
            - number of links
            - number of short TIs without proximate deletion bridges
            - links with copy number gain
            - start and end locations
         */

        List<String> fragmentArms = Lists.newArrayList();
        Map<String,Integer> originArms = Maps.newHashMap(); // chr-arm and consistency (1/-1, 0 if matched, or other if invalid)

        int chainCount = cluster.getChains().size();
        int inconsistentChains = 0;
        int shortTiCount = 0;
        int longTiCount = 0;

        for(final SvChain chain : cluster.getChains())
        {
            boolean chainConsistent = chain.isConsistent();

            cluster.getMetrics().ChainedLength += chain.getLength(false);

            for(final LinkedPair pair : chain.getLinkedPairs())
            {
                final SvVarData first = pair.first();

                if(pair.first().isSglBreakend() || pair.second().isSglBreakend())
                    continue;

                final String chrArm = first.getBreakend(pair.firstLinkOnStart()).getChrArm();

                if(pair.baseLength() <= SHORT_TI_LENGTH)
                {
                    ++shortTiCount;

                    if(!fragmentArms.contains(chrArm))
                        fragmentArms.add(chrArm);
                }
                else
                {
                    ++longTiCount;
                }
            }

            if(!chainConsistent)
            {
                ++inconsistentChains;
            }

            for(int se = SE_START; se <= SE_END; ++se)
            {
                final SvBreakend chainBreakend = chain.getOpenBreakend(se);

                if(chainBreakend != null)
                {
                    addBreakendToArmConsistency(originArms, chainBreakend);
                }
            }
        }

        final List<SvVarData> unlinkedSVs = cluster.getUnlinkedSVs();

        for(final SvVarData var : unlinkedSVs)
        {
            for(int se = SE_START; se <= SE_END; ++se)
            {
                if(se == SE_END && var.isSglBreakend())
                    continue;

                addBreakendToArmConsistency(originArms, var.getBreakend(se));
            }
        }

        fragmentArms = fragmentArms.stream().filter(x -> !originArms.containsKey(x)).collect(Collectors.toList());

        int consistentArmCount = (int)originArms.values().stream().filter(x -> x >= -1 && x <= 1).count();

        // determine complex / active arms as those with more than 1 local topology event not including TI-only segments
        Map<String,Integer> armClusterCounts = Maps.newHashMap();

        for(final ArmCluster armCluster : cluster.getArmClusters())
        {
            if(armCluster.getType() == ArmClusterType.TI_ONLY)
                continue;

            Integer count = armClusterCounts.get(armCluster.getChrArm());

            if(count == null)
            {
                armClusterCounts.put(armCluster.getChrArm(), 1);
            }
            else
            {
                armClusterCounts.put(armCluster.getChrArm(), count + 1);
            }
        }

        int complexArms = (int)armClusterCounts.values().stream().filter(x -> x > 1).count();

        ClusterMetrics metrics = cluster.getMetrics();
        metrics.OriginArms = originArms.size();
        metrics.FragmentArms = fragmentArms.size();
        metrics.ConsistentArms = consistentArmCount;
        metrics.ComplexArms = complexArms;

        if(cluster.getSvCount() > 2)
        {
            boolean isComplex = cluster.requiresReplication() || !cluster.getFoldbacks().isEmpty();
            boolean isComplete = (inconsistentChains == 0) && (consistentArmCount == originArms.size());

            LNX_LOGGER.debug("cluster({}) {} chains({} incons={}) arms({} frag={} origin={} consis={}) tiCount(short={} long={})",
                    cluster.id(), isComplete ? "COMPLETE" : "incomplete", chainCount, inconsistentChains,
                    cluster.getArmGroups().size(), fragmentArms.size(), originArms.size(), consistentArmCount,
                    shortTiCount, longTiCount);

            if(isComplete && !isComplex && longTiCount > 0 && cluster.getResolvedType() == COMPLEX)
            {
                // chromothripsis is currently defined as fully chained simple cluster
                // but needs to take into account the copy number gain / loss compared with the surrounding chromatid

                cluster.addAnnotation(CLUSTER_ANNOT_SHATTERING);
            }
        }
    }

    private static void addBreakendToArmConsistency(final Map<String,Integer> armConsistencyMap, final SvBreakend breakend)
    {
        final String chrArm = breakend.getChrArm();

        Integer armConsistency = armConsistencyMap.get(chrArm);
        if(armConsistency == null)
        {
            armConsistencyMap.put(chrArm, calcConsistency(breakend));
            return;
        }

        if(armConsistency == INCONSISTENT_ARM)
            return;

        if(armConsistency == CONSISTENT_ARM || calcConsistency(breakend) == armConsistency)
        {
            armConsistencyMap.put(chrArm, INCONSISTENT_ARM);
        }
        else
        {
            armConsistencyMap.put(chrArm, CONSISTENT_ARM);
        }
    }

    public static void annotateReplicationBeforeRepair(final SvCluster cluster)
    {
        if(cluster.getSvCount() == 1 || cluster.getSvCount() > 100)
            return;

        if(cluster.getResolvedType() != COMPLEX && cluster.getResolvedType() != SIMPLE_GRP)
            return;

        // criteria:
        // - uniform JCN within a chain
        // - at least one overlapping TI
        // - copy number gain between the chain ends, not explained by another chain or SV

        for(SvChain chain : cluster.getChains())
        {
            if(!chain.isConsistent() || chain.isClosedLoop())
                continue;

            final SvBreakend chainStart = chain.getOpenBreakend(true);
            final SvBreakend chainEnd = chain.getOpenBreakend(false);

            if(chainStart == null || chainEnd == null || !chainEnd.getChrArm().equals(chainStart.getChrArm()))
                continue;

            // check for any repeated SV
            if(chain.hasRepeatedSV())
                continue;

            final ChainMetrics chainMetrics = extractChainMetrics(chain);

            // require internal TIs with gain, which mandates chain ends are on the same arm
            if(chainMetrics.ChainEndsAway != 1 || chainMetrics.InternalTICnGain == 0 || chainMetrics.OverlappingTIs == 0)
                continue;

            long internalGainLength = getCopyNumberGainLength(cluster, chain);

            if(internalGainLength == 0)
                continue;

            long chainRange = abs(chainEnd.position() - chainStart.position());
            double chainEndsCN = max(chainStart.copyNumber(), chainEnd.copyNumber());
            double gainPercent = internalGainLength / (double) chainRange;

            if(gainPercent >= 0.05 && chainRange >= 1000)
            {
                cluster.addAnnotation(CLUSTER_ANNOT_REP_REPAIR);

                LNX_LOGGER.debug("cluster({} desc={} type={}) chain({} cn={}) rep-repair: links({} internal={} withGain={} overlaps={}) length(chain={} gain={} perc={})",
                        cluster.id(), cluster.getDesc(), cluster.getResolvedType(), chain.id(), formatJcn(chainEndsCN),
                        chain.getLinkCount(), chainMetrics.InternalTIs, chainMetrics.InternalTICnGain, chainMetrics.OverlappingTIs,
                        chainRange, internalGainLength, String.format("%.2f", gainPercent));
            }
        }
    }

    private static long getCopyNumberGainLength(final SvCluster cluster, final SvChain chain)
    {
        final SvBreakend chainStart = chain.getOpenBreakend(true);
        final SvBreakend chainEnd = chain.getOpenBreakend(false);
        int lowerIndex = min(chainStart.getClusterChrPosIndex(), chainEnd.getClusterChrPosIndex());
        int upperIndex = max(chainStart.getClusterChrPosIndex(), chainEnd.getClusterChrPosIndex());
        double chainEndsCN = max(chainStart.copyNumber(), chainEnd.copyNumber());
        long internalGainLength = 0;

        // sum of segments of CN gain between the chain ends
        final List<SvBreakend> breakendList = cluster.getChrBreakendMap().get(chainStart.chromosome());

        int prevPosition = 0;
        double prevCN = 0;
        boolean inSegment = false;
        double netPloidy = 0;

        for(int i = lowerIndex + 1; i < upperIndex - 1; ++i)
        {
            SvBreakend breakend = breakendList.get(i);

            // if a breakend from another chain is encountered then exit due to the uncertainty around chaining
            if(!chain.getSvList().contains(breakend.getSV()))
                return 0;

            if(!inSegment)
            {
                if(breakend.orientation() == 1)
                {
                    if(prevPosition > 0)
                    {
                        if(breakend.copyNumber() > chainEndsCN && !copyNumbersEqual(breakend.copyNumber(), chainEndsCN))
                        {
                            internalGainLength += breakend.position() - prevPosition;
                        }
                    }

                    continue;
                }

                inSegment = true;
                prevPosition = breakend.position();
                prevCN = breakend.copyNumber();
                netPloidy = breakend.jcn();
            }
            else
            {
                int segmentLength = breakend.position() - prevPosition;

                if(breakend.orientation() == -1)
                {
                    // another breakend increasing CN - record segment CN up to this point
                    if(prevCN > chainEndsCN && !copyNumbersEqual(prevCN, chainEndsCN))
                    {
                        internalGainLength += segmentLength;
                    }

                    // move position on
                    prevPosition = breakend.position();
                    prevCN = breakend.copyNumber();
                    netPloidy += breakend.jcn();
                }
                else
                {
                    if(breakend.copyNumber() > chainEndsCN && !copyNumbersEqual(breakend.copyNumber(), chainEndsCN))
                    {
                        internalGainLength += segmentLength;
                    }

                    netPloidy -= breakend.jcn();
                    prevPosition = breakend.position();

                    if(copyNumbersEqual(netPloidy, 0))
                        inSegment = false;
                }
            }
        }

        return internalGainLength;
    }

}
