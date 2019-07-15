package com.hartwig.hmftools.linx.analysis;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.SGL;
import static com.hartwig.hmftools.linx.chaining.LinkFinder.getMinTemplatedInsertionLength;
import static com.hartwig.hmftools.linx.analysis.SvClassification.isFilteredResolvedType;
import static com.hartwig.hmftools.linx.analysis.SvUtilities.CHROMOSOME_ARM_P;
import static com.hartwig.hmftools.linx.analysis.SvUtilities.CHROMOSOME_ARM_Q;
import static com.hartwig.hmftools.linx.analysis.SvUtilities.appendStr;
import static com.hartwig.hmftools.linx.analysis.SvUtilities.calcConsistency;
import static com.hartwig.hmftools.linx.analysis.SvUtilities.copyNumbersEqual;
import static com.hartwig.hmftools.linx.cn.CnDataLoader.CN_SEG_DATA_CN_AFTER;
import static com.hartwig.hmftools.linx.cn.CnDataLoader.CN_SEG_DATA_CN_BEFORE;
import static com.hartwig.hmftools.linx.cn.CnDataLoader.CN_SEG_DATA_MAP_AFTER;
import static com.hartwig.hmftools.linx.cn.CnDataLoader.CN_SEG_DATA_MAP_BEFORE;
import static com.hartwig.hmftools.linx.types.SvBreakend.DIRECTION_CENTROMERE;
import static com.hartwig.hmftools.linx.types.SvCluster.CLUSTER_ANNOT_SHATTERING;
import static com.hartwig.hmftools.linx.types.SvCluster.isSpecificCluster;
import static com.hartwig.hmftools.linx.types.SvLinkedPair.LOCATION_TYPE_EXTERNAL;
import static com.hartwig.hmftools.linx.types.SvLinkedPair.LOCATION_TYPE_INTERNAL;
import static com.hartwig.hmftools.linx.types.SvLinkedPair.LOCATION_TYPE_REMOTE;
import static com.hartwig.hmftools.linx.types.SvVarData.SE_END;
import static com.hartwig.hmftools.linx.types.SvVarData.SE_START;
import static com.hartwig.hmftools.linx.types.SvVarData.isStart;
import static com.hartwig.hmftools.linx.types.SvaConstants.NO_DB_MARKER;
import static com.hartwig.hmftools.linx.types.SvaConstants.SHORT_TI_LENGTH;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.linx.cn.CnDataLoader;
import com.hartwig.hmftools.linx.types.SvArmGroup;
import com.hartwig.hmftools.linx.types.SvBreakend;
import com.hartwig.hmftools.linx.cn.SvCNData;
import com.hartwig.hmftools.linx.chaining.SvChain;
import com.hartwig.hmftools.linx.types.SvCluster;
import com.hartwig.hmftools.linx.types.SvLinkedPair;
import com.hartwig.hmftools.linx.types.SvVarData;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

// post-clustering and chaining routines for annotating clusters, chains and links

public class ClusterAnnotations
{
    private static final Logger LOGGER = LogManager.getLogger(ClusterAnnotations.class);

    public static final String ALL_ANNOTATIONS = "ALL";
    public static final String DOUBLE_MINUTES = "DM";
    public static final String FOLDBACK_MATCHES = "FBM";

    public static boolean runAnnotation(final String annotationsList, final String annotation)
    {
        if(annotationsList.isEmpty())
            return false;

        if(annotationsList.contains(ALL_ANNOTATIONS))
            return true;

        return annotationsList.contains(annotation);
    }

    public static void annotateTemplatedInsertions(final List<SvCluster> clusters, final Map<String, List<SvBreakend>> chrBreakendMap)
    {
        /* for each TI cache if:
            - it has a DB on either or both sides in the same cluster
            - local, the distance to the next SV in the cluster
            - is enclosed by the chain ends or remote if clear
            - how many other TIs it overlaps
         */

        for(final SvCluster cluster : clusters)
        {
            if(cluster.getChains().isEmpty())
                continue;

            isSpecificCluster(cluster);

            // gather up start and end arms from each chain, to determine origin arms for the cluster
            List<String> startEndArms = Lists.newArrayList();

            for(final SvChain chain : cluster.getChains())
            {
                for (int be1 = SE_START; be1 <= SE_END; ++be1)
                {
                    final SvBreakend chainEnd = chain.getOpenBreakend(isStart(be1));
                    if (chainEnd == null)
                        continue;

                    if (!startEndArms.contains(chainEnd.getChrArm()))
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
                        chainEndsCN = chainStartCN;
                }

                for(final SvLinkedPair pair : chain.getLinkedPairs())
                {
                    if(pair.first().type() == SGL || pair.second().type() == SGL)
                        continue;

                    SvBreakend pairLowerBe = pair.getBreakend(true);
                    SvBreakend pairUpperBe = pair.getBreakend(false);
                    double pairCN = (pairLowerBe.copyNumber() + pairUpperBe.copyNumber()) * 0.5;
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

                            if (pairLowerBe.position() >= chainLowerBe.position() - lowerBuffer
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

                    List<SvLinkedPair> uniqueOverlaps = Lists.newArrayList();
                    int overlapCount = 0;

                    for(final SvLinkedPair otherPair : chain.getLinkedPairs())
                    {
                        if(pair == otherPair)
                            continue;

                        if(uniqueOverlaps.stream().anyMatch(x -> x.matches(otherPair)))
                            continue;

                        if(!pair.chromosome().equals(otherPair.chromosome()) || pair.length() <= otherPair.length())
                            continue;

                        uniqueOverlaps.add(otherPair);

                        // make note of this overlap if the overlap distance exceeds the short TI length

                        if((pair.getBreakend(true).position() > otherPair.getBreakend(false).position())
                        || (pair.getBreakend(false).position() < otherPair.getBreakend(true).position()))
                        {
                            continue;
                        }

                        long overlapStart = max(pair.getBreakend(true).position(), otherPair.getBreakend(true).position());
                        long overlapEnd = min(pair.getBreakend(false).position(), otherPair.getBreakend(false).position());

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

                    // closest DB info
                    SvLinkedPair dbFirst = pair.first().getDBLink(pair.firstLinkOnStart());
                    SvLinkedPair dbSecond = pair.second().getDBLink(pair.secondLinkOnStart());

                    pair.setDBLenFirst(dbFirst != null ? dbFirst.length() : NO_DB_MARKER);
                    pair.setDBLenSecond(dbSecond != null ? dbSecond.length() : NO_DB_MARKER);

                    // whether this pair is on the same arm as the chain ends
                    if(startEndArms.contains(pair.firstBreakend().getChrArm()))
                    {
                        pair.setOnArmOfOrigin(true);
                    }
                }
            }
        }
    }

    private static int getTraversedSvCount(final SvCluster cluster, final List<SvBreakend> breakendList, int lowerIndex, int upperIndex)
    {
        String traversedInfo = getTraversedSvData(cluster, breakendList, lowerIndex, upperIndex);

        if(traversedInfo.isEmpty())
            return 0;

        String[] items = traversedInfo.split(";");
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
            LOGGER.error("invalid indices({} & {}) vs breakend list size({})", lowerIndex, upperIndex, breakendList.size());
            return "";
        }

        String traversedInfo = "";

        for (int i = lowerIndex + 1; i <= upperIndex - 1; ++i)
        {
            final SvBreakend breakend = breakendList.get(i);
            final SvCluster otherCluster = breakend.getCluster();

            if (otherCluster == cluster || otherCluster.isResolved())
                continue;

            traversedInfo = appendStr(traversedInfo,
                    String.format("%d %.2f", breakend.orientation(), breakend.getSV().copyNumberChange(breakend.usesStart())),
                    ';');
        }

        return traversedInfo;
    }

    private static int NEXT_SV_DISTANCE = 0;
    private static int NEXT_CLUSTERED_SV_DISTANCE = 1;

    private static int[] getNextClusterSVData(final SvCluster cluster, final List<SvBreakend> breakendList, final SvLinkedPair pair)
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
                int distance = (int)(lowerBreakend.position() - breakend.position());

                if (breakend.getCluster() == cluster)
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

                int distance = (int)(breakend.position() - upperBreakend.position());

                if (breakend.getCluster() == cluster)
                {
                    if (nextSvData[NEXT_CLUSTERED_SV_DISTANCE] == -1 || distance < nextSvData[NEXT_CLUSTERED_SV_DISTANCE])
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

    private static int CHAIN_TI_COUNT = 0;
    private static int CHAIN_TI_DB_COUNT = 1;
    private static int CHAIN_TI_SHORT_COUNT = 2;
    private static int CHAIN_TI_ASMB_COUNT = 3;

    public static void annotateChainedClusters(final SvCluster cluster, long proximityCutoff)
    {
        if(cluster.isResolved() || cluster.getChains().isEmpty())
            return;

        boolean isComplex = cluster.requiresReplication() || !cluster.getFoldbacks().isEmpty();
        boolean isIncomplete = !cluster.isFullyChained(false) || cluster.getTypeCount(SGL) > 0;

        // skip simple chained clusters
        if(cluster.getArmCount() == 1 && cluster.getSvCount() == 2)
            return;

        /* data to gather for each arm in the chain
            - number of links
            - number of short TIs without proximate deletion bridges
            - links with copy number gain
            - start and end locations
         */

        List<String> originArms = Lists.newArrayList();
        List<String> fragmentArms = Lists.newArrayList();

        int chainCount = cluster.getChains().size();
        int unlinkedSvCount = cluster.getUnlinkedSVs().size();
        int inconsistentChains = 0;
        int repeatedChainEndArms = 0;

        // isSpecificCluster(cluster);

        List<String> chainEndArms = Lists.newArrayList();

        for (final SvChain chain : cluster.getChains())
        {
            if (!chain.isConsistent())
            {
                ++inconsistentChains;
                continue;
            }

            final SvBreakend firstBreakend = chain.getOpenBreakend(true);
            final SvBreakend lastBreakend = chain.getOpenBreakend(false);

            final String startChrArm = firstBreakend != null ? firstBreakend.getChrArm() : "";
            final String endChrArm = lastBreakend != null ? lastBreakend.getChrArm() : "";

            if(!startChrArm.isEmpty() && !chainEndArms.contains(startChrArm))
                chainEndArms.add(startChrArm);
            else
                ++repeatedChainEndArms;

            if(!startChrArm.equals(endChrArm))
            {
                if (!endChrArm.isEmpty() && !chainEndArms.contains(endChrArm))
                    chainEndArms.add(startChrArm);
                else
                    ++repeatedChainEndArms;
            }

            Map<String, int[]> armDataMap = new HashMap();

            if(!startChrArm.isEmpty())
                armDataMap.put(startChrArm, new int[CHAIN_TI_ASMB_COUNT+1]);

            if(!endChrArm.isEmpty())
                armDataMap.put(endChrArm, new int[CHAIN_TI_ASMB_COUNT+1]);

            int shortTICount = 0;
            long chainLinkLength = 0;

            for(final SvLinkedPair pair : chain.getLinkedPairs())
            {
                final SvVarData first = pair.first();

                if(pair.first().type() == SGL || pair.second().type() == SGL)
                    continue;

                chainLinkLength += pair.length();

                final String chrArm = first.getBreakend(pair.firstLinkOnStart()).getChrArm();

                int[] armData = armDataMap.get(chrArm);

                if(armData == null)
                {
                    armData = new int[CHAIN_TI_ASMB_COUNT+1];
                    armDataMap.put(chrArm, armData);
                }

                ++armData[CHAIN_TI_COUNT];

                if((pair.getDBLenFirst() > NO_DB_MARKER && pair.getDBLenFirst() <= proximityCutoff)
                || (pair.getDBLenSecond() > NO_DB_MARKER && pair.getDBLenSecond() <= proximityCutoff))
                {
                    ++armData[CHAIN_TI_DB_COUNT];
                }

                if(pair.length() <= SHORT_TI_LENGTH)
                {
                    ++armData[CHAIN_TI_SHORT_COUNT];
                    ++shortTICount;

                    if(pair.isAssembled())
                        ++armData[CHAIN_TI_ASMB_COUNT];
                }
            }
        }

        cluster.setArmData(originArms.size(), fragmentArms.size());

        int armGroupCount = cluster.getArmGroups().size();

        final List<SvVarData> unlinkedRemoteSVs = cluster.getUnlinkedRemoteSVs();

        int inconsistentArmCount = 0;

        for(final SvArmGroup armGroup : cluster.getArmGroups())
        {
            if(!armGroup.isConsistent())
            {
                ++inconsistentArmCount;
                continue;
            }

            for (final SvVarData var : unlinkedRemoteSVs)
            {
                if (armGroup.getSVs().contains(var))
                {
                    ++inconsistentArmCount;
                    continue;
                }
            }
        }

        boolean isComplete = (inconsistentChains == 0) && (repeatedChainEndArms == 0) && (unlinkedSvCount == 0);

        LOGGER.debug("cluster({}) {} chains({} incons={}) chainEnds(arms={} repeats={}) unlinkedSVs({} armCount({} incons={}))",
                cluster.id(), isComplete ? "COMPLETE" : "incomplete",
                chainCount, inconsistentChains, chainEndArms.size(), repeatedChainEndArms,
                unlinkedSvCount, armGroupCount, inconsistentArmCount);

        if(isComplete)
        {
            // chromothripsis is currently defined as fully chained simple cluster
            // but needs to take into account the copy number gain / loss compared with the surrounding chromatid
            if(!isComplex)
            {
                cluster.addAnnotation(CLUSTER_ANNOT_SHATTERING);
            }
        }
    }

    public static void annotateFoldbacks(final List<SvCluster> clusters)
    {
        // now foldbacks are known, add other annotations about them
        // FIXME: the foldback info is not being set on both SVs for chained foldbacks, nor on the correct end
        for(final SvCluster cluster : clusters)
        {
            final Map<String, List<SvBreakend>> chrBreakendMap = cluster.getChrBreakendMap();

            for (final Map.Entry<String, List<SvBreakend>> entry : chrBreakendMap.entrySet())
            {
                int chromosomeFoldbackCount = 0;
                List<SvVarData> chainedFoldbacks = Lists.newArrayList(); // to avoid double counting

                for(final SvBreakend breakend : entry.getValue())
                {
                    if (!breakend.getSV().isFoldback())
                        continue;

                    // only process each SV's foldback once where both breakends are part of it
                    if(!breakend.getSV().isChainedFoldback())
                    {
                        if(!breakend.usesStart())
                            continue;
                    }
                    else
                    {
                        final SvVarData otherSv = breakend.getSV().getChainedFoldbackSv();
                        if(chainedFoldbacks.contains(otherSv))
                            continue;

                        chainedFoldbacks.add(breakend.getSV());
                    }

                    ++chromosomeFoldbackCount;
                }

                if(chromosomeFoldbackCount == 0)
                    continue;

                chainedFoldbacks.clear();
                int foldbackIndex = 0;

                for(final SvBreakend breakend : entry.getValue())
                {
                    // isSpecificSV(breakend.getSV());

                    if(!breakend.getSV().isFoldback())
                        continue;

                    if(!breakend.getSV().isChainedFoldback())
                    {
                        if(!breakend.usesStart())
                            continue;
                    }
                    else
                    {
                        final SvVarData otherSv = breakend.getSV().getChainedFoldbackSv();
                        if(chainedFoldbacks.contains(otherSv))
                            continue;

                        chainedFoldbacks.add(breakend.getSV());
                    }

                    String existingInfo = breakend.getSV().getFoldbackInfo(breakend.usesStart());

                    if(existingInfo.isEmpty())
                        existingInfo = breakend.getSV().getFoldbackInfo(!breakend.usesStart());

                    long armLength = SvUtilities.getChromosomalArmLength(breakend.chromosome(), breakend.arm());
                    double positionPercent;
                    int foldbackRank;

                    if(breakend.arm() == CHROMOSOME_ARM_P)
                    {
                        foldbackRank = foldbackIndex;
                        positionPercent = breakend.position() / (double)armLength;
                    }
                    else
                    {
                        foldbackRank = chromosomeFoldbackCount - foldbackIndex - 1;
                        long chromosomeLength = SvUtilities.CHROMOSOME_LENGTHS.get(breakend.chromosome());
                        long centromere = chromosomeLength - armLength;
                        positionPercent = 1 - (breakend.position() - centromere) / (double)armLength;
                    }

                    ++foldbackIndex;

                    String foldbackInfo = String.format("%s;%s;%d;%.4f",
                            existingInfo, breakend.direction(), foldbackRank, positionPercent);

                    for(int be = SE_START; be <= SE_END; ++be)
                    {
                        boolean isStart = isStart(be);

                        if(breakend.getSV().getFoldbackBreakend(isStart) == null)
                            continue;

                        breakend.getSV().setFoldbackInfo(isStart, foldbackInfo);
                    }
                }
            }
        }
    }

    public static void findIncompleteFoldbackCandidates(final String sampleId, final SvCluster cluster,
            final Map<String, List<SvBreakend>> chrBreakendMap, final CnDataLoader cnAnalyser)
    {
        // for each cluster with incomplete foldbacks, search for candidate clusters which could resolve it

        // for now just focus on single foldback clusters
        if(cluster.getFoldbacks().isEmpty())
            return;

        if(cluster.isResolved()) // eg LowQual
            return;

        int foldbackCount = 0;
        double maxFoldbackPloidy = 0;

        List<SvVarData> clusterFoldbacks = cluster.getFoldbacks();

        Map<SvArmGroup, List<SvVarData>> armGroupFoldbacks = new HashMap();

        for(SvArmGroup armGroup : cluster.getArmGroups())
        {
            final String chromosome = armGroup.chromosome();
            final String arm = armGroup.arm();

            List<SvVarData> chainedFoldbacks = Lists.newArrayList(); // to avoid double counting

            List<SvVarData> foldbacks = Lists.newArrayList();

            // first divvy up foldbacks into their chromosomal arms
            for (SvVarData var : clusterFoldbacks)
            {
                SvBreakend breakend = null;

                if (!var.isChainedFoldback())
                {
                    breakend = var.getBreakend(true);
                }
                else
                {
                    // only process one of the 2 SVs involed in a chained foldback
                    final SvVarData otherSv = var.getChainedFoldbackSv();
                    if (chainedFoldbacks.contains(otherSv))
                        continue;

                    chainedFoldbacks.add(var);
                    breakend = var.getChainedFoldbackBreakend();
                }

                // skip this foldback if it's not on the current arm group
                if (!breakend.chromosome().equals(chromosome) || !breakend.arm().equals(arm))
                {
                    continue;
                }

                if (foldbacks.isEmpty())
                {
                    armGroupFoldbacks.put(armGroup, foldbacks);
                }

                // build up an ordered list of the foldbacks, telomere to centromere
                int index = 0;
                while (index < foldbacks.size())
                {
                    SvVarData otherFoldback = foldbacks.get(index);

                    SvBreakend otherBreakend = !otherFoldback.isChainedFoldback()
                            ? otherFoldback.getBreakend(true) : otherFoldback.getChainedFoldbackBreakend();

                    if (arm == CHROMOSOME_ARM_P && breakend.position() < otherBreakend.position())
                        break;
                    else if (arm == CHROMOSOME_ARM_Q && breakend.position() > otherBreakend.position())
                        break;

                    ++index;
                }

                foldbacks.add(index, var);

                maxFoldbackPloidy = max(maxFoldbackPloidy, breakend.getSV().ploidyMin());
                ++foldbackCount;
            }
        }

        boolean isMultiArm = armGroupFoldbacks.size() > 1;

        final Map<String, List<SvCNData>> chrCopyNumberDataMap = cnAnalyser.getChrCnDataMap();

        // now look at each arm with foldbacks independently
        for(Map.Entry<SvArmGroup, List<SvVarData>> entry : armGroupFoldbacks.entrySet())
        {
            final SvArmGroup armGroup = entry.getKey();
            List<SvVarData> foldbacks = entry.getValue();

            final String chromosome = armGroup.chromosome();
            final String arm = armGroup.arm();

            // find the outermost foldbacks to check against all other opposing cluster breakends
            SvBreakend breakendLower = null; // lowest facing up
            SvBreakend breakendUpper = null; // highest facing down

            String foldbackIds = "";
            String foldbackPloidies = "";
            String foldbackOrientations = "";
            String foldbackChained = "";

            for (SvVarData var : foldbacks)
            {
                SvBreakend breakend = !var.isChainedFoldback() ? var.getBreakend(true) : var.getChainedFoldbackBreakend();

                if (breakend.orientation() == 1 && (breakendUpper == null || breakend.position() > breakendUpper.position()))
                {
                    breakendUpper = breakend;
                }
                else if (breakend.orientation() == -1 && (breakendLower == null || breakend.position() < breakendLower.position()))
                {
                    breakendLower = breakend;
                }

                String direction = breakend.direction();
                foldbackOrientations = appendStr(foldbackOrientations, direction, ';');
                double ploidy = direction == DIRECTION_CENTROMERE ? abs(var.ploidyMin()) : -abs(var.ploidyMin());
                foldbackPloidies = appendStr(foldbackPloidies, String.format("%.2f", ploidy), ';');
                foldbackIds = appendStr(foldbackIds, breakend.getSV().id(), ';');
                foldbackChained = appendStr(foldbackChained, Boolean.toString(var.isChainedFoldback()), ';');
            }

            List<SvCNData> cnDataList = chrCopyNumberDataMap.get(chromosome);

            long telomereEndPos = 0;
            long centromereEndPos = 0;
            double telomereEndCN = 0;
            double telomereEndMap = 0;
            double centromereEndCN = 0;
            double centromereEndMap = 0;
            double telomereMinFacingPloidy = Double.NaN;
            double centromereMinFacingPloidy = Double.NaN;

            List<SvBreakend> clusterBreakendList = cluster.getChrBreakendMap().get(chromosome);

            for(int i = 0; i < clusterBreakendList.size(); ++i)
            {
                SvBreakend breakend = clusterBreakendList.get(i);
                SvBreakend nextBreakend = i < clusterBreakendList.size() - 1 ? clusterBreakendList.get(i+1) : null;

                if(breakend.arm() != arm)
                {
                    if(arm == CHROMOSOME_ARM_P)
                        break;
                    else
                        continue;
                }

                // take CN data from the breakends closest to the telomere and centromere
                if((arm == CHROMOSOME_ARM_P && telomereEndPos == 0)
                || (arm == CHROMOSOME_ARM_Q && nextBreakend == null))
                {
                    telomereEndPos = breakend.position();

                    SvCNData cnData = breakend.getSV().getCopyNumberData(breakend.usesStart(), arm == CHROMOSOME_ARM_P);

                    if(cnData != null)
                    {
                        telomereEndCN = cnData.CopyNumber;
                        telomereEndMap = cnData.majorAllelePloidy();
                    }
                }

                if((arm == CHROMOSOME_ARM_P && (nextBreakend == null || nextBreakend.arm() != arm))
                || (arm == CHROMOSOME_ARM_Q && centromereEndPos == 0))
                {
                    centromereEndPos = breakend.position();

                    SvCNData cnData = breakend.getSV().getCopyNumberData(breakend.usesStart(), arm == CHROMOSOME_ARM_Q);

                    if(cnData != null)
                    {
                        centromereEndCN = cnData.CopyNumber;
                        centromereEndMap = cnData.majorAllelePloidy();
                    }
                }
            }

            double[] netPloidies = calcNetCopyNumberChangeAcrossCluster(clusterBreakendList, arm, false);

            telomereMinFacingPloidy = netPloidies[0];
            centromereMinFacingPloidy = netPloidies[1];

            // get centromere & telomere data
            double[] centromereCNData = cnAnalyser.getCentromereCopyNumberData(chromosome, arm.equals(CHROMOSOME_ARM_P));
            double telomereCN = 0;
            double telomereMAP = 0;

            if(cnDataList != null && !cnDataList.isEmpty())
            {
                SvCNData telemoreData = arm == CHROMOSOME_ARM_P ? cnDataList.get(0) : cnDataList.get(cnDataList.size() - 1);
                telomereCN = telemoreData.CopyNumber;
                telomereMAP = telemoreData.majorAllelePloidy();
            }

            // now find all clusters with opposing breakends
            final List<SvBreakend> allBreakendList = chrBreakendMap.get(chromosome);

            // cache max opposing and net ploidy for each opposing cluster
            List<SvCluster> processedClusters = Lists.newArrayList();
            int opposingClusterCount = 0;
            String allClusterInfo = "";

            for (int i = 0; i <= 1; ++i)
            {
                SvBreakend fbBreakend = (i == 0) ? breakendLower : breakendUpper;

                if (fbBreakend == null)
                    continue;

                int index = fbBreakend.getChrPosIndex();

                while (true)
                {
                    if (fbBreakend.orientation() == 1) // walk in the direction the foldback faces
                        --index;
                    else
                        ++index;

                    if (index < 0 || index >= allBreakendList.size())
                        break;

                    final SvBreakend nextBreakend = allBreakendList.get(index);
                    if (nextBreakend.arm() != fbBreakend.arm())
                        break;

                    if (nextBreakend.orientation() == fbBreakend.orientation())
                        continue;

                    final SvCluster nextCluster = nextBreakend.getCluster();

                    if (nextCluster == cluster || processedClusters.contains(nextCluster))
                        continue;

                    processedClusters.add(nextCluster);

                    if (nextCluster.isResolved())
                        continue;

                    // found an opposing non-simple cluster, gather up details about it
                    // max opposing CN min poidy and net ploidy
                    final List<SvBreakend> nextClusterBreakends = nextCluster.getChrBreakendMap().get(chromosome);

                    double maxOpposingPloidy = 0;
                    double netPloidy = 0;

                    for (final SvBreakend otherBreakend : nextClusterBreakends)
                    {
                        if (otherBreakend.arm() != fbBreakend.arm())
                            continue;

                        netPloidy += otherBreakend.getSV().ploidyMin() * otherBreakend.orientation();

                        // take the max opposing breakend for lower and upper FB breakends
                        if (breakendLower != null && otherBreakend.orientation() == 1
                                && otherBreakend.position() > breakendLower.position())
                        {
                            maxOpposingPloidy = max(maxOpposingPloidy, otherBreakend.getSV().ploidyMin());
                        }
                        else if (breakendUpper != null && otherBreakend.orientation() == -1
                                && otherBreakend.position() < breakendUpper.position())
                        {
                            maxOpposingPloidy = max(maxOpposingPloidy, otherBreakend.getSV().ploidyMin());
                        }
                    }

                    ++opposingClusterCount;

                    String clusterInfo = String.format("%d/%d/%.2f/%.2f",
                            nextCluster.id(), nextCluster.getSvCount(), netPloidy, maxOpposingPloidy);

                    allClusterInfo = appendStr(allClusterInfo, clusterInfo, ';');
                }
            }

            // SampleId,ClusterId,ClusterCount,ClusterDesc,FoldbackCount,
            String infoStr = String.format("%s,%d,%d,%s,%d",
                    sampleId, cluster.id(), cluster.getSvCount(), cluster.getDesc(), foldbackCount);

            // IsMultiArm,Chromosome,Arm,ArmSvCount,ArmFoldbackCount,FbIds,FbOrientations,FbPloidies,FbChainedTypes,
            infoStr += String.format(",%s,%s,%s,%d,%d,%s,%s,%s,%s",
                    isMultiArm, chromosome, arm, armGroup.getSVs().size(), foldbacks.size(),
                    foldbackIds, foldbackOrientations, foldbackPloidies, foldbackChained);

            // arm cluster SV data:
            // TeloEndPos,CentroEndPos,TeloEndCN,CentroEndCN,TeloEndMAP,CentroEndMAP,TeloMinFacingPloidy,CentroMinFacingPloidy
            infoStr += String.format(",%d,%d,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f",
                    telomereEndPos, centromereEndPos, telomereEndCN, centromereEndCN, telomereEndMap, centromereEndMap,
                    telomereMinFacingPloidy, centromereMinFacingPloidy);

            // TeleCN,TeloMAP,PreCentroCN,PreCentroMAP,PostCentroCN,PostCentroMAP
            infoStr += String.format(",%.2f,%.2f,%.2f,%.2f,%.2f,%.2f",
                    telomereCN, telomereMAP, centromereCNData[CN_SEG_DATA_CN_BEFORE], centromereCNData[CN_SEG_DATA_MAP_BEFORE],
                    centromereCNData[CN_SEG_DATA_CN_AFTER], centromereCNData[CN_SEG_DATA_MAP_AFTER]);

            // MaxFoldbackPloidy,OpposingClusterCount,OpposingClusterInfo
            infoStr += String.format(",%.2f,%d,%s",
                    maxFoldbackPloidy, opposingClusterCount, allClusterInfo);

            LOGGER.info("INCONSIST_FBS: {}", infoStr);
        }
    }

    public static double[] calcNetCopyNumberChangeAcrossCluster(final List<SvBreakend> breakendList, final String arm,
            boolean breakOnOtherCluster)
    {
        // find the net copy number change facing both the telomere and centromere
        double telomereMinFacingPloidy = 0;
        double centromereMinFacingPloidy = 0;

        for(int i = 0; i <= 1; ++i) // first forwards, then backwards through the breakends
        {
            byte currentOrientation = (i == 0) ? (byte)1 : (byte)-1;
            boolean traverseUp = (i == 0);

            long prevPos = 0;
            boolean seenNegOrient = false;
            boolean seenPosOrient = false;
            double netPloidy = 0;
            double maxPloidy = 0;

            int index = traverseUp ? 0 : breakendList.size() - 1;
            int prevFullListIndex = -1;
            while(index >= 0 && index < breakendList.size())
            {
                SvBreakend breakend = breakendList.get(index);

                // optionally break as soon as this cluster skips over another cluster's breakend
                if(breakOnOtherCluster && prevFullListIndex > -1 && abs(breakend.getChrPosIndex() - prevFullListIndex) > 1)
                    break;

                if (breakend.arm() == arm)
                {
                    // TEMP: skip multiple breakends
                    boolean skipBreakend = false;
                    if (prevPos == breakend.position())
                    {
                        if ((breakend.orientation() == 1 && seenPosOrient)
                        || (breakend.orientation() == -1 && seenNegOrient))
                        {
                            LOGGER.debug("skipping multiple breakend({}) with same orientation", breakend.toString());
                            skipBreakend = true;
                        }
                        else
                        {
                            if (breakend.orientation() == 1)
                                seenPosOrient = true;
                            else
                                seenNegOrient = true;
                        }
                    }
                    else
                    {
                        prevPos = breakend.position();
                        seenPosOrient = false;
                        seenNegOrient = false;
                    }

                    if (!skipBreakend)
                    {
                        if (breakend.orientation() == currentOrientation)
                            netPloidy += breakend.copyNumberChange();
                        else
                            netPloidy -= breakend.copyNumberChange();

                        // find the maximum net copy number change
                        maxPloidy = max(maxPloidy, netPloidy);
                    }
                }
                else
                {
                    // early exit if into next arm
                    if(traverseUp && arm == CHROMOSOME_ARM_P)
                        break;
                    else if(!traverseUp && arm == CHROMOSOME_ARM_Q)
                        break;
                }

                if(traverseUp)
                    ++index;
                else
                    --index;
            }

            boolean facesTelomere = (currentOrientation == 1) == (arm == CHROMOSOME_ARM_P);

            if(facesTelomere)
            {
                telomereMinFacingPloidy = maxPloidy;
            }
            else
            {
                centromereMinFacingPloidy = maxPloidy;
            }
        }

        double[] results = {telomereMinFacingPloidy, centromereMinFacingPloidy};
        return results;
    }


}
