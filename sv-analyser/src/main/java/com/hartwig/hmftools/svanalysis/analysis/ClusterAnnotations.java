package com.hartwig.hmftools.svanalysis.analysis;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.round;

import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.BND;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.SGL;
import static com.hartwig.hmftools.svanalysis.analysis.ClusterAnalyser.SHORT_TI_LENGTH;
import static com.hartwig.hmftools.svanalysis.analysis.LinkFinder.NO_DB_MARKER;
import static com.hartwig.hmftools.svanalysis.analysis.SvUtilities.copyNumbersEqual;
import static com.hartwig.hmftools.svanalysis.types.SvChain.getRepeatedSvSequence;
import static com.hartwig.hmftools.svanalysis.types.SvCluster.CLUSTER_ANNONTATION_CT;
import static com.hartwig.hmftools.svanalysis.types.SvCluster.CLUSTER_ANNONTATION_DM;
import static com.hartwig.hmftools.svanalysis.types.SvCluster.RESOLVED_TYPE_DEL_EXT_TI;
import static com.hartwig.hmftools.svanalysis.types.SvCluster.RESOLVED_TYPE_DUP_EXT_TI;
import static com.hartwig.hmftools.svanalysis.types.SvLinkedPair.LINK_TYPE_TI;
import static com.hartwig.hmftools.svanalysis.types.SvVarData.SVI_END;
import static com.hartwig.hmftools.svanalysis.types.SvVarData.SVI_START;
import static com.hartwig.hmftools.svanalysis.types.SvVarData.isStart;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.svanalysis.types.SvArmCluster;
import com.hartwig.hmftools.svanalysis.types.SvArmGroup;
import com.hartwig.hmftools.svanalysis.types.SvBreakend;
import com.hartwig.hmftools.svanalysis.types.SvChain;
import com.hartwig.hmftools.svanalysis.types.SvCluster;
import com.hartwig.hmftools.svanalysis.types.SvLinkedPair;
import com.hartwig.hmftools.svanalysis.types.SvVarData;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

// post-clustering and chaining routines for annotating clusters, chains and links

public class ClusterAnnotations
{
    private static final Logger LOGGER = LogManager.getLogger(ClusterAnnotations.class);


    public static void annotateClusterArmSegments(final SvCluster cluster)
    {
        if (cluster.isResolved())
            return;

        if (cluster.getCount() < 6)
            return;

        // isSpecificCluster(cluster);

        final Map<String, List<SvBreakend>> chrBreakendMap = cluster.getChrBreakendMap();

        for (final Map.Entry<String, List<SvBreakend>> entry : chrBreakendMap.entrySet())
        {
            List<SvBreakend> breakendList = entry.getValue();

            if (breakendList.size() < 6)
                continue;

            // first establish the lowest copy number segments
            int bndCount = 0;
            boolean hasConsecutiveOrientations = false;
            double lowestCopyNumber = -1;

            List<SvArmCluster> localSegments = Lists.newArrayList();
            SvArmCluster armCuster = null;

            for (int i = 0; i < breakendList.size() - 1; ++i)
            {
                final SvBreakend breakend = breakendList.get(i);

                if(!hasConsecutiveOrientations && i < breakendList.size() - 1)
                {
                    hasConsecutiveOrientations = (breakend.orientation() != breakendList.get(i+1).orientation());
                }

                double lowCopyNumber = breakend.getCopyNumber(false);

                if (breakend.orientation() == 1 && (lowestCopyNumber == -1 || lowCopyNumber < lowestCopyNumber))
                {
                    // a new lowest segment, which will be the first original DSB
                    // all previous breakends will be put into the first arm cluster
                    lowestCopyNumber = lowCopyNumber;

                    localSegments.clear();

                    // add all previous breakends to the first arm cluster
                    armCuster = new SvArmCluster(localSegments.size(), cluster, breakend.chromosome(), breakend.arm());

                    for (int j = 0; j <= i; ++j)
                        armCuster.addBreakend(breakendList.get(j));

                    localSegments.add(armCuster);
                    armCuster = null;
                }
                else if (breakend.orientation() == 1 && copyNumbersEqual(lowestCopyNumber, lowCopyNumber))
                {
                    // end of a segment since next CN equals the low - add the last breakend
                    if(armCuster == null)
                    {
                        // probably an indication of failed clustering, so work around it
                        armCuster = new SvArmCluster(localSegments.size(), cluster, breakend.chromosome(), breakend.arm());
                        localSegments.add(armCuster);
                    }

                    armCuster.addBreakend(breakend);
                    armCuster = null;
                }
                else
                {
                    // either continuation of an existing segment or start of a new one
                    if (armCuster == null)
                    {
                        armCuster = new SvArmCluster(localSegments.size(), cluster, breakend.chromosome(), breakend.arm());
                        localSegments.add(armCuster);
                    }

                    armCuster.addBreakend(breakend);
                }

                if (breakend.getSV().type() == BND)
                    ++bndCount;
            }

            if (bndCount > 0 || !hasConsecutiveOrientations)
                continue;

            for (final SvArmCluster armCluster : localSegments)
            {
                List<SvBreakend> acBreakendList = armCluster.getBreakends();

                if(acBreakendList.size() < 4)
                    continue;

                int backwardCount = (int)acBreakendList.stream().filter(x -> x.orientation() == 1).count();
                int forwardCount = (int)acBreakendList.stream().filter(x -> x.orientation() == -1).count();

                // form into mutually exclusive pairs based on copy number change match
                List<SvBreakend> unlinkedBreakends = Lists.newArrayList(acBreakendList);
                List<SvLinkedPair> pairs = formPossibleLinkedPairsByShortest(unlinkedBreakends);

                LOGGER.debug("cluster({}) armCluster({} : {}_{}) count({} fwd={} bak={}) pairs({}) unlinked({})",
                        cluster.id(), armCluster.id(), armCluster.chromosome(), armCluster.arm(),
                        acBreakendList.size(), forwardCount, backwardCount, pairs.size(), unlinkedBreakends.size());
            }
        }
    }

    private static List<SvLinkedPair> formPossibleLinkedPairsByShortest(List<SvBreakend> unlinkedBreakends)
    {
        List<SvLinkedPair> pairs = Lists.newArrayList();

        while(!unlinkedBreakends.isEmpty())
        {
            boolean pairFound = false;

            for(int i = 0; i < unlinkedBreakends.size() - 2; ++i)
            {
                final SvBreakend be1 = unlinkedBreakends.get(i);
                final SvBreakend be2 = unlinkedBreakends.get(i+1);

                if(be1.orientation() != -1 || be2.orientation() != 1)
                    continue;

                // pair off these if the CN matches
                if (!copyNumbersEqual(be1.copyNumberChange(), be2.copyNumberChange()))
                    continue;

                SvLinkedPair pair = SvLinkedPair.from(be1, be2, LINK_TYPE_TI);
                pairs.add(pair);

                pairFound = true;
                unlinkedBreakends.remove(i+1); // higher index removed first
                unlinkedBreakends.remove(i);
                break;
            }

            if(!pairFound)
                break;
        }

        return pairs;
    }

    private static List<SvLinkedPair> formPossibleLinkedPairsConsecutively(List<SvBreakend> unlinkedBreakends)
    {
        List<SvLinkedPair> pairs = Lists.newArrayList();

        int j = 0;
        while(j < unlinkedBreakends.size())
        {
            final SvBreakend be1 = unlinkedBreakends.get(j);

            if(be1.orientation() == 1)
            {
                ++j;
                continue;
            }

            boolean linkFound = false;
            for(int k = j+1; k < unlinkedBreakends.size(); ++k)
            {
                final SvBreakend be2 = unlinkedBreakends.get(k);

                if (be2.orientation() == -1)
                    continue;

                // pair off these if the CN matches
                if (!copyNumbersEqual(be1.copyNumberChange(), be2.copyNumberChange()))
                    continue;

                SvLinkedPair pair = SvLinkedPair.from(be1, be2, LINK_TYPE_TI);
                pairs.add(pair);

                linkFound = true;
                unlinkedBreakends.remove(k); // higher index removed first
                unlinkedBreakends.remove(j);
                break;

            }

            if(linkFound)
                continue;
            else
                ++j;
        }

        return pairs;
    }

    public static void annotateTemplatedInsertions(final List<SvCluster> clusters, final Map<String, List<SvBreakend>> chrBreakendMap)
    {
        /* work out:
            - whether a TI has a DB on either or both sides in the same cluster
            - number of chained assembled TIs in a row
            - if local, the distance to the next SV in the cluster
         */

        for(final SvCluster cluster : clusters)
        {
            if(cluster.getChains().isEmpty())
                continue;

            // isSpecificCluster(cluster);

            // gather up start and end arms from each chain
            List<String> startEndArms = Lists.newArrayList();

            for(final SvChain chain : cluster.getChains())
            {
                for (int be1 = SVI_START; be1 <= SVI_END; ++be1)
                {
                    boolean isFirst = isStart(be1);

                    final SvBreakend chainEnd = chain.getChainEndSV(isFirst).getBreakend(chain.chainEndOpenOnStart(isFirst));
                    if (chainEnd == null)
                        continue;

                    if (!startEndArms.contains(chainEnd.getChrArm()))
                        startEndArms.add(chainEnd.getChrArm());
                }
            }

            for(final SvChain chain : cluster.getChains())
            {
                for(final SvLinkedPair pair : chain.getLinkedPairs())
                {
                    if(pair.linkType() != LINK_TYPE_TI)
                        continue;

                    if(pair.first().type() == SGL || pair.second().type() == SGL)
                        continue;

                    // find closest SV in this cluster
                    final SvVarData first = pair.first();
                    final SvVarData second = pair.second();
                    SvBreakend firstBreakend = first.getBreakend(pair.firstLinkOnStart());
                    final List<SvBreakend> breakendList = chrBreakendMap.get(firstBreakend.chromosome());

                    int[] nextSVData = getNextClusterSVData(cluster, breakendList, pair);
                    pair.setNextSVData(nextSVData[NEXT_SV_DISTANCE], nextSVData[NEXT_SV_TRAVERSED_COUNT]);

                    pair.setTraversedSVCount(getTraversedSvCount(cluster, breakendList,
                            pair.getBreakend(true).getChrPosIndex(), pair.getBreakend(false).getChrPosIndex()));

                    SvLinkedPair dbFirst = first.getDBLink(pair.firstLinkOnStart());
                    SvLinkedPair dbSecond = pair.second().getDBLink(pair.secondLinkOnStart());

                    pair.setDBLenFirst(dbFirst != null ? dbFirst.length() : NO_DB_MARKER);
                    pair.setDBLenSecond(dbSecond != null ? dbSecond.length() : NO_DB_MARKER);

                    if(startEndArms.contains(firstBreakend.getChrArm()))
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
            return "";

        String traversedInfo = "";

        for (int i = lowerIndex + 1; i <= upperIndex - 1; ++i)
        {
            final SvBreakend breakend = breakendList.get(i);
            final SvCluster otherCluster = breakend.getSV().getCluster();

            if (otherCluster == cluster || otherCluster.isResolved())
                continue;

            if(!traversedInfo.isEmpty())
                traversedInfo += ";";

            traversedInfo += String.format("%d %.2f", breakend.orientation(), breakend.getSV().copyNumberChange(breakend.usesStart()));
        }

        return traversedInfo;
    }

    private static int NEXT_SV_DISTANCE = 0;
    private static int NEXT_SV_TRAVERSED_COUNT = 1;

    private static int[] getNextClusterSVData(final SvCluster cluster, final List<SvBreakend> breakendList, final SvLinkedPair pair)
    {
        // walk forward and backwards from this pair to the closest SV in the same cluster
        // counting the number of non-trivial SVs traversed in the process
        int[] nextSvData = {-1, 0};

        SvBreakend firstBreakend = pair.first().getBreakend(pair.firstLinkOnStart());
        SvBreakend secondBreakend = pair.second().getBreakend(pair.secondLinkOnStart());
        int lowerIndex = min(firstBreakend.getChrPosIndex(), secondBreakend.getChrPosIndex());
        int upperIndex = max(firstBreakend.getChrPosIndex(), secondBreakend.getChrPosIndex());
        int svsTraversed = 0;

        if(lowerIndex > 0)
        {
            for(int i = lowerIndex - 1; i >= 0; --i)
            {
                final SvBreakend breakend = breakendList.get(i);

                if (breakend.getSV().getCluster() == cluster)
                {
                    final SvBreakend refBreakend = breakendList.get(lowerIndex);
                    nextSvData[NEXT_SV_DISTANCE] = (int)(refBreakend.position() - breakend.position());
                    nextSvData[NEXT_SV_TRAVERSED_COUNT] = svsTraversed;
                    break;
                }
                else if(!breakend.getSV().getCluster().isResolved())
                {
                    ++svsTraversed;
                }
            }
        }

        if(upperIndex < breakendList.size() - 1)
        {
            svsTraversed = 0;

            for(int i = upperIndex + 1; i < breakendList.size(); ++i)
            {
                final SvBreakend breakend = breakendList.get(i);

                if (breakend.getSV().getCluster() == cluster)
                {
                    final SvBreakend refBreakend = breakendList.get(upperIndex);
                    long distance = breakend.position() - refBreakend.position();

                    if (nextSvData[NEXT_SV_DISTANCE] == -1 || distance < nextSvData[NEXT_SV_DISTANCE])
                    {
                        nextSvData[NEXT_SV_TRAVERSED_COUNT] = svsTraversed;
                        nextSvData[NEXT_SV_DISTANCE] = (int)distance;
                    }

                    break;
                }
                else if(!breakend.getSV().getCluster().isResolved())
                {
                    ++svsTraversed;
                }
            }
        }

        return nextSvData;
    }

    public static void analyseOverlappingTIs(final String sampleId, final Map<String, List<SvBreakend>> chrBreakendMap)
    {
        // look for any TIs where one falls within the bounds of the other, as a sign of replication after shattering
        for(final Map.Entry<String, List<SvBreakend>> entry : chrBreakendMap.entrySet())
        {
            List<SvBreakend> breakendList = entry.getValue();

            for(int i = 1; i < breakendList.size() - 3; ++i)
            {
                final SvBreakend beStart1 = breakendList.get(i);
                final SvBreakend beStart2 = breakendList.get(i + 1);

                if(beStart1.orientation() == 1 || beStart2.orientation() == 1)
                    continue;

                final SvBreakend beEnd1 = breakendList.get(i + 2);
                final SvBreakend beEnd2 = breakendList.get(i + 3);

                if(beEnd1.orientation() == -1 || beEnd2.orientation() == -1)
                    continue;

                // now the orientations are correct, check if these 2 sets of breakends form one TI enclosed by the other

                SvLinkedPair pair1 = beStart1.getSV().getLinkedPair(beStart1.usesStart());
                SvLinkedPair pair2 = beStart2.getSV().getLinkedPair(beStart2.usesStart());

                boolean pair1Matched = pair1 != null && pair1.hasBreakend(beEnd2);
                boolean pair2Matched = pair2 != null && pair2.hasBreakend(beEnd1);

                boolean sameCluster = beStart1.getSV().getCluster() == beStart2.getSV().getCluster()
                        && beStart1.getSV().getCluster() == beEnd1.getSV().getCluster()
                        && beStart1.getSV().getCluster() == beEnd2.getSV().getCluster();

                long startGap = beStart2.position() - beStart1.position();

                if((pair1Matched && pair2Matched) || (pair1 != null && pair2 != null))
                {
                    LOGGER.info("sample({}) cluster({}) matched TIs: outer({}) inner({}) gap(start={} middle={} end={})",
                            sampleId, sameCluster ? beStart1.getSV().getCluster().id() : "multiple",
                            pair1.toString(), pair2.toString(), beStart2.position() - beStart1.position(),
                            beEnd2.position() - beStart2.position(), beEnd2.position() - beEnd1.position());

                    long endGap = pair1.length() - pair2.length() - startGap;

                    if(!pair1Matched || !pair2Matched)
                    {
                        if (endGap < 0)
                            continue;
                    }

                    LOGGER.info("sample({}) cluster({}) enclosed TIs: outer({}) inner({}) gap(start={} middle={} end={})",
                            sampleId, sameCluster ? beStart1.getSV().getCluster().id() : "multiple",
                            pair1.toString(), pair2.toString(), startGap, pair2.length(), endGap);

                    LOGGER.info("ENCLOSED_TI: {},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{}",
                            sampleId, sameCluster ? beStart1.getSV().getCluster().id() : "multiple",
                            pair1.first().chromosome(pair1.firstLinkOnStart()), pair1Matched, pair1Matched,
                            pair1.first().id(), pair1.first().position(pair1.firstLinkOnStart()),
                            pair1.second().id(), pair1.second().position(pair1.secondLinkOnStart()),
                            pair2.first().id(), pair2.first().position(pair2.firstLinkOnStart()),
                            pair2.second().id(), pair2.second().position(pair2.secondLinkOnStart()),
                            pair1.length(), pair2.length(), startGap, endGap);

                }
                else
                {
                    // breakend form 2 TIs with one enclosed by pairings don't match for some reason
                    LOGGER.info("ENCLOSED_TI: {},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{}",
                            sampleId, sameCluster ? beStart1.getSV().getCluster().id() : "multiple",
                            beStart1.chromosome(), pair1Matched, pair1Matched,
                            beStart1.getSV().id(), beStart1.position(), beStart2.getSV().id(), beStart2.position(),
                            beEnd1.getSV().id(), beEnd1.position(), beEnd2.getSV().id(), beEnd2.position(),
                            beEnd2.position() - beStart1.position(), beEnd1.position() - beStart2.position(),
                            startGap, beEnd2.position() - beEnd1.position());
                }
            }
        }
    }

    private static int CHAIN_TI_COUNT = 0;
    private static int CHAIN_TI_DB_COUNT = 1;
    private static int CHAIN_TI_SHORT_COUNT = 2;
    private static int CHAIN_TI_ASMB_COUNT = 3;

    public static void classifyChainedClusters(final SvCluster cluster, long proximityCutoff)
    {
        if(cluster.isResolved())
            return;

        boolean isComplex = cluster.hasReplicatedSVs() || !cluster.getFoldbacks().isEmpty();
        boolean isIncomplete = !cluster.isFullyChained() || cluster.getTypeCount(SGL) > 0;

        // skip simple chained clusters
        if(cluster.getArmCount() == 1 && cluster.getCount() == 2)
            return;

        // isSpecificCluster(cluster);

        /* data to gather for each arm in the chain
            - number of links
            - number of short TIs without proximate deletion bridges
            - links with copy number gain
            - start and end locations
         */

        boolean allChainsConsistent = true;
        List<String> originArms = Lists.newArrayList();
        List<String> fragmentArms = Lists.newArrayList();

        for (final SvChain chain : cluster.getChains())
        {
            if(!chain.isConsistent())
            {
                allChainsConsistent = false;
                continue;
            }

            Map<String, int[]> armDataMap = new HashMap();

            final SvBreakend firstBreakend = chain.getFirstSV().getBreakend(chain.firstLinkOpenOnStart());
            final SvBreakend lastBreakend = chain.getLastSV().getBreakend(chain.lastLinkOpenOnStart());

            final String startChrArm = firstBreakend != null ? firstBreakend.getChrArm() : "";
            final String endChrArm = lastBreakend != null ? lastBreakend.getChrArm() : "";

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

            // check for synthetic DELs and DUPs from longer chains
            if(allChainsConsistent && !isIncomplete && !isComplex
                    && cluster.getChains().size() == 1 && chain.getLinkCount() == shortTICount
                    && firstBreakend.getChrArm().equals(lastBreakend.getChrArm())
                    && firstBreakend.orientation() != lastBreakend.orientation())
            {
                long syntheticLength = abs(lastBreakend.position() - firstBreakend.position());
                long avgLinkLength = round(chainLinkLength/chain.getLinkCount());

                cluster.setSynDelDupData(syntheticLength, avgLinkLength);

                if((firstBreakend.position() < lastBreakend.position() && firstBreakend.orientation() == 1)
                        || (lastBreakend.position() < firstBreakend.position() && lastBreakend.orientation() == 1))
                {
                    cluster.setResolved(false, RESOLVED_TYPE_DEL_EXT_TI);
                }
                else
                {
                    cluster.setResolved(false, RESOLVED_TYPE_DUP_EXT_TI);
                }

                LOGGER.debug("cluster({}) chainLinks({}) synLen({}) avgTILen({}) marked as {}",
                        cluster.id(), chain.getLinkCount(), syntheticLength, avgLinkLength, cluster.getResolvedType());
            }

            String chainInfo = startChrArm + "-" + endChrArm;

            for (Map.Entry<String,int[]> entry : armDataMap.entrySet())
            {
                final String chrArm = entry.getKey();
                final int[] armData = entry.getValue();

                int linkCount = armData[CHAIN_TI_COUNT];
                int dbCount = armData[CHAIN_TI_DB_COUNT];
                int shortCount = armData[CHAIN_TI_SHORT_COUNT];
                int assembledCount = armData[CHAIN_TI_ASMB_COUNT];

                boolean isOrigin = false;
                if(chrArm.equals(startChrArm) || chrArm.equals(endChrArm))
                {
                    isOrigin = true;
                }
                else if(shortCount > dbCount)
                {
                    isOrigin = false;
                }

                if(isOrigin)
                {
                    // if this arm exists twice already, then more than 2 chains end up on the same arm which is invalid
                    long chrArmCount = originArms.stream().filter(x -> x.equals(chrArm)).count();

                    if(chrArmCount >= 2)
                        allChainsConsistent = false;

                    originArms.add(chrArm);

                    if(fragmentArms.contains(chrArm))
                        fragmentArms.remove(chrArm);
                }
                else if(!isOrigin && !fragmentArms.contains(chrArm))
                {
                    fragmentArms.add(chrArm);
                }

                chainInfo += String.format(" %s %s: LK=%d DB=%d SH=%d AS=%d",
                        isOrigin ? "O" : "F", chrArm, linkCount, dbCount, shortCount, assembledCount);
            }

            chain.setDetails(chainInfo);

            // LOGGER.debug("cluster({}) chain({}) {}", cluster.id(), chain.id(), chainInfo);
        }

        cluster.setArmData(originArms.size(), fragmentArms.size());

        if(allChainsConsistent && !isIncomplete && !isComplex)
        {
            cluster.addAnnotation(String.format("%s", CLUSTER_ANNONTATION_CT));
        }
    }

    public static void reportClusterNeoChromosomes(final SvCluster cluster)
    {
        if(cluster.getChains().isEmpty())
            return;

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
        }

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
                cluster.getChains().size(), inconsistentChains, chainEndArms.size(), repeatedChainEndArms,
                unlinkedSvCount, armGroupCount, inconsistentArmCount);

        if(isComplete)
        {
            cluster.addAnnotation(String.format("COMPLETE"));
        }
    }

    private static double DOUBLE_MINUTE_PLOIDY_THRESHOLD = 8;
    private static double DOUBLE_MINUTE_PLOIDY_GAP_RATIO = 3;

    public static void reportDoubleMinutes(final SvCluster cluster, final Map<String, List<SvBreakend>> chrBreakendMap)
    {
        // order SVs in descending ploidy order
        List<Double> ploidyList = Lists.newArrayList();
        List<SvVarData> indexSvList = Lists.newArrayList();
        boolean hasHighPloidy = false;

        for(final SvVarData var : cluster.getSVs())
        {
            if(var.isReplicatedSv())
                continue;

            double ploidy = var.getSvData().ploidy();
            int i = 0;
            for(; i < ploidyList.size(); ++i)
            {
                Double otherPloidy = ploidyList.get(i);
                if(ploidy > otherPloidy)
                    break;
            }

            ploidyList.add(i, ploidy);
            indexSvList.add(i, var);

            if(ploidy >= DOUBLE_MINUTE_PLOIDY_THRESHOLD)
                hasHighPloidy = true;
        }

        if(!hasHighPloidy)
            return;

        boolean isPotentialDM = false;

        // keep track of any chromosome with high ploidy as an indication of false DMs and/or under-clustering
        List<String> highPloidyChromosomes = Lists.newArrayList();

        int svsAboveThreshold = 0;
        double minPloidyAboveThreshold = 0;
        for(int i = 0; i < ploidyList.size(); ++i)
        {
            Double ploidy = ploidyList.get(i);

            if(ploidy < DOUBLE_MINUTE_PLOIDY_THRESHOLD)
                break;

            ++svsAboveThreshold;
            minPloidyAboveThreshold = ploidy;

            final SvVarData var = indexSvList.get(i);

            for(int be = SVI_START; be <= SVI_END; ++be)
            {
                boolean useStart = isStart(be);

                if(var.isNullBreakend() && !useStart)
                    continue;

                if(!highPloidyChromosomes.contains(var.chromosome(useStart)))
                    highPloidyChromosomes.add(var.chromosome(useStart));
            }

            // check vs next
            if(i == ploidyList.size() - 1)
            {
                LOGGER.debug(String.format("cluster(%s count=%d) DM highPloidyCount(%d chr=%d) currentSV(%s) ploidy(%.2f) with no others",
                        cluster.id(), cluster.getUniqueSvCount(), svsAboveThreshold, highPloidyChromosomes.size(), var.posId(), ploidy));
                isPotentialDM = true;
                break;
            }
            else
            {
                double nextPloidy = ploidyList.get(i+1);

                if(nextPloidy * DOUBLE_MINUTE_PLOIDY_GAP_RATIO < ploidy)
                {
                    LOGGER.debug(String.format("cluster(%s count=%d) DM highPloidyCount(%d chr=%d) currentSV(%s) ploidy(%.2f) vs next(%.3f)",
                            cluster.id(), cluster.getUniqueSvCount(), svsAboveThreshold, highPloidyChromosomes.size(),
                            indexSvList.get(i).posId(), ploidy, nextPloidy));
                    isPotentialDM = true;
                    break;
                }
            }
        }

        if(isPotentialDM)
        {
            // check for high ploidy in other variants on the relevant chromosomes
            boolean otherClustersHaveHighPloidy = false;

            for(final String chromsome : highPloidyChromosomes)
            {
                final List<SvBreakend> breakendList = chrBreakendMap.get(chromsome);

                if(breakendList == null)
                    continue;

                for(final SvBreakend breakend : breakendList)
                {
                    if(breakend.getSV().getCluster() == cluster)
                        continue;

                    if(breakend.getSV().getSvData().ploidy() * DOUBLE_MINUTE_PLOIDY_GAP_RATIO >= minPloidyAboveThreshold)
                    {
                        otherClustersHaveHighPloidy = true;
                        break;
                    }
                }

                if(otherClustersHaveHighPloidy)
                    break;
            }

            final String dmAnnotation = otherClustersHaveHighPloidy ? CLUSTER_ANNONTATION_DM + "_Unclear" : CLUSTER_ANNONTATION_DM;
            cluster.addAnnotation(dmAnnotation);
        }

        /*
        if(cluster.isResolved() || !cluster.isFullyChained() || cluster.getChains().get(0).getLinkCount() <= 2)
            return;

        if(!cluster.getChains().get(0).couldFormLoop())
            return;

        final SvChain chain = cluster.getChains().get(0);

        // check for high copy number within this loop
        boolean hasHighCN = false;
        double cnTotal = 0;
        double maxCN = 0;
        double ploidyTotal = 0;

        for(final SvVarData var : chain.getSvList())
        {
            ploidyTotal += var.getSvData().ploidy();
            cnTotal += (var.copyNumber(true) + var.copyNumber(false)) * 0.5;
            maxCN = max(maxCN, var.copyNumber(true));
            maxCN = max(maxCN, var.copyNumber(false));

            if(var.copyNumber(true) >= DOUBLE_MINUTE_COPY_NUMBER_THRESHOLD || var.copyNumber(false) >= DOUBLE_MINUTE_COPY_NUMBER_THRESHOLD)
            {
                hasHighCN = true;
            }
        }

        if(hasHighCN)
        {
            LOGGER.debug(String.format("sample(%s) cluster(%d) chain(%d) links(%d) closed loop, copyNumber(avg=%.1f max=%.2f) ploidy(%.1f)",
                    sampleId, cluster.id(), chain.id(), chain.getLinkCount(),
                    cnTotal/chain.getSvList().size(), maxCN, ploidyTotal/chain.getSvList().size()));

            cluster.addAnnotation(CLUSTER_ANNONTATION_DM);
        }
        */
    }

    public static void findChainRepeatedSegments(final String sampleId, final SvCluster cluster, final SvChain chain)
    {
        if(!chain.hasReplicatedSVs())
            return;

        List<SvVarData> replicatedSVs = Lists.newArrayList();

        final List<SvVarData> svList = chain.getSvList();

        for (int i = 0; i < svList.size(); ++i)
        {
            final SvVarData var1 = svList.get(i);

            if(replicatedSVs.contains(var1))
                continue;

            for (int j = i + 1; j < svList.size(); ++j)
            {
                final SvVarData var2 = svList.get(j);

                if (!var1.equals(var2, true))
                    continue;

                replicatedSVs.add(var1);

                // look for repeated sections forwards or backwards from this point
                List<SvVarData> forwardRepeats = getRepeatedSvSequence(svList, i, j, true);

                boolean forwardSequence = false;

                if(!forwardRepeats.isEmpty())
                {
                    forwardSequence = true;
                    replicatedSVs.addAll(forwardRepeats);
                }
                else
                {
                    forwardSequence = false;
                    forwardRepeats = getRepeatedSvSequence(svList, i, j, false);
                }

                if(!forwardRepeats.isEmpty())
                {
                    replicatedSVs.addAll(forwardRepeats);

                    forwardRepeats.set(0, var1);

                    String svIds = var1.id();
                    for(int k = 1; k < forwardRepeats.size(); ++k)
                        svIds += ";" + forwardRepeats.get(k).id();

                    if(forwardRepeats.size() >= 4)
                    {
                        LOGGER.debug("sample({}) cluster({}) chain({}) {} sequence of {} SVs starting at index({}:{}) SV({})",
                                sampleId, cluster.id(), chain.id(), forwardSequence ? "forward" : "reverse",
                                forwardRepeats.size(), i, j, var1.id());

                        // ClusterId,ChainId,SequenceCount,VarIds,MatchDirection
                        LOGGER.debug("CF_REPEAT_SEQ: {},{},{},{},{},{}",
                                sampleId, cluster.id(), chain.id(), forwardRepeats.size(), svIds, forwardSequence);
                    }

                    break;
                }

                // no sequence found
            }
        }
    }

}
