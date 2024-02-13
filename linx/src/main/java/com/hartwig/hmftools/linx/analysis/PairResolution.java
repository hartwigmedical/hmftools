package com.hartwig.hmftools.linx.analysis;

import static java.lang.Math.abs;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.linx.LinxConfig.LNX_LOGGER;
import static com.hartwig.hmftools.linx.analysis.ClusterClassification.getSyntheticGapLength;
import static com.hartwig.hmftools.linx.analysis.ClusterClassification.getSyntheticLength;
import static com.hartwig.hmftools.linx.analysis.ClusterClassification.getSyntheticTiLength;
import static com.hartwig.hmftools.linx.analysis.SvUtilities.NO_LENGTH;
import static com.hartwig.hmftools.linx.analysis.SvUtilities.calcConsistency;
import static com.hartwig.hmftools.linx.chaining.ChainJcnLimits.jcnMatch;
import static com.hartwig.hmftools.linx.chaining.ChainUtils.reverseSectionOnBreakend;
import static com.hartwig.hmftools.linx.types.LinxConstants.MIN_SIMPLE_DUP_DEL_CUTOFF;
import static com.hartwig.hmftools.linx.types.LinxConstants.SHORT_TI_LENGTH;
import static com.hartwig.hmftools.linx.types.ResolvedType.DEL_TI;
import static com.hartwig.hmftools.linx.types.ResolvedType.DUP_TI;
import static com.hartwig.hmftools.linx.types.ResolvedType.FB_INV_PAIR;
import static com.hartwig.hmftools.linx.types.ResolvedType.NONE;
import static com.hartwig.hmftools.linx.types.ResolvedType.RECIP_INV;
import static com.hartwig.hmftools.linx.types.ResolvedType.RECIP_INV_DEL_DUP;
import static com.hartwig.hmftools.linx.types.ResolvedType.RECIP_INV_DUPS;
import static com.hartwig.hmftools.linx.types.ResolvedType.RECIP_TRANS;
import static com.hartwig.hmftools.linx.types.ResolvedType.RECIP_TRANS_DEL_DUP;
import static com.hartwig.hmftools.linx.types.ResolvedType.RECIP_TRANS_DUPS;
import static com.hartwig.hmftools.linx.types.ResolvedType.RESOLVED_FOLDBACK;
import static com.hartwig.hmftools.linx.types.ResolvedType.UNBAL_TRANS_TI;
import static com.hartwig.hmftools.linx.types.SvCluster.CLUSTER_ANNOT_DM;

import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.linx.chaining.ChainJcnLimits;
import com.hartwig.hmftools.linx.chaining.SvChain;
import com.hartwig.hmftools.linx.cn.LohEvent;
import com.hartwig.hmftools.linx.types.DbPair;
import com.hartwig.hmftools.linx.types.LinkedPair;
import com.hartwig.hmftools.linx.types.ResolvedType;
import com.hartwig.hmftools.linx.types.SvBreakend;
import com.hartwig.hmftools.linx.types.SvCluster;
import com.hartwig.hmftools.linx.types.SvVarData;

import org.jetbrains.annotations.NotNull;

public class PairResolution
{
    public static boolean isClusterPairType(final ResolvedType resolvedType)
    {
        return isClusterReciprocalType(resolvedType) || isClusterTemplatedInsertionType(resolvedType);
    }

    public static boolean isClusterReciprocalType(final ResolvedType resolvedType)
    {
        return (resolvedType == RECIP_TRANS || resolvedType == RECIP_TRANS_DEL_DUP || resolvedType == RECIP_TRANS_DUPS
                || resolvedType == RECIP_INV || resolvedType == RECIP_INV_DEL_DUP || resolvedType == RECIP_INV_DUPS
                || resolvedType == RESOLVED_FOLDBACK || resolvedType == FB_INV_PAIR);
    }

    public static boolean isClusterTemplatedInsertionType(final ResolvedType resolvedType)
    {
        return (resolvedType == DEL_TI || resolvedType == DUP_TI || resolvedType == UNBAL_TRANS_TI);
    }

    public static void classifyPairClusters(
            final SvCluster cluster, long longDelThreshold, long longDupThreshold, final Map<String,List<SvBreakend>> chrBreakendMap,
            boolean isGermline)
    {
        // classifies 2-event clusters based on the types of breakends and their orientations, and whether a long TI exists
        // treat existing chains as SVs - ie with 2 breakends from the open ends
        if(cluster.getChains().size() > 2 || cluster.getSglBreakendCount() > 0)
            return;

        if(cluster.hasAnnotation(CLUSTER_ANNOT_DM))
            return;

        // establish the nature of the breakends
        // first reduce SVs and/or chains to a set of breakends and record TI length(s)
        SvBreakend startBe1 = null;
        SvBreakend endBe1 = null;
        SvBreakend startBe2 = null;
        SvBreakend endBe2 = null;

        LinkedPair longTiLink = null;
        boolean uniformJcn = false;

        // establish the characteristics of this cluster:
        // - how many chains - 0, 1 or 2
        // - does it have a long TI
        // - what are the orientations and arms of its open breakends (from unlinked SVs or chain ends)

        List<SvChain> clusterChains = Lists.newArrayList(cluster.getChains());
        SvVarData unchainedSv = !cluster.getUnlinkedSVs().isEmpty() ? cluster.getUnlinkedSVs().get(0) : null;
        int unchainedSvCount = cluster.getUnlinkedSVs().size();
        List<LinkedPair> longTiLinks = Lists.newArrayList();
        boolean existingChainModified = false;

        if(!clusterChains.isEmpty())
        {
            longTiLinks.addAll(cluster.getChains().get(0).getLinkedPairs().stream()
                    .filter(x -> x.baseLength() > SHORT_TI_LENGTH).collect(Collectors.toList()));

            if(cluster.getChains().size() == 2)
            {
                longTiLinks.addAll(cluster.getChains().get(1).getLinkedPairs().stream()
                        .filter(x -> x.baseLength() > SHORT_TI_LENGTH).collect(Collectors.toList()));
            }

            if(longTiLinks.size() > 1)
                return;

            if(longTiLinks.size() == 1)
                longTiLink = longTiLinks.get(0);
        }

        boolean isSingleChain = cluster.isFullyChained(false) && clusterChains.size() == 1;

        // first look for a single chain with only short TIs - this can be classified as a synthetic DEL or DUP
        if(isSingleChain && longTiLink == null)
        {
            //  a single consistent chain with only short TIs may be resolved as a synthetic DEL or DUP or something else,
            classifySyntheticDelDups(cluster, longDelThreshold, longDupThreshold);
            return;
        }

        if(!isGermline && longTiLinks.stream().anyMatch(x -> isInterruptedTI(x, chrBreakendMap)))
            return;

        // otherwise test whether a single chain be split at the long TI to make 2 chains and a balance translocation
        if(isSingleChain && longTiLink != null)
        {
            // go into regular 2 break logic after first breaking the chain at this long or internal TI
            SvChain chain = clusterChains.get(0);
            clusterChains.clear();
            existingChainModified = true;

            final List<LinkedPair> pairs = chain.getLinkedPairs();

            if(pairs.size() > 1)
            {
                SvChain newChain = new SvChain(chain.id());

                if(pairs.get(0) == longTiLink || pairs.get(pairs.size() - 1) == longTiLink)
                {
                    unchainedSvCount = 1;
                    unchainedSv = pairs.get(0) == longTiLink ? pairs.get(0).first() : pairs.get(pairs.size() - 1).second();

                    for(LinkedPair pair : pairs)
                    {
                        if(pair != longTiLink)
                        {
                            newChain.addLink(pair, false);
                        }
                    }

                    clusterChains.add(newChain);
                }
                else
                {
                    for(LinkedPair pair : pairs)
                    {
                        if(pair == longTiLink)
                        {
                            clusterChains.add(newChain);

                            newChain = new SvChain(chain.id() + 1);
                            continue;
                        }

                        newChain.addLink(pair, false);
                    }

                    clusterChains.add(newChain);
                }
            }
        }

        if(clusterChains.isEmpty())
        {
            if(cluster.getSvCount() != 2)
                return;

            final SvVarData var1 = cluster.getSV(0);
            final SvVarData var2 = cluster.getSV(1);

            startBe1 = var1.getBreakend(true);
            endBe1 = var1.getBreakend(false);
            startBe2 = var2.getBreakend(true);
            endBe2 = var2.getBreakend(false);
            uniformJcn = ChainJcnLimits.jcnMatch(var1, var2);
        }
        else if(cluster.isFullyChained(false) && clusterChains.size() == 2)
        {
            SvChain chain1 = clusterChains.get(0);
            SvChain chain2 = clusterChains.get(1);

            startBe1 = chain1.getOpenBreakend(true);
            endBe1 = chain1.getOpenBreakend(false);
            startBe2 = chain2.getOpenBreakend(true);
            endBe2 = chain2.getOpenBreakend(false);
            uniformJcn = jcnMatch(chain1.jcn(), chain1.jcnUncertainty(), chain2.jcn(), chain2.jcnUncertainty());
        }
        else if(clusterChains.size() == 1 && unchainedSvCount == 1)
        {
            // check for a single SV and a chain
            SvVarData var = unchainedSv;
            SvChain chain = clusterChains.get(0);

            startBe1 = chain.getOpenBreakend(true);
            endBe1 = chain.getOpenBreakend(false);

            startBe2 = var.getBreakend(true);
            endBe2 = var.getBreakend(false);

            uniformJcn = jcnMatch(var.jcn(), var.jcnUncertainty(), chain.jcn(), chain.jcnUncertainty());
        }
        else
        {
            return;
        }

        if(!isGermline && areInterruptedBreakends(startBe1, endBe1, startBe2, endBe2, chrBreakendMap))
            return;

        if(startBe1.chromosome().equals(endBe1.chromosome()) && startBe2.chromosome().equals(endBe2.chromosome())
        && startBe1.chromosome().equals(startBe2.chromosome()))
        {
            if(!isSingleChain)
                return;

            // same arm events - first 2 INVs
            if(startBe1.orientation() == endBe1.orientation() && startBe2.orientation() == endBe2.orientation()
            && startBe1.orientation() != startBe2.orientation())
            {
                classifyInversionPairClusters(
                        cluster, longDelThreshold, longDupThreshold,
                        startBe1, endBe1, startBe2, endBe2, uniformJcn, longTiLink);
            }
        }
        else
        {
            // check for translocation type events
            classifyTranslocationPairClusters(
                    cluster, longDelThreshold, longDupThreshold, startBe1, endBe1, startBe2, endBe2,
                    uniformJcn, longTiLink, existingChainModified ? clusterChains: null);
        }
    }

    private static void classifySyntheticDelDups(final SvCluster cluster, long longDelThreshold, long longDupThreshold)
    {
        if(!cluster.isFullyChained(true) || cluster.getChains().size() != 1 || cluster.getSglBreakendCount() > 0)
            return;

        SvChain chain = cluster.getChains().get(0);

        final SvBreakend chainStart = chain.getOpenBreakend(true);
        final SvBreakend chainEnd = chain.getOpenBreakend(false);

        if(!chainStart.chromosome().equals(chainEnd.chromosome()) || chainStart.arm() != chainEnd.arm())
            return;

        if(chainStart.orientation() == chainEnd.orientation())
            return;

        boolean faceAway = (chainStart.position() < chainEnd.position()) == (chainStart.orientation() == 1);

        long totalChainLength = chain.getLength(false);
        long syntheticLength = abs(chainStart.position() - chainEnd.position());

        ResolvedType resolvedType = faceAway ? ResolvedType.DEL : ResolvedType.DUP;

        LNX_LOGGER.debug("cluster({}) chain(links=({} len={} synLen({}) marked as {}",
                cluster.id(), chain.getLinkCount(), totalChainLength, syntheticLength, resolvedType);

        boolean withinLongThreshold = resolvedType == ResolvedType.DEL ?
                syntheticLength < longDelThreshold : syntheticLength < longDupThreshold;

        boolean resolved = withinLongThreshold;
        cluster.setResolved(resolved, resolvedType);
        setPairLengthData(cluster);
    }

    private static void setPairLengthData(final SvCluster cluster)
    {
        // format is: SyntheticLength (DEL or DUP), TI or DB or OverlapLength, SV Gap Length (if applicable)
        if(cluster.getChains().size() != 1)
            return;

        long syntheticLength = getSyntheticLength(cluster);
        long otherLength = getSyntheticTiLength(cluster);
        long gapLength = getSyntheticGapLength(cluster);

        cluster.addAnnotation(String.format("PairLen=%d;%d;%d", syntheticLength, otherLength, gapLength));
    }

    private static boolean isLohBoundedTi(final LinkedPair pair)
    {
        for(final LohEvent lohEvent : pair.first().getCluster().getLohEvents())
        {
            if(lohEvent.getBreakend(true) == pair.firstBreakend() || lohEvent.getBreakend(true) == pair.secondBreakend())
                return true;
            else if(lohEvent.getBreakend(false) == pair.firstBreakend() || lohEvent.getBreakend(false) == pair.secondBreakend())
                return true;
        }

        return false;
    }

    private static void classifyTranslocationPairClusters(
            SvCluster cluster, long longDelThreshold, long longDupThreshold,
            SvBreakend startBe1, SvBreakend endBe1, SvBreakend startBe2, SvBreakend endBe2,
            boolean uniformJcn, final LinkedPair longestTiPair,
            List<SvChain> rearrangedChains)
    {
        // first check consistency
        int consistency = calcConsistency(startBe1) + calcConsistency(endBe1) + calcConsistency(startBe2) + calcConsistency(endBe2);

        if(consistency != 0)
            return;

        // reciprocal translocations require each of the breakend pairs to start and end on the same pair of arms
        SvBreakend arm1Be1 = null;
        SvBreakend arm2Be1 = null;
        SvBreakend arm1Be2 = null;
        SvBreakend arm2Be2 = null;

        if(!startBe1.getChrArm().equals(endBe1.getChrArm()) && !startBe2.getChrArm().equals(endBe2.getChrArm()))
        {
            arm1Be1 = startBe1;
            arm2Be1 = endBe1;

            if(startBe2.getChrArm().equals(arm1Be1.getChrArm()) && endBe2.getChrArm().equals(arm2Be1.getChrArm()))
            {
                arm1Be2 = startBe2;
                arm2Be2 = endBe2;
            }
            else if(endBe2.getChrArm().equals(arm1Be1.getChrArm()) && startBe2.getChrArm().equals(arm2Be1.getChrArm()))
            {
                arm1Be2 = endBe2;
                arm2Be2 = startBe2;
            }
        }

        if(arm1Be1 == null || arm1Be2 == null || arm2Be1 == null || arm2Be2 == null)
        {
            // check for a synthetic BND with a long TI
            if(longestTiPair != null && cluster.getChains().size() == 1 && cluster.getUnlinkedSVs().isEmpty())
            {
                final SvChain chain = cluster.getChains().get(0);

                if(!chain.getOpenBreakend(SE_START).getChrArm().equals(chain.getOpenBreakend(SE_END).getChrArm()))
                {
                    cluster.setResolved(false, UNBAL_TRANS_TI);
                }
            }

            return;
        }

        /* establish configuration:
            - Deletion bridge on both arms
            - Facing on one arm.  Deletion Bridge on other arm
            - Facing on 2 arms
            - 1 arm shared only. Facing or deletion bridge
         */

        // check for linked, short DBs at both ends
        DbPair arm1Db1 = arm1Be1.getDBLink();
        DbPair arm2Db1 = arm2Be1.getDBLink();

        boolean arm1MatchingDB = arm1Db1 != null && arm1Db1 == arm1Be2.getDBLink();
        boolean arm2MatchingDB = arm2Db1 != null && arm2Db1 == arm2Be2.getDBLink();

        // basic reciprocal translocation
        ResolvedType resolvedType;
        boolean isResolved = true;
        if(arm1MatchingDB && arm2MatchingDB)
        {
            resolvedType = RECIP_TRANS;
            isResolved = (arm1Db1.length() <= longDelThreshold && arm2Db1.length() <= longDelThreshold);

            cluster.addAnnotation(String.format("PairLen=%d;%d;%d", arm1Db1.length(), arm2Db1.length(), NO_LENGTH));
        }
        else
        {
            if(longestTiPair == null)
                return;

            // set synthetic lengths prior to any chain reconfiguration
            setPairLengthData(cluster);

            boolean lohBoundedTi = isLohBoundedTi(longestTiPair); // length no longer a classifying factor

            if(!arm1MatchingDB && !arm2MatchingDB)
            {
                resolvedType = lohBoundedTi && uniformJcn ? DUP_TI : RECIP_TRANS_DUPS;
            }
            else
            {
                resolvedType = lohBoundedTi && uniformJcn ? DEL_TI : RECIP_TRANS_DEL_DUP;
            }

            // test the overlap and DB lengths to determine whether this cluster is resolved (ie protected from clustering)
            if(longestTiPair.baseLength() > longDupThreshold)
            {
                isResolved = false;
            }
            else
            {
                final SvBreakend nonLongTiBe1 = longestTiPair.hasBreakend(startBe1) ? endBe1 : startBe1;
                final SvBreakend nonLongTiBe2 = longestTiPair.hasBreakend(startBe2) ? endBe2 : startBe2;
                long breakendDistance = abs(nonLongTiBe1.position() - nonLongTiBe2.position());

                if(resolvedType == DEL_TI || resolvedType == RECIP_TRANS_DEL_DUP)
                {
                    if(breakendDistance > longDelThreshold)
                        isResolved = false;
                }
                else
                {
                    // other breakends overlap in another DUP
                    if(breakendDistance > longDupThreshold)
                        isResolved = false;
                }
            }
        }

        LNX_LOGGER.debug("cluster({}) longestTI({}) resolvedType({}) isResolved({})",
                cluster.id(), longestTiPair != null ? longestTiPair.baseLength() : "none", resolvedType, isResolved);

        if(rearrangedChains != null)
        {
            LNX_LOGGER.debug("cluster({}) splitting existing chain into {} for resolvedType({})",
                    cluster.id(), rearrangedChains.size(), cluster.getResolvedType());

            final SvChain existingChain = cluster.getChains().get(0);
            cluster.dissolveLinksAndChains();

            for(SvChain newChain : rearrangedChains)
            {
                newChain.setJcnData(existingChain.jcn(), existingChain.jcnUncertainty());
                cluster.addChain(newChain, true);
            }
        }

        cluster.setResolved(isResolved, resolvedType);
    }

    private static void classifyInversionPairClusters(
            SvCluster cluster, long longDelThreshold, long longDupThreshold,
            SvBreakend startBe1, SvBreakend endBe1, SvBreakend startBe2, SvBreakend endBe2,
            boolean uniformJcn, @NotNull final LinkedPair longestTiPair)
    {
        /* establish configuration:
        1. FB_INV_PAIR - facing foldbacks - +ve INV positions 3,4, -ve INV positions 1,2

        2. Two DELS - 2x DBs, +ve INV positions 1,3, -ve INV positions 2,4
        - RECIP_INV if both DBs are present otherwise leave as-is

        3. One INV encloses the other - has 1x DB, 1x TI which must be long (>1K) since otherwise would be shard and just a synthetic DEL,
         +ve INV positions 1,3, -ve INV positions 2,4
        - if the TI is LOH bounded and the INVs uniform ploidy (now JCN), this is a DEL_TI
        - if the inner INV length < 100K, this is a RESOLVED_FOLDBACK
        - otherwise is a RECIP_INV_DEL_DUP, and the chain is rearranged so the long TI is formed

        4. Overlapping INVs with 2 facing breakends - 2x TIs, +ve INV positions 2,4, -ve INV positions 1,3
        - DOUBLE_MINUTE if satisfies DM rules
        - if the TI is LOH bounded and the INVs uniform ploidy), this is a DUP_TI
        - otherwise is a RECIP_INV_DUP2, and the chain is rearranged so the long TI is formed
        */

        // assign relative positions for easier comparison
        SvBreakend lowerBe1 = startBe1.position() < endBe1.position() ? startBe1 : endBe1;
        SvBreakend upperBe1 = lowerBe1 == startBe1 ? endBe1 : startBe1;
        SvBreakend lowerBe2 = startBe2.position() < endBe2.position() ? startBe2 : endBe2;
        SvBreakend upperBe2 = lowerBe2 == startBe2 ? endBe2 : startBe2;

        // 1. Facing Inversions
        if(upperBe1.position() < lowerBe2.position() || upperBe2.position() < lowerBe1.position())
        {
            // no overlap and facing away is unclassified
            boolean inversionsFacing = (upperBe1.position() < lowerBe2.position()) == (upperBe1.orientation() == -1);

            if(!inversionsFacing)
                return;

            // inversions are facing but not overlapping
            cluster.setResolved(false, FB_INV_PAIR);
            return;
        }

        // check for linked, short DBs at both ends
        DbPair lowerDb1 = lowerBe1.getDBLink();
        DbPair upperDb1 = upperBe1.getDBLink();

        boolean lowerMatchingDB = lowerDb1 != null && lowerDb1 == lowerBe2.getDBLink();
        boolean upperMatchingDB = upperDb1 != null && upperDb1 == upperBe2.getDBLink();

        // 2. Reciprocal INV case
        if(lowerMatchingDB && upperMatchingDB)
        {
            // basic reciprocal inversion
            boolean isResolved = (lowerDb1.length() <= longDelThreshold && upperDb1.length() <= longDelThreshold);
            cluster.setResolved(isResolved, RECIP_INV);
            cluster.addAnnotation(String.format("PairLen=%d;%d;%d", lowerDb1.length(), upperDb1.length(), NO_LENGTH));
            return;
        }
        else if((lowerBe1.orientation() == 1 && lowerBe1.position() < lowerBe2.position() && upperBe1.position() < upperBe2.position())
        || (lowerBe2.orientation() == 1 && lowerBe2.position() < lowerBe1.position() && upperBe2.position() < upperBe1.position()))
        {
            // still possibly a reciprocal INV but interrupted by other breakends
            return;
        }

        // set prior to any chain reconfiguration
        setPairLengthData(cluster);

        boolean lohBoundedTi = isLohBoundedTi(longestTiPair);

        ResolvedType resolvedType = NONE;

        // 3. One INV encloses the other - the DEL scenario
        if((lowerBe1.position() < lowerBe2.position() && upperBe1.position() > upperBe2.position())
        || (lowerBe2.position() < lowerBe1.position() && upperBe2.position() > upperBe1.position()))
        {
            if(lohBoundedTi && uniformJcn)
            {
                resolvedType = DEL_TI;
            }
            else
            {
                long innerInversionLength = min(upperBe1.position() - lowerBe1.position(), upperBe2.position() - lowerBe2.position());

                if(innerInversionLength < MIN_SIMPLE_DUP_DEL_CUTOFF)
                {
                    resolvedType = RESOLVED_FOLDBACK;
                }
                else
                {
                    // re-arrange the links so that the inner inversion links using its outer breakend to make a longer link
                    SvChain chain = cluster.getChains().get(0);

                    // work out which is the inner INV then switch around the end of it which is part of the long TI
                    SvBreakend breakendToSwitch = null;

                    if(lowerBe1.position() < lowerBe2.position())
                    {
                        // pair 2 is the inner pair
                        breakendToSwitch = longestTiPair.hasBreakend(lowerBe2) ? lowerBe2 : upperBe2;
                    }
                    else
                    {
                        breakendToSwitch = longestTiPair.hasBreakend(lowerBe1) ? lowerBe1 : upperBe1;
                    }

                    if(reverseSectionOnBreakend(chain, breakendToSwitch))
                    {
                        LNX_LOGGER.debug("cluster({}) reconfigured chain:", cluster.id());
                        chain.logLinks();
                    }

                    resolvedType = RECIP_INV_DEL_DUP;
                }
            }
        }
        else
        {
            // otherwise the inner breakends overlap
            if(lohBoundedTi && uniformJcn)
            {
                resolvedType = DUP_TI;
            }
            else
            {
                // flip the chain on the inner breakend which is in the long TI
                SvChain chain = cluster.getChains().get(0);
                SvBreakend breakendToSwitch = null;

                if(lowerBe1.position() > lowerBe2.position())
                {
                    // lower breakend 1 is one of the inner breakends
                    breakendToSwitch = longestTiPair.hasBreakend(lowerBe1) ? lowerBe1 : longestTiPair.getOtherBreakend(lowerBe1);
                }
                else
                {
                    breakendToSwitch = longestTiPair.hasBreakend(lowerBe2) ? lowerBe2 : longestTiPair.getOtherBreakend(lowerBe1);
                }

                if(reverseSectionOnBreakend(chain, breakendToSwitch))
                {
                    LNX_LOGGER.debug("cluster({}) reconfigured chain:", cluster.id());
                    chain.logLinks();
                }

                resolvedType = RECIP_INV_DUPS;
            }
        }

        boolean isResolved = false;

        // test DEL and DUP lengths vs thresholds to determine whether the cluster is protected
        long longestTiLength = cluster.getChains().get(0).getLinkedPairs().stream().mapToLong(LinkedPair::positionDistance).max().getAsLong();
        final SvBreakend chainStart = cluster.getChains().get(0).getOpenBreakend(true);
        final SvBreakend chainEnd = cluster.getChains().get(0).getOpenBreakend(false);
        long nonTiDistance = abs(chainStart.position() - chainEnd.position());

        if(resolvedType == RECIP_INV_DUPS || resolvedType == DUP_TI)
        {
            isResolved = longestTiLength <= longDupThreshold && nonTiDistance <= longDupThreshold;
        }
        else if(resolvedType == RECIP_INV_DEL_DUP || resolvedType == DEL_TI)
        {
            isResolved = longestTiLength <= longDupThreshold && nonTiDistance <= longDelThreshold;
        }
        else
        {
            isResolved = false;
        }

        LNX_LOGGER.debug("cluster({}) longestTI({}) resolvedType({}) isResolved({})",
                cluster.id(), longestTiPair != null ? longestTiPair.baseLength() : "none", resolvedType, isResolved);

        cluster.setResolved(isResolved, resolvedType);
    }

    private static boolean isInterruptedTI(final LinkedPair pair, final Map<String,List<SvBreakend>> chrBreakendMap)
    {
        return isInterruptedBreakendPair(pair.getBreakend(true), pair.getBreakend(false), chrBreakendMap);
    }

    private static boolean areInterruptedBreakends(
            final SvBreakend startBe1, final SvBreakend endBe1, final SvBreakend startBe2, final SvBreakend endBe2,
            final Map<String,List<SvBreakend>> chrBreakendMap)
    {
        if(startBe1.getChrArm().equals(startBe2.getChrArm()) && isInterruptedBreakendPair(startBe1, startBe2, chrBreakendMap))
            return true;
        else if(startBe1.getChrArm().equals(endBe2.getChrArm()) && isInterruptedBreakendPair(startBe1, endBe2, chrBreakendMap))
            return true;
        else if(endBe1.getChrArm().equals(startBe2.getChrArm()) && isInterruptedBreakendPair(endBe1, startBe2, chrBreakendMap))
            return true;
        else if(endBe1.getChrArm().equals(endBe2.getChrArm()) && isInterruptedBreakendPair(endBe1, endBe2, chrBreakendMap))
            return true;

        return false;
    }

    private static boolean isInterruptedBreakendPair(
            final SvBreakend be1, final SvBreakend be2, final Map<String,List<SvBreakend>> chrBreakendMap)
    {
        // look between these breakends for any non-simple, non-contained SV
        final List<SvBreakend> breakendList = chrBreakendMap.get(be1.chromosome());

        final SvBreakend lowerBe = be1.position() < be2.position() ? be1 : be2;
        final SvBreakend upperBe = lowerBe == be1 ? be2 : be1;

        int startIndex = lowerBe.getChrPosIndex();
        int endIndex = upperBe.getChrPosIndex();

        for(int i = startIndex + 1; i < endIndex; ++i)
        {
            final SvVarData var = breakendList.get(i).getSV();

            if(var == be1.getSV() || var == be2.getSV() || var.getCluster() == be1.getSV().getCluster())
                continue;

            if(var.getCluster().getSvCount() > 1)
                return true;

            if(!var.isSimpleType())
                return true;

            // lastly confirm the the simple SV is confined to this segment
            if(var.getBreakend(true).getChrPosIndex() < startIndex || var.getBreakend(false).getChrPosIndex() > endIndex)
                return true;
        }

        return false;
    }

}
