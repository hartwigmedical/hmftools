package com.hartwig.hmftools.linx.analysis;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.round;

import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.BND;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.DEL;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.DUP;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.INV;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.SGL;
import static com.hartwig.hmftools.linx.analysis.SvUtilities.calcConsistency;
import static com.hartwig.hmftools.linx.chaining.ChainPloidyLimits.calcPloidyUncertainty;
import static com.hartwig.hmftools.linx.chaining.ChainPloidyLimits.ploidyMatch;
import static com.hartwig.hmftools.linx.types.ResolvedType.COMPLEX;
import static com.hartwig.hmftools.linx.types.ResolvedType.DEL_TI;
import static com.hartwig.hmftools.linx.types.ResolvedType.DUP_TI;
import static com.hartwig.hmftools.linx.types.ResolvedType.FB_INV_PAIR;
import static com.hartwig.hmftools.linx.types.ResolvedType.NONE;
import static com.hartwig.hmftools.linx.types.ResolvedType.PAIR_OTHER;
import static com.hartwig.hmftools.linx.types.ResolvedType.RECIP_INV;
import static com.hartwig.hmftools.linx.types.ResolvedType.RECIP_INV_DEL_DUP;
import static com.hartwig.hmftools.linx.types.ResolvedType.RECIP_INV_DUPS;
import static com.hartwig.hmftools.linx.types.ResolvedType.RECIP_TRANS;
import static com.hartwig.hmftools.linx.types.ResolvedType.RECIP_TRANS_DEL_DUP;
import static com.hartwig.hmftools.linx.types.ResolvedType.RECIP_TRANS_DUPS;
import static com.hartwig.hmftools.linx.types.ResolvedType.RESOLVED_FOLDBACK;
import static com.hartwig.hmftools.linx.types.SvaConstants.SHORT_TI_LENGTH;

import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.linx.chaining.SvChain;
import com.hartwig.hmftools.linx.cn.LohEvent;
import com.hartwig.hmftools.linx.types.ResolvedType;
import com.hartwig.hmftools.linx.types.SvBreakend;
import com.hartwig.hmftools.linx.types.SvCluster;
import com.hartwig.hmftools.linx.types.SvLinkedPair;
import com.hartwig.hmftools.linx.types.SvVarData;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class PairResolution
{
    private static final Logger LOGGER = LogManager.getLogger(PairResolution.class);

    public static boolean isClusterPairType(final ResolvedType resolvedType)
    {
        return (resolvedType == RECIP_TRANS || resolvedType == RECIP_TRANS_DEL_DUP || resolvedType == RECIP_TRANS_DUPS
                || resolvedType == RECIP_INV || resolvedType == RECIP_INV_DEL_DUP || resolvedType == RECIP_INV_DUPS
                || resolvedType == DUP_TI || resolvedType == DEL_TI
                || resolvedType == FB_INV_PAIR || resolvedType == RESOLVED_FOLDBACK || resolvedType == FB_INV_PAIR);
    }

    public static boolean isLohBoundedTi(final SvLinkedPair pair)
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

    public static void classifyPairClusters(
            SvCluster cluster, boolean isFinal, int proximityThreshold, long longDelThreshold, long longDupThreshold)
    {
        // classifies 2-event clusters based on the types of breakends and their orientations
        // treat existing chains as SVs - ie with 2 breakends from the open ends
        if(cluster.getChains().size() > 2 || cluster.getSglBreakendCount() > 0)
            return;

        // establish the nature of the breakends
        // first reduce SVs and/or chains to a set of breakends and record TI length(s)
        SvBreakend startBe1 = null;
        SvBreakend endBe1 = null;
        SvBreakend startBe2 = null;
        SvBreakend endBe2 = null;

        SvLinkedPair longTiLink = null;
        boolean uniformPloidy = false;

        // 2x translocations can exist in the form of 2 BNDs,

        List<SvChain> clusterChains = Lists.newArrayList(cluster.getChains());
        SvVarData unchainedSv = !cluster.getUnlinkedSVs().isEmpty() ? cluster.getUnlinkedSVs().get(0) : null;
        int unchainedSvCount = cluster.getUnlinkedSVs().size();
        boolean existingChainModified = false;

        if(!clusterChains.isEmpty())
        {
            List<SvLinkedPair> longTiLinks = cluster.getChains().get(0).getLinkedPairs().stream()
                    .filter(x -> x.length() > SHORT_TI_LENGTH).collect(Collectors.toList());

            if(cluster.getChains().size() == 2)
            {
                longTiLinks.addAll(cluster.getChains().get(1).getLinkedPairs().stream()
                        .filter(x -> x.length() > SHORT_TI_LENGTH).collect(Collectors.toList()));
            }

            if(longTiLinks.size() > 1)
                return;

            if(longTiLinks.size() == 1)
                longTiLink = longTiLinks.get(0);
        }

        boolean isSingleChain = cluster.isFullyChained(false) && clusterChains.size() == 1;

        // first handle the single chain case either handling it as a synthetic with only short TIs, or by splitting it
        if(isSingleChain && longTiLink == null)
        {
            //  a single chain with only short TIs may be resovled as a synthetic DEL or DUP or something else,
            // but not a pair-type cluster
            classifySyntheticDelDups(cluster, longDelThreshold, longDupThreshold);
            return;
        }

            // first handle the single chain case either handling it as a synthetic with only short TIs, or by splitting it
        if(isSingleChain)
        {
            // go into regular 2 break logic after first breaking the chain at this long TI
            SvChain chain = clusterChains.get(0);
            clusterChains.clear();
            existingChainModified = true;

            final List<SvLinkedPair> pairs = chain.getLinkedPairs();

            if(pairs.size() > 1)
            {
                SvChain newChain = new SvChain(chain.id());

                if (pairs.get(0) == longTiLink || pairs.get(pairs.size() - 1) == longTiLink)
                {
                    unchainedSvCount = 1;
                    unchainedSv = pairs.get(0) == longTiLink ? pairs.get(0).first() : pairs.get(pairs.size() - 1).second();

                    for (SvLinkedPair pair : pairs)
                    {
                        if (pair != longTiLink)
                        {
                            newChain.addLink(pair, false);
                        }
                    }

                    clusterChains.add(newChain);
                }
                else
                {
                    for (SvLinkedPair pair : pairs)
                    {
                        if (pair == longTiLink)
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
            uniformPloidy = ploidyMatch(var1.ploidy(), var1.ploidyUncertainty(), var2.ploidy(), var2.ploidyUncertainty());
        }
        else if(cluster.isFullyChained(false) && clusterChains.size() == 2)
        {
            SvChain chain1 = clusterChains.get(0);
            SvChain chain2 = clusterChains.get(1);

            startBe1 = chain1.getOpenBreakend(true);
            endBe1 = chain1.getOpenBreakend(false);
            startBe2 = chain2.getOpenBreakend(true);
            endBe2 = chain2.getOpenBreakend(false);
            uniformPloidy = ploidyMatch(chain1.ploidy(), chain1.ploidyUncertainty(), chain2.ploidy(), chain2.ploidyUncertainty());
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

            uniformPloidy = ploidyMatch(var.ploidy(), var.ploidyUncertainty(), chain.ploidy(), chain.ploidyUncertainty());
        }
        else
        {
            return;
        }

        if(startBe1.getChrArm().equals(endBe1.getChrArm()) && startBe2.getChrArm().equals(endBe2.getChrArm())
        && startBe1.getChrArm().equals(startBe2.getChrArm()))
        {
            if(!isSingleChain)
                return;

            // same arm events - first 2 INVs
            if(startBe1.orientation() == endBe1.orientation() && startBe2.orientation() == endBe2.orientation()
            && startBe1.orientation() != startBe2.orientation())
            {
                classifyInversionPairClusters(
                        cluster, isFinal, proximityThreshold, longDupThreshold,
                        startBe1, endBe1, startBe2, endBe2, uniformPloidy, longTiLink);
            }
            else if(startBe1.orientation() != endBe1.orientation() && startBe2.orientation() != endBe2.orientation())
            {
                // next a DUP and DEL
                boolean firstIsDel = (startBe1.orientation() == 1) == (startBe1.position() < endBe1.position());
                boolean secondIsDel = (startBe2.orientation() == 1) == (startBe2.position() < endBe2.position());

                if(secondIsDel != firstIsDel)
                {
                    classifyDelDupPairClusters(cluster, startBe1, endBe1, startBe2, endBe2, uniformPloidy);
                }
            }
        }
        else
        {
            // check for translocation type events
            classifyTranslocationPairClusters(
                    cluster, isFinal, proximityThreshold, longDupThreshold, startBe1, endBe1, startBe2, endBe2,
                    uniformPloidy, longTiLink, existingChainModified ? clusterChains: null);
        }
    }

    private static void classifySyntheticDelDups(SvCluster cluster, long longDelThreshold, long longDupThreshold)
    {
        if(!cluster.isFullyChained(true) || cluster.getChains().size() != 1 || cluster.getSglBreakendCount() > 0)
            return;

        SvChain chain = cluster.getChains().get(0);

        final SvBreakend startBreakend = chain.getOpenBreakend(true);
        final SvBreakend endBreakend = chain.getOpenBreakend(false);

        if(!startBreakend.chromosome().equals(endBreakend.chromosome()) || startBreakend.arm() != endBreakend.arm())
            return;

        if(startBreakend.orientation() == endBreakend.orientation())
            return;

        boolean faceAway = (startBreakend.position() < endBreakend.position()) == (startBreakend.orientation() == 1);

        int totalChainLength = chain.getLength(false);
        long syntheticLength = abs(startBreakend.position() - endBreakend.position());

        ResolvedType resolvedType = faceAway ? ResolvedType.DEL : ResolvedType.DUP;

        LOGGER.debug("cluster({}) chain(links=({} len={} synLen({}) marked as {}",
                cluster.id(), chain.getLinkCount(), totalChainLength, syntheticLength, resolvedType);

        boolean withinLongThreshold = false;

        if(resolvedType == ResolvedType.DEL)
            withinLongThreshold = syntheticLength < longDelThreshold;
        else if(resolvedType == ResolvedType.DUP)
            withinLongThreshold = syntheticLength < longDupThreshold;

        boolean resolved = withinLongThreshold;
        cluster.setResolved(resolved, resolvedType);
    }

    private static void classifyTranslocationPairClusters(
            SvCluster cluster, boolean isFinal, int proximityThreshold, long longDupThreshold,
            SvBreakend startBe1, SvBreakend endBe1, SvBreakend startBe2, SvBreakend endBe2,
            boolean uniformPloidy, final SvLinkedPair longestTiPair,
            List<SvChain> rearrangedChains)
    {
        // first check consistency
        int consistency = calcConsistency(startBe1) + calcConsistency(endBe1) + calcConsistency(startBe2) + calcConsistency(endBe2);

        if(consistency != 0)
            return;

        // reciprocal translocations require each of the breakend pairs to start and end on the same pair of arms

        if(startBe1.getChrArm().equals(endBe1.getChrArm()) || startBe2.getChrArm().equals(endBe2.getChrArm()))
        {
            // one chain start and ends on the same arm, the other doesn't
            cluster.setResolved(false, PAIR_OTHER);
            return;
        }

        // assign to arms for easier comparison
        SvBreakend arm1Be1 = startBe1;
        SvBreakend arm2Be1 = endBe1;
        SvBreakend arm1Be2 = null;
        SvBreakend arm2Be2 = null;

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
        else
        {
            cluster.setResolved(false, PAIR_OTHER);
            return;
        }

        /* establish configuration:
            - Deletion bridge on both arms
            - Facing on one arm.  Deletion Bridge on other arm
            - Facing on 2 arms
            - 1 arm shared only. Facing or deletion bridge
         */

        // check for linked, short DBs at both ends
        SvLinkedPair arm1Db1 = arm1Be1.getDBLink();
        SvLinkedPair arm2Db1 = arm2Be1.getDBLink();

        boolean arm1MatchingDB = arm1Db1 != null && arm1Db1 == arm1Be2.getDBLink();
        boolean arm2MatchingDB = arm2Db1 != null && arm2Db1 == arm2Be2.getDBLink();

        // basic reciprocal translocation
        ResolvedType resolvedType;
        boolean isResolved = false;
        if(arm1MatchingDB && arm2MatchingDB)
        {
            if(!isFinal && (arm1Db1.length() > proximityThreshold || arm2Db1.length() > proximityThreshold))
                return;

            resolvedType = RECIP_TRANS;
            isResolved = true;
        }
        else
        {
            if (longestTiPair == null)
                return;

            boolean longOrLohBoundedTi = isLohBoundedTi(longestTiPair) || longestTiPair.length() > longDupThreshold;

            if (!arm1MatchingDB && !arm2MatchingDB)
            {
                resolvedType = longOrLohBoundedTi && uniformPloidy ? DUP_TI : RECIP_TRANS_DUPS;
            }
            else
            {
                resolvedType = longOrLohBoundedTi && uniformPloidy ? DEL_TI : RECIP_TRANS_DEL_DUP;
            }

            isResolved = (resolvedType == RECIP_TRANS_DUPS || resolvedType == RECIP_TRANS_DEL_DUP);
        }

        LOGGER.debug("cluster({}) longestTI({}) marked as {}",
                cluster.id(), longestTiPair != null ? longestTiPair.length() : "none", resolvedType);

        if(rearrangedChains != null)
        {
            LOGGER.debug("cluster({}) splitting existing chain into {} for resolvedType({})",
                    cluster.id(), rearrangedChains.size(), cluster.getResolvedType());

            final SvChain existingChain = cluster.getChains().get(0);
            cluster.dissolveLinksAndChains();

            for(SvChain newChain : rearrangedChains)
            {
                newChain.setPloidyData(existingChain.ploidy(), existingChain.ploidyUncertainty());
                cluster.addChain(newChain, true);
            }
        }

        cluster.setResolved(isResolved, resolvedType);
    }

    private static void classifyInversionPairClusters(
            SvCluster cluster, boolean isFinal, int proximityThreshold, long longDupThreshold,
            SvBreakend startBe1, SvBreakend endBe1, SvBreakend startBe2, SvBreakend endBe2,
            boolean uniformPloidy, @NotNull final SvLinkedPair longestTiPair)
    {
        /* establish configuration:
            - Facing, no overlap
            - Facing, inner breakends overlap
            - Facing, outer breakends overlap
            - INV encloses INV
            - INV face away
         */

        // assign relative positions for easier comparison
        SvBreakend lowerBe1 = startBe1.position() < endBe1.position() ? startBe1 : endBe1;
        SvBreakend upperBe1 = lowerBe1 == startBe1 ? endBe1 : startBe1;
        SvBreakend lowerBe2 = startBe2.position() < endBe2.position() ? startBe2 : endBe2;
        SvBreakend upperBe2 = lowerBe2 == startBe2 ? endBe2 : startBe2;

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
        SvLinkedPair lowerDb1 = lowerBe1.getDBLink();
        SvLinkedPair upperDb1 = upperBe1.getDBLink();

        boolean lowerMatchingDB = lowerDb1 != null && lowerDb1 == lowerBe2.getDBLink();
        boolean upperMatchingDB = upperDb1 != null && upperDb1 == upperBe2.getDBLink();

        // basic reciprocal inversion
        if(lowerMatchingDB && upperMatchingDB)
        {
            if(!isFinal && (lowerDb1.length() > proximityThreshold || upperDb1.length() > proximityThreshold))
                return;

            cluster.setResolved(true, RECIP_INV);
            return;
        }

        boolean longOrLohBoundedTi = isLohBoundedTi(longestTiPair) || longestTiPair.length() > longDupThreshold;

        if((lowerBe1.position() < lowerBe2.position() && upperBe1.position() > upperBe2.position())
        || (lowerBe2.position() < lowerBe1.position() && upperBe2.position() > upperBe1.position()))
        {
            // one INV encloses the other - the DEL scenario
            if(longOrLohBoundedTi && uniformPloidy)
            {
                cluster.setResolved(false, DEL_TI);
            }
            else
            {
                long innerInversionLength = min(upperBe1.position() - lowerBe1.position(), upperBe2.position() - lowerBe2.position());

                if (innerInversionLength < 100000)
                {
                    cluster.setResolved(false, RESOLVED_FOLDBACK);
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

                    chain.reverseSectionOnBreakend(breakendToSwitch);

                    LOGGER.debug("cluster({}) reconfigured chain:", cluster.id());
                    chain.logLinks();

                    cluster.setResolved(false, RECIP_INV_DEL_DUP);
                }
            }
        }
        else
        {
            // otherwise the inner breakends overlap
            if(longOrLohBoundedTi && uniformPloidy)
            {
                cluster.setResolved(false, DUP_TI);
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

                chain.reverseSectionOnBreakend(breakendToSwitch);

                LOGGER.debug("cluster({}) reconfigured chain:", cluster.id());
                chain.logLinks();

                cluster.setResolved(false, RECIP_INV_DUPS);
            }
        }

        LOGGER.debug("cluster({}) longestTI({}) marked as {}",
                cluster.id(), longestTiPair.length(), cluster.getResolvedType());
    }

    public static void classifyDelDupPairClusters(
            SvCluster cluster,  SvBreakend startBe1, SvBreakend endBe1, SvBreakend startBe2, SvBreakend endBe2, boolean uniformPloidy)
    {
        if(!uniformPloidy)
            return;

        SvBreakend lowerBe1 = startBe1.position() < endBe1.position() ? startBe1 : endBe1;
        SvBreakend upperBe1 = lowerBe1 == startBe1 ? endBe1 : startBe1;
        SvBreakend lowerBe2 = startBe2.position() < endBe2.position() ? startBe2 : endBe2;
        SvBreakend upperBe2 = lowerBe2 == startBe2 ? endBe2 : startBe2;

        if(upperBe1.position() < lowerBe2.position() || upperBe2.position() < lowerBe1.position())
        {
            // ignore a DUP and DEL which don't overlap
            return;
        }

        if(lowerBe1.position() < lowerBe2.position() && upperBe1.position() > upperBe2.position())
        {
            // check the DEL isn't enclosing the DUP
            if(lowerBe1.orientation() == 1)
                return;

            cluster.setResolved(false, DUP_TI);
        }
        else if(lowerBe2.position() < lowerBe1.position() && upperBe2.position() > upperBe1.position())
        {
            if(lowerBe2.orientation() == 1)
                return;

            cluster.setResolved(false, DUP_TI);
        }
        else
        {
            // overlapping
            cluster.setResolved(false, DEL_TI);
        }
    }

}
