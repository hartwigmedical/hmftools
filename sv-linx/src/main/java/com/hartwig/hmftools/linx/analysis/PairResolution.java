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
import static com.hartwig.hmftools.linx.chaining.ChainPloidyLimits.calcPloidyUncertainty;
import static com.hartwig.hmftools.linx.chaining.ChainPloidyLimits.ploidyMatch;
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
        // handles synthetic DELs and DUPs from a single chain with no long TIs, and other paired configurations from either 2 SVs,
        // a chain and an SV or 2 chains
        if(cluster.getChains().size() > 2 || cluster.getTypeCount(SGL) > 0)
            return;

        // establish the nature of the breakends
        // first reduce BNDs and/or chains to a set of breakends and record TI length(s)
        SvBreakend startBe1 = null;
        SvBreakend endBe1 = null;
        SvBreakend startBe2 = null;
        SvBreakend endBe2 = null;

        SvLinkedPair longTiLink = null;
        boolean uniformPloidy = false;

        if(!cluster.getChains().isEmpty())
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

        // first handle the single chain case either handling it as a synthetic with only short TIs, or by splitting it
        if(cluster.isFullyChained(false) && cluster.getChains().size() == 1)
        {
            if(longTiLink == null)
            {
                classifySyntheticDelDups(cluster, longDelThreshold, longDupThreshold);

                // may alternatively be another synthetic type, so don't classify as PAIR OTHER
                return;
            }
            else
            {
                // go into regular 2 break logic by effectively breaking the chain at this long TI
                SvChain chain = cluster.getChains().get(0);

                cluster.dissolveLinksAndChains();

                final List<SvLinkedPair> pairs = chain.getLinkedPairs();

                if(pairs.size() > 1)
                {
                    SvChain newChain = new SvChain(chain.id());
                    newChain.setPloidyData(chain.ploidy(), chain.ploidyUncertainty());

                    if (pairs.get(0) == longTiLink || pairs.get(pairs.size() - 1) == longTiLink)
                    {
                        for (SvLinkedPair pair : pairs)
                        {
                            if (pair != longTiLink)
                            {
                                newChain.addLink(pair, false);
                            }
                        }

                        cluster.addChain(newChain, false);
                    }
                    else
                    {
                        for (SvLinkedPair pair : pairs)
                        {
                            if (pair == longTiLink)
                            {
                                cluster.addChain(newChain, false);

                                newChain = new SvChain(chain.id() + 1);
                                newChain.setPloidyData(chain.ploidy(), chain.ploidyUncertainty());
                                continue;
                            }

                            newChain.addLink(pair, false);
                        }

                        cluster.addChain(newChain, false);
                    }
                }
            }
        }

        if(cluster.getChains().isEmpty())
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
        else if(cluster.isFullyChained(false) && cluster.getChains().size() == 2)
        {
            SvChain chain1 = cluster.getChains().get(0);
            SvChain chain2 = cluster.getChains().get(1);

            startBe1 = chain1.getOpenBreakend(true);
            endBe1 = chain1.getOpenBreakend(false);
            startBe2 = chain2.getOpenBreakend(true);
            endBe2 = chain2.getOpenBreakend(false);
            uniformPloidy = ploidyMatch(chain1.ploidy(), chain1.ploidyUncertainty(), chain2.ploidy(), chain2.ploidyUncertainty());
        }
        else if(cluster.getChains().size() == 1 && cluster.getUnlinkedSVs().size() == 1)
        {
            // check for a single SV and a chain
            SvVarData var = cluster.getUnlinkedSVs().get(0);
            SvChain chain = cluster.getChains().get(0);

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
                    cluster, isFinal, proximityThreshold, longDupThreshold,
                    startBe1, endBe1, startBe2, endBe2, uniformPloidy, longTiLink);
        }

        // if(cluster.getResolvedType() == NONE)
        //    cluster.setResolved(false, PAIR_OTHER);
    }

    public static void classifySyntheticDelDups(SvCluster cluster, long longDelThreshold, long longDupThreshold)
    {
        if(!cluster.isFullyChained(true) || cluster.getChains().size() != 1 || cluster.getTypeCount(SGL) > 0)
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
        long longestTILength = chain.getLinkedPairs().stream().mapToLong(x -> x.length()).max().getAsLong();
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
            boolean uniformPloidy, final SvLinkedPair longestTiPair)
    {
        // first reduce BNDs and/or chains to a set of breakends and record TI length(s)
        if(startBe1.getChrArm().equals(endBe1.getChrArm()) || startBe2.getChrArm().equals(endBe2.getChrArm()))
            return;

        /* establish configuration:
            - Deletion bridge on both arms
            - Facing on one arm.  Deletion Bridge on other arm
            - Facing on 2 arms
            - 1 arm shared only. Facing or deletion bridge
         */

        // check for 3 different arms being involved
        List<String> arms = Lists.newArrayList(startBe1.getChrArm(), endBe1.getChrArm());
        if(!arms.contains(startBe2.getChrArm()) || !arms.contains(endBe2.getChrArm()))
        {
            return;
        }

        // assign to arms for easier comparison
        SvBreakend arm1Be1 = startBe1;
        SvBreakend arm2Be1 = endBe1;
        SvBreakend arm1Be2 = startBe2.getChrArm().equals(arm1Be1.getChrArm()) ? startBe2 : endBe2;
        SvBreakend arm2Be2 = startBe2.getChrArm().equals(arm2Be1.getChrArm()) ? startBe2 : endBe2;

        // check for linked, short DBs at both ends
        SvLinkedPair arm1Db1 = arm1Be1.getDBLink();
        SvLinkedPair arm2Db1 = arm2Be1.getDBLink();

        boolean arm1MatchingDB = arm1Db1 != null && arm1Db1 == arm1Be2.getDBLink();
        boolean arm2MatchingDB = arm2Db1 != null && arm2Db1 == arm2Be2.getDBLink();

        // basic reciprocal translocation
        if(arm1MatchingDB && arm2MatchingDB)
        {
            if(!isFinal && (arm1Db1.length() > proximityThreshold || arm2Db1.length() > proximityThreshold))
                return;

            cluster.setResolved(true, RECIP_TRANS);
            return;
        }

        if(longestTiPair == null)
            return;

        ResolvedType resolvedType;
        boolean longOrLohBoundedTi = isLohBoundedTi(longestTiPair) || longestTiPair.length() > longDupThreshold;

        if(!arm1MatchingDB && !arm2MatchingDB)
        {
            resolvedType = longOrLohBoundedTi && uniformPloidy ? DUP_TI : RECIP_TRANS_DUPS;
        }
        else
        {
            resolvedType = longOrLohBoundedTi && uniformPloidy ? DEL_TI : RECIP_TRANS_DEL_DUP;
        }

        LOGGER.debug("cluster({}) longestTI({}) marked as {}",
                cluster.id(), longestTiPair.length(), resolvedType);

        cluster.setResolved(false, resolvedType);
    }

    private static void classifyInversionPairClusters(
            SvCluster cluster, boolean isFinal, int proximityThreshold, long longDupThreshold,
            SvBreakend startBe1, SvBreakend endBe1, SvBreakend startBe2, SvBreakend endBe2,
            boolean uniformPloidy, final SvLinkedPair longestTiPair)
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

        if(longestTiPair == null)
            return;

        boolean longOrLohBoundedTi = isLohBoundedTi(longestTiPair) || longestTiPair.length() > longDupThreshold;

        if((lowerBe1.position() < lowerBe2.position() && upperBe1.position() > upperBe2.position())
        || (lowerBe2.position() < lowerBe1.position() && upperBe2.position() > upperBe1.position()))
        {
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
                    cluster.setResolved(false, RECIP_INV_DEL_DUP);
                }
            }
        }
        else
        {
            if(longOrLohBoundedTi && uniformPloidy)
            {
                cluster.setResolved(false, DUP_TI);
            }
            else
            {
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


    public static void markInversionPairClusters(SvCluster cluster, boolean isFinal, int proximityThreshold)
    {
        if(!cluster.isFullyChained(false) || cluster.getChains().size() != 1)
            return;

        SvChain chain = cluster.getChains().get(0);

        int totalLinks = 0;

        SvLinkedPair tiPair = null;

        for (SvLinkedPair pair : chain.getLinkedPairs())
        {
            if(tiPair == null)
            {
                tiPair = pair;
            }
            else if(pair.length() > tiPair.length())
            {
                if(tiPair.length() > SHORT_TI_LENGTH)
                    return; // can only be one long pair

                tiPair = pair;
            }


            ++totalLinks;
        }

        if(tiPair == null)
            return;

        SvBreakend tiStart = tiPair.getBreakend(true);
        SvBreakend tiEnd = tiPair.getBreakend(false);

        SvBreakend chainStart = chain.getOpenBreakend(true);
        SvBreakend chainEnd = chain.getOpenBreakend(false);

        // check for opposite orientations

        if(chainStart.orientation() == chainEnd.orientation())
            return;

        // check for linked, short DBs at both ends
        SvLinkedPair startDB = tiStart.getSV().getDBLink(tiStart.usesStart());
        SvLinkedPair endDB = tiEnd.getSV().getDBLink(tiEnd.usesStart());

        if(startDB == null || endDB == null)
            return;

        if(!isFinal && (startDB.length() > proximityThreshold || endDB.length() > proximityThreshold))
            return;

        long syntheticLength = 0;
        if ((startDB.hasBreakend(chainStart) && endDB.hasBreakend(chainEnd)) || (startDB.hasBreakend(chainEnd) && endDB.hasBreakend(chainStart)))
        {
            syntheticLength = abs(chainStart.position() - chainEnd.position());
        }
        else
        {
            return;
        }

        // require the large TI to be at least 50% of the length of teh
        if (tiPair.length() < 0.5 * syntheticLength)
            return;

        ResolvedType resolvedType = RECIP_INV;

        LOGGER.debug("cluster({}) chain(links=({} longestTI={}) synLen({}) marked as {}",
                cluster.id(), totalLinks, tiPair.length(), syntheticLength, resolvedType);

        cluster.setResolved(false, resolvedType);
    }

}
