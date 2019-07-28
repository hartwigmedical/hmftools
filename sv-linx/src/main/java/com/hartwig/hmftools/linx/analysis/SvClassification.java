package com.hartwig.hmftools.linx.analysis;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.round;

import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.BND;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.DEL;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.DUP;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.INS;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.INV;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.SGL;
import static com.hartwig.hmftools.linx.analysis.ClusterAnalyser.SMALL_CLUSTER_SIZE;
import static com.hartwig.hmftools.linx.analysis.SvUtilities.copyNumbersEqual;
import static com.hartwig.hmftools.linx.chaining.LinkFinder.getMinTemplatedInsertionLength;
import static com.hartwig.hmftools.linx.types.ResolvedType.COMPLEX;
import static com.hartwig.hmftools.linx.types.ResolvedType.DUP_BE;
import static com.hartwig.hmftools.linx.types.ResolvedType.FB_INV_PAIR;
import static com.hartwig.hmftools.linx.types.ResolvedType.INF;
import static com.hartwig.hmftools.linx.types.ResolvedType.LINE;
import static com.hartwig.hmftools.linx.types.ResolvedType.LOW_VAF;
import static com.hartwig.hmftools.linx.types.ResolvedType.NONE;
import static com.hartwig.hmftools.linx.types.ResolvedType.PAIR_OTHER;
import static com.hartwig.hmftools.linx.types.ResolvedType.RECIP_DUPS;
import static com.hartwig.hmftools.linx.types.ResolvedType.RECIP_DUP_DEL;
import static com.hartwig.hmftools.linx.types.ResolvedType.RECIP_INV;
import static com.hartwig.hmftools.linx.types.ResolvedType.RECIP_TRANS;
import static com.hartwig.hmftools.linx.types.ResolvedType.SGL_PAIR_DEL;
import static com.hartwig.hmftools.linx.types.ResolvedType.SGL_PAIR_DUP;
import static com.hartwig.hmftools.linx.types.ResolvedType.SGL_PAIR_INS;
import static com.hartwig.hmftools.linx.types.ResolvedType.SIMPLE_GRP;
import static com.hartwig.hmftools.linx.types.ResolvedType.UNBAL_TRANS;
import static com.hartwig.hmftools.linx.types.SvaConstants.MIN_DEL_LENGTH;
import static com.hartwig.hmftools.linx.types.SvaConstants.SHORT_TI_LENGTH;

import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.variant.structural.StructuralVariantType;
import com.hartwig.hmftools.linx.types.ResolvedType;
import com.hartwig.hmftools.linx.types.SvBreakend;
import com.hartwig.hmftools.linx.chaining.SvChain;
import com.hartwig.hmftools.linx.types.SvCluster;
import com.hartwig.hmftools.linx.types.SvLinkedPair;
import com.hartwig.hmftools.linx.types.SvVarData;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class SvClassification
{
    // super category for an SV or cluster
    public static final String SUPER_TYPE_SIMPLE = "SIMPLE";
    public static final String SUPER_TYPE_INSERTION = "INSERTION";
    public static final String SUPER_TYPE_BREAK_PAIR = "BREAK_PAIR";
    public static final String SUPER_TYPE_COMPLEX = "COMPLEX";
    public static final String SUPER_TYPE_INCOMPLETE = "INCOMPLETE";
    public static final String SUPER_TYPE_ARTIFACT = "ARTIFACT";

    private static final Logger LOGGER = LogManager.getLogger(SvClassification.class);


    public static String getSuperType(SvCluster cluster)
    {
        ResolvedType resolvedType = cluster.getResolvedType();

        if (resolvedType == LINE)
            return SUPER_TYPE_INSERTION;

        if (resolvedType.isSimple())
            return SUPER_TYPE_SIMPLE;

        if(isFilteredResolvedType(resolvedType))
            return SUPER_TYPE_ARTIFACT;

        if(resolvedType == FB_INV_PAIR || resolvedType == RECIP_INV
        || resolvedType == RECIP_TRANS || resolvedType == UNBAL_TRANS
        || resolvedType == RECIP_DUPS || resolvedType == RECIP_DUP_DEL)
        {
            return SUPER_TYPE_BREAK_PAIR;
        }

        if(isIncompleteType(resolvedType) || resolvedType == PAIR_OTHER)
        {
            return SUPER_TYPE_INCOMPLETE;
        }

        return SUPER_TYPE_COMPLEX;
    }

    public static boolean isSimpleSingleSV(final SvCluster cluster)
    {
        return cluster.getSvCount() == 1 && cluster.getSV(0).isSimpleType();
    }

    public static boolean isIncompleteType(final ResolvedType resolvedType)
    {
        return (resolvedType == ResolvedType.INV || resolvedType == ResolvedType.SGL || resolvedType == INF);
    }

    public static boolean isSyntheticType(SvCluster cluster)
    {
        ResolvedType resolvedType = cluster.getResolvedType();

        if(resolvedType.isSimple())
        {
            if(cluster.getTypeCount(SGL) == 2)
                return false;
            else
                return cluster.getSvCount() > 1;
        }

        // will be more when look for synthetic INVs, SGLs and translocations
        if(resolvedType == RECIP_INV || resolvedType == RECIP_TRANS
        || resolvedType == RECIP_DUPS || resolvedType == RECIP_DUP_DEL
        || resolvedType == FB_INV_PAIR)
        {
            return cluster.getSvCount() > 2;
        }

        if(resolvedType == UNBAL_TRANS)
            return cluster.getSvCount() > 1;

        if(isIncompleteType(resolvedType))
            return cluster.getSvCount() > 1;

        return false;
    }

    public static boolean isFilteredResolvedType(final ResolvedType resolvedType)
    {
        return resolvedType.equals(DUP_BE) || resolvedType.equals(LOW_VAF);
    }

    public static void setClusterResolvedState(
            SvCluster cluster, boolean isFinal, long longDelThreshold, long longDupThreshold, int proximityThreshold)
    {
        if(cluster.getResolvedType() != NONE)
            return;

        if(cluster.hasLinkingLineElements())
        {
            // skip further classification for now
            cluster.setResolved(true, LINE);
            return;
        }

        if (smallSimpleSVs(cluster))
        {
            if(cluster.getSvCount() == 2 && (cluster.getTypeCount(DEL) + cluster.getTypeCount(DUP) == 2) && cluster.getTypeCount(DEL) != 2)
            {
                markSyntheticDelDups(cluster, longDelThreshold, longDupThreshold);

                if(cluster.getResolvedType() != NONE)
                    return;
            }

            boolean hasLongSVs = false;
            for(SvVarData var : cluster.getSVs())
            {
                if((var.type() == DEL && var.length() >= longDelThreshold) || (var.type() == DUP && var.length() >= longDupThreshold))
                {
                    hasLongSVs = true;
                    break;
                }
            }

            if(cluster.getSvCount() == 1)
            {
                StructuralVariantType type = cluster.getSV(0).type();

                if(type == DEL)
                    cluster.setResolved(!hasLongSVs, ResolvedType.DEL);
                else if(type == DUP)
                    cluster.setResolved(!hasLongSVs, ResolvedType.DUP);
                else if(type == INS)
                    cluster.setResolved(!hasLongSVs, ResolvedType.INS);
            }
            else
            {
                cluster.setResolved(!hasLongSVs, SIMPLE_GRP);
            }

            return;
        }

        markSyntheticTypes(cluster, longDelThreshold, longDupThreshold, proximityThreshold);

        if(cluster.getResolvedType() != NONE)
            return;

        if(isFinal)
        {
            if(cluster.getSvCount() == 1)
            {
                StructuralVariantType type = cluster.getSV(0).type();

                if(type == BND)
                    cluster.setResolved(false, UNBAL_TRANS);
                else if(type == INV)
                    cluster.setResolved(false, ResolvedType.INV);
                else if(type == SGL && cluster.getSV(0).isNoneSegment())
                    cluster.setResolved(false, INF);
                else if(type == SGL)
                    cluster.setResolved(false, ResolvedType.SGL);

                return;
            }

            markSyntheticIncompletes(cluster);

            if(cluster.getResolvedType() != NONE)
                return;

            if(cluster.getSvCount() == 2)
            {
                cluster.setResolved(false, PAIR_OTHER);
                return;
            }

            cluster.setResolved(false, COMPLEX);
        }
    }

    private static boolean smallSimpleSVs(final SvCluster cluster)
    {
        if(cluster.getSvCount() > SMALL_CLUSTER_SIZE)
            return false;

        for(final SvVarData var : cluster.getSVs())
        {
            if(!var.isSimpleType())
                return false;
        }

        return true;
    }

    public static void markSyntheticTypes(SvCluster cluster, long longDelThreshold, long longDupThreshold, int proximityThreshold)
    {
        if(cluster.getTypeCount(SGL) == 0)
        {
            markSyntheticReciprocalInversion(cluster, proximityThreshold);

            if (cluster.getResolvedType() != NONE)
                return;

            markSyntheticReciprocalTranslocation(cluster, proximityThreshold);

            if (cluster.getResolvedType() != NONE)
                return;

            markSyntheticDelDups(cluster, longDelThreshold, longDupThreshold);
        }
        else if (cluster.getSvCount() == 2 && cluster.isConsistent() && cluster.getTypeCount(SGL) == 2)
        {
            final SvVarData sgl1 = cluster.getSV(0);
            final SvVarData sgl2 = cluster.getSV(1);
            ResolvedType resolvedType = markSinglePairResolvedType(sgl1, sgl2);

            if(resolvedType != NONE)
            {
                long length = abs(sgl1.position(true) - sgl2.position(true));
                cluster.setResolved(true, resolvedType);
                cluster.setSyntheticData(length, 0);
            }
        }
    }

    public static void markSyntheticDelDups(SvCluster cluster, long longDelThreshold, long longDupThreshold)
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

        int totalChainLength = 0;
        int longTICount = 0;
        long longestTILength = 0;
        SvLinkedPair longestPair = null;
        for(SvLinkedPair pair : chain.getLinkedPairs())
        {
            if(pair.length() > SHORT_TI_LENGTH)
            {
                ++longTICount;

                if(longTICount > 1)
                    return;
            }

            if(pair.length() > longestTILength)
            {
                longestTILength = pair.length();
                longestPair = pair;
            }

            totalChainLength += pair.length();
        }

        long avgTiLength = round(totalChainLength / (double)chain.getLinkCount());
        long syntheticLength = abs(startBreakend.position() - endBreakend.position());

        ResolvedType resolvedType = NONE;

        if(longTICount == 0)
        {
            resolvedType = faceAway ? ResolvedType.DEL : ResolvedType.DUP;
        }
        else if(longTICount == 1)
        {
            // skip DEL-DUP constituents due to likely being false positives
            if(cluster.getSvCount() == 2 && cluster.getTypeCount(DEL) + cluster.getTypeCount(DUP) == 2)
                return;

            if(!faceAway && syntheticLength > SHORT_TI_LENGTH)
            {
                resolvedType = RECIP_DUPS;
            }
            else if(faceAway)
            {
                resolvedType = RECIP_DUP_DEL;
            }
        }

        if(resolvedType == NONE)
            return;

        LOGGER.debug("cluster({}) chain(links=({} len={} tiLen(longest={} avg={}) synLen({}) marked as {}",
                cluster.id(), chain.getLinkCount(), totalChainLength, longestTILength, avgTiLength, syntheticLength, resolvedType);

        boolean withinLongThreshold = false;

        if(resolvedType == ResolvedType.DEL)
            withinLongThreshold = syntheticLength < longDelThreshold;
        else if(resolvedType == ResolvedType.DUP || resolvedType == RECIP_DUPS)
            withinLongThreshold = syntheticLength < longDupThreshold;
        else
            withinLongThreshold = syntheticLength < max(longDelThreshold,longDupThreshold);

        boolean resolved = withinLongThreshold && (longTICount == 0);
        cluster.setResolved(resolved, resolvedType);
        cluster.setSyntheticData(syntheticLength, longestTILength);

        if(resolvedType == RECIP_DUPS)
        {
            classifyReciprocalDupPairs(cluster, chain, longestPair);
        }
    }

    private static void classifyReciprocalDupPairs(SvCluster cluster, SvChain chain, SvLinkedPair longestPair)
    {
        SvBreakend startBreakend = chain.getOpenBreakend(true);
        SvBreakend endBreakend = chain.getOpenBreakend(false);

        long longestTILength = 0;
        long syntheticLength = 0;

        if(!longestPair.chromosome().equals(startBreakend.chromosome()))
        {
            // remove the long chain link on the assumption that these breakends aren't joined,
            // leaving 2 separate chains or SVs
            cluster.getChains().clear();

            SvChain newChain = new SvChain(0);

            for (SvLinkedPair pair : chain.getLinkedPairs())
            {
                if (newChain == null)
                    newChain = new SvChain(1);

                if (pair.length() <= SHORT_TI_LENGTH)
                {
                    newChain.addLink(pair, false);
                    longestTILength = max(pair.length(), longestTILength);
                }
                else
                {
                    // skip the long TI and cache the chain
                    if (newChain != null && newChain.getLinkCount() > 0)
                    {
                        cluster.addChain(newChain, false);
                        newChain = null;
                    }
                }
            }

            if (newChain != null && newChain.getLinkCount() > 0)
                cluster.addChain(newChain, false);

            LOGGER.debug("cluster({}) split reciprocal DUP pair into 2 chains between chrs({} & {})",
                    cluster.id(), longestPair.chromosome(), startBreakend.chromosome());
        }
        else
        {
            SvBreakend chainLowerBe = startBreakend.position() < endBreakend.position() ? startBreakend : endBreakend;
            SvBreakend chainUpperBe = startBreakend == chainLowerBe ? endBreakend : startBreakend;

            SvBreakend lowerTiBe = longestPair.getBreakend(true);
            SvBreakend upperTiBe = longestPair.getBreakend(false);

            if(lowerTiBe.position() > chainLowerBe.position() && upperTiBe.position() < chainUpperBe.position())
            {
                cluster.setResolved(false, FB_INV_PAIR);
                return;
            }

            // change the chaining configuration to make a longer TI between the opposing breakends

            // artificially create the new TI
            SvLinkedPair newLink = null;

            if(chainLowerBe.orientation() == -1 && upperTiBe.orientation() == 1 && chainLowerBe.position() < upperTiBe.position()
            && chainLowerBe.getSV() != upperTiBe.getSV())
            {
                newLink = SvLinkedPair.from(chainLowerBe, upperTiBe);
            }
            else if(chainUpperBe.orientation() == 1 && lowerTiBe.orientation() == -1 && chainUpperBe.position() > lowerTiBe.position()
            && chainUpperBe.getSV() != lowerTiBe.getSV())
            {
                newLink = SvLinkedPair.from(lowerTiBe, chainUpperBe);
            }
            else
            {
                LOGGER.debug("cluster({}) unresolvable RECIP DUP PAIR", cluster.id());
                cluster.setResolved(false, COMPLEX);
                cluster.setSyntheticData(0, 0);
                return;
            }

            SvChain newChain = new SvChain(0);
            newChain.addLink(newLink, true);

            longestTILength = newLink.length();

            List<SvLinkedPair> linksToAdd = chain.getLinkedPairs().stream()
                    .filter(x -> x.length() <= SHORT_TI_LENGTH)
                    .collect(Collectors.toList());

            int index = 0;
            int iterations = 0;
            int linksCount = linksToAdd.size();
            while(!linksToAdd.isEmpty())
            {
                if(index >= linksToAdd.size())
                    index = 0;

                SvLinkedPair pair = linksToAdd.get(index);

                longestTILength = max(pair.length(), longestTILength);

                if(newChain.canAddLinkedPairToStart(pair))
                {
                    newChain.addLink(pair, true);
                    linksToAdd.remove(index);
                }
                else if(newChain.canAddLinkedPairToEnd(pair))
                {
                    newChain.addLink(pair, false);
                    linksToAdd.remove(index);
                }
                else
                {
                    ++index;
                }

                ++iterations;

                if(iterations > 2 * linksCount)
                    break;
            }

            if(!linksToAdd.isEmpty())
            {
                LOGGER.warn("cluster({}) failed to build new RECIP DUP PAIR chain from links({})", cluster.id(), linksCount);
                cluster.setResolved(false, COMPLEX);
                cluster.setSyntheticData(0, 0);
                return;
            }

            cluster.getChains().clear();

            if (newChain != null && newChain.getLinkCount() > 0)
                cluster.addChain(newChain, false);

            syntheticLength = abs(newChain.getOpenBreakend(true).position() - newChain.getOpenBreakend(false).position());

            LOGGER.debug("cluster({}) reconfigured reciprocal DUP pair with newLink({})", cluster.id(), newLink.length());
        }

        cluster.setResolved(false, RECIP_DUPS);
        cluster.setSyntheticData(syntheticLength, longestTILength);
    }

    public static void markSyntheticReciprocalInversion(SvCluster cluster, int proximityThreshold)
    {
        if(!cluster.isFullyChained(false) || cluster.getChains().size() != 1)
            return;

        SvChain chain = cluster.getChains().get(0);

        int totalLinks = 0;

        // check chains are both short and have the correct orientations
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

        SvBreakend startBe1 = tiPair.getBreakend(true);
        SvBreakend endBe1 = tiPair.getBreakend(false);

        SvBreakend startBe2 = chain.getOpenBreakend(true);
        SvBreakend endBe2 = chain.getOpenBreakend(false);

        // check for same arm and opposite orientations
        if(!startBe1.getChrArm().equals(startBe2.getChrArm()))
            return;

        if(!startBe2.getChrArm().equals(endBe2.getChrArm()) || startBe2.orientation() == endBe2.orientation())
            return;

        // check for linked, short DBs at both ends
        SvLinkedPair startDB = startBe1.getSV().getDBLink(startBe1.usesStart());
        SvLinkedPair endDB = endBe1.getSV().getDBLink(endBe1.usesStart());

        if(startDB == null || startDB.length() > proximityThreshold || endDB == null || endDB.length() > proximityThreshold)
            return;

        long syntheticLength = 0;
        if ((startDB.hasBreakend(startBe2) && endDB.hasBreakend(endBe2)) || (startDB.hasBreakend(endBe2) && endDB.hasBreakend(startBe2)))
        {
            syntheticLength = abs(startBe2.position() - endBe2.position());
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
        cluster.setSyntheticData(syntheticLength, tiPair.length());
    }

    public static void markSyntheticReciprocalTranslocation(SvCluster cluster, int proximityThreshold)
    {
        // can be formed from 1 BNDs, 1 chain and a BND, or 2 chains

        if(cluster.getChains().size() > 2)
            return;

        SvBreakend startBe1 = null;
        SvBreakend endBe1 = null;
        SvBreakend startBe2 = null;
        SvBreakend endBe2 = null;

        long longestTILength = 0;
        int totalLinks = 0;

        if(cluster.getChains().isEmpty())
        {
            if(cluster.getSvCount() != 2 || cluster.getTypeCount(BND) != 2)
                return;

            final SvVarData var1 = cluster.getSV(0);
            final SvVarData var2 = cluster.getSV(1);

            startBe1 = var1.getBreakend(true);
            endBe1 = var1.getBreakend(false);
            startBe2 = var2.getBreakend(true);
            endBe2 = var2.getBreakend(false);

        }
        else if(cluster.isFullyChained(false) && cluster.getChains().size() == 2)
        {
            SvChain chain1 = cluster.getChains().get(0);
            SvChain chain2 = cluster.getChains().get(1);

            // check chains are both short and have the correct orientations
            for(SvChain chain : cluster.getChains())
            {
                for (SvLinkedPair pair : chain.getLinkedPairs())
                {
                    if (pair.length() > SHORT_TI_LENGTH)
                        return;

                    longestTILength = max(pair.length(), longestTILength);
                    ++totalLinks;
                }
            }

            startBe1 = chain1.getOpenBreakend(true);
            endBe1 = chain1.getOpenBreakend(false);
            startBe2 = chain2.getOpenBreakend(true);
            endBe2 = chain2.getOpenBreakend(false);
        }
        else if(cluster.getChains().size() == 1 && cluster.getUnlinkedSVs().size() == 1)
        {
            // check for a single SV and a chain
            SvVarData var = cluster.getUnlinkedSVs().get(0);
            if(var.type() != BND)
                return;

            SvChain chain = cluster.getChains().get(0);

            for (SvLinkedPair pair : chain.getLinkedPairs())
            {
                if (pair.length() > SHORT_TI_LENGTH)
                    return;

                longestTILength = max(pair.length(), longestTILength);
                ++totalLinks;
            }

            startBe1 = chain.getOpenBreakend(true);
            endBe1 = chain.getOpenBreakend(false);

            startBe2 = var.getBreakend(true);
            endBe2 = var.getBreakend(false);
        }
        else
        {
            return;
        }

        if(startBe1.getChrArm().equals(endBe1.getChrArm()) || startBe2.getChrArm().equals(endBe2.getChrArm()))
            return;

        // check for linked, short DBs at both ends
        SvLinkedPair startDB = startBe1.getSV().getDBLink(startBe1.usesStart());
        SvLinkedPair endDB = endBe1.getSV().getDBLink(endBe1.usesStart());

        if(startDB == null || startDB.length() > proximityThreshold || endDB == null || endDB.length() > proximityThreshold)
            return;

        if(startBe1.getChrArm().equals(startBe2.getChrArm()) && startBe1.orientation() != startBe2.orientation()
                && endBe1.getChrArm().equals(endBe2.getChrArm()) && endBe1.orientation() != endBe2.orientation())
        {
            // start to start
            if (!startDB.hasBreakend(startBe2) || !endDB.hasBreakend(endBe2))
                return;

            // if(arePairedDeletionBridges(var1, var2))

        }
        else if(startBe1.getChrArm().equals(endBe2.getChrArm()) && startBe1.orientation() != endBe2.orientation()
                && endBe1.getChrArm().equals(startBe2.getChrArm()) && endBe1.orientation() != startBe2.orientation())
        {
            // start to end
            if (!startDB.hasBreakend(endBe2) || !endDB.hasBreakend(startBe2))
                return;
        }
        else
        {
            return;
        }

        ResolvedType resolvedType = RECIP_TRANS;

        LOGGER.debug("cluster({}) chain(links=({} longestTI={}) marked as {}",
                cluster.id(), totalLinks, longestTILength, resolvedType);

        cluster.setResolved(true, resolvedType);
        cluster.setSyntheticData(0, longestTILength);
    }

    public static void markSyntheticIncompletes(SvCluster cluster)
    {
        // look for chains of short TIs which when reduced form a SGL, INF INV or unbalanced TRANS (ie a BND)
        if (cluster.getSvCount() > 5 || !cluster.isFullyChained(false) || cluster.getChains().size() != 1)
            return;

        SvChain chain = cluster.getChains().get(0);

        // test the chain for short TIs only
        int totalChainLength = 0;
        int longTICount = 0;
        long longestTILength = 0;
        for (SvLinkedPair pair : chain.getLinkedPairs())
        {
            if (pair.length() > SHORT_TI_LENGTH)
                return;

            longestTILength = max(pair.length(), longestTILength);
            totalChainLength += pair.length();
        }

        // first look for chains ending in SGLs or INFs
        ResolvedType resolvedType = NONE;

        if((chain.getLastSV().type() == SGL && chain.getLastSV().isNoneSegment())
        || (chain.getFirstSV().type() == SGL && chain.getFirstSV().isNoneSegment()))
        {
            resolvedType = INF;
        }
        else if(chain.getLastSV().type() == SGL || chain.getFirstSV().type() == SGL)
        {
            resolvedType = ResolvedType.SGL;
        }
        else
        {
            final SvBreakend startBreakend = chain.getOpenBreakend(true);
            final SvBreakend endBreakend = chain.getOpenBreakend(false);

            if (!startBreakend.chromosome().equals(endBreakend.chromosome()) || startBreakend.arm() != endBreakend.arm())
            {
                resolvedType = UNBAL_TRANS;
            }
            else
            {
                resolvedType = ResolvedType.INV;
            }
        }

        if (resolvedType == NONE)
            return;

        LOGGER.debug("cluster({}) chain(links=({} len={} tiLen({}) marked as {}",
                cluster.id(), chain.getLinkCount(), totalChainLength, longestTILength, resolvedType);

        cluster.setResolved(false, resolvedType);
        cluster.setSyntheticData(0, longestTILength);
    }

    public static ResolvedType markSinglePairResolvedType(final SvVarData sgl1, final SvVarData sgl2)
    {
        if(sgl1.sglToCentromereOrTelomere() || sgl2.sglToCentromereOrTelomere())
            return NONE;

        final SvBreakend breakend1 = sgl1.getBreakend(true);
        final SvBreakend breakend2 = sgl2.getBreakend(true);

        // to form a simple del or dup, they need to have different orientations
        if(breakend1.orientation() == breakend2.orientation())
            return NONE;

        // check copy number consistency
        double cn1 = sgl2.copyNumberChange(true);
        double cn2 = sgl1.copyNumberChange(true);

        if(!copyNumbersEqual(cn1, cn2))
            return NONE;

        boolean breakendsFace = (breakend1.position() < breakend2.position() && breakend1.orientation() == -1)
                || (breakend2.position() < breakend1.position() && breakend2.orientation() == -1);

        long length = abs(breakend1.position() - breakend2.position());

        ResolvedType resolvedType = NONE;

        if(breakendsFace)
        {
            // a DUP if breakends are further than the anchor distance away, else an INS
            int minTiLength = getMinTemplatedInsertionLength(breakend1, breakend2);
            if(length >= minTiLength)
                resolvedType = SGL_PAIR_DUP;
            else
                resolvedType = SGL_PAIR_INS;
        }
        else
        {
            // a DEL if the breakends are further than the min DEL length, else an INS
            if(length >= MIN_DEL_LENGTH)
                resolvedType = SGL_PAIR_DEL;
            else
                resolvedType = SGL_PAIR_INS;
        }

        // mark these differently from those formed from normal SVs
        return resolvedType;
    }


}
