package com.hartwig.hmftools.svanalysis.analysis;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.round;

import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.BND;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.DEL;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.DUP;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.INS;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.INV;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.SGL;
import static com.hartwig.hmftools.svanalysis.analysis.ClusterAnalyser.SMALL_CLUSTER_SIZE;
import static com.hartwig.hmftools.svanalysis.analysis.SvClusteringMethods.markSinglePairResolvedType;
import static com.hartwig.hmftools.svanalysis.types.SvaConstants.SHORT_TI_LENGTH;

import com.hartwig.hmftools.common.variant.structural.StructuralVariant;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantType;
import com.hartwig.hmftools.svanalysis.types.SvBreakend;
import com.hartwig.hmftools.svanalysis.types.SvChain;
import com.hartwig.hmftools.svanalysis.types.SvCluster;
import com.hartwig.hmftools.svanalysis.types.SvLinkedPair;
import com.hartwig.hmftools.svanalysis.types.SvVarData;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class SvClassification
{
    // Resolved types - all the standard usual SV types plus
    public static final String RESOLVED_TYPE_NONE = "NONE";
    public static final String RESOLVED_TYPE_DUP_BE = "DUP_BE";
    public static final String RESOLVED_TYPE_POLY_G_C = "POLY_G_C";
    public static final String RESOLVED_TYPE_LINE = "LINE";

    public static final String RESOLVED_TYPE_DEL = "DEL";
    public static final String RESOLVED_TYPE_DUP = "DUP";
    public static final String RESOLVED_TYPE_INF = "INF";
    public static final String RESOLVED_TYPE_SGL = "SGL";
    public static final String RESOLVED_TYPE_INV = "INV";
    public static final String RESOLVED_TYPE_INS = "INS";
    public static final String RESOLVED_TYPE_RECIPROCAL_TRANS = "RECIP_TRANS";
    public static final String RESOLVED_TYPE_RECIPROCAL_INV = "RECIP_INV";
    public static final String RESOLVED_TYPE_UNBALANCED_TRANS = "UNBAL_TRANS";
    public static final String RESOLVED_TYPE_RECIPROCAL_DUP_PAIR = "RECIP_DUPS";
    public static final String RESOLVED_TYPE_RECIPROCAL_DUP_DEL = "RECIP_DUP_DEL";
    public static final String RESOLVED_TYPE_COMPLEX = "COMPLEX";
    public static final String RESOLVED_TYPE_PAIR_OTHER = "PAIR_OTHER";
    public static final String RESOLVED_TYPE_SIMPLE_GRP = "SIMPLE_GRP";
    public static final String RESOLVED_TYPE_FB_INV_PAIR = "FB_INV_PAIR";

    // public static final String RESOLVED_TYPE_SIMPLE_SV = "SIMPLE";
    // public static final String RESOLVED_TYPE_SYNTH_DEL = "SYNTH_DEL";
    // public static final String RESOLVED_TYPE_SYNTH_DUP = "SYNTH_DUP";
    // public static final String RESOLVED_TYPE_SGL_PAIR_INS = "SGL_PAIR_INS";
    // public static final String RESOLVED_TYPE_SGL_PAIR_DEL = "SGL_PAIR_DEL";
    // public static final String RESOLVED_TYPE_SGL_PAIR_DUP = "SGL_PAIR_DUP";
    // public static final String RESOLVED_TYPE_SGL_PLUS_INCONSISTENT = "SGL_BND_INV";

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
        String resolvedType = cluster.getResolvedType();

        if (resolvedType == RESOLVED_TYPE_LINE)
            return SUPER_TYPE_INSERTION;

        if (isSimpleType(resolvedType))
            return SUPER_TYPE_SIMPLE;

        if(isFilteredResolvedType(resolvedType))
            return SUPER_TYPE_ARTIFACT;

        // should just leave synthetic types
        if(isSyntheticSimpleType(resolvedType))
            return SUPER_TYPE_SIMPLE;

        if(resolvedType == RESOLVED_TYPE_FB_INV_PAIR || resolvedType == RESOLVED_TYPE_RECIPROCAL_INV
        || resolvedType == RESOLVED_TYPE_RECIPROCAL_TRANS || resolvedType == RESOLVED_TYPE_UNBALANCED_TRANS
        || resolvedType == RESOLVED_TYPE_RECIPROCAL_DUP_PAIR || resolvedType == RESOLVED_TYPE_RECIPROCAL_DUP_DEL)
        {
            return SUPER_TYPE_BREAK_PAIR;
        }

        if(resolvedType == RESOLVED_TYPE_INV || resolvedType == RESOLVED_TYPE_SGL || resolvedType == RESOLVED_TYPE_INF
        || resolvedType == RESOLVED_TYPE_PAIR_OTHER)
        {
            return SUPER_TYPE_INCOMPLETE;
        }

        return SUPER_TYPE_COMPLEX;
    }

    public static boolean isSyntheticSimpleType(final String resolvedType)
    {
        return isSimpleType(resolvedType);
    }

    public static boolean isSimpleSingleSV(final SvCluster cluster)
    {
        return cluster.getSvCount() == 1 && isSimpleType(cluster.getResolvedType());
    }

    public static boolean isSimpleType(final String resolvedType)
    {
        if(resolvedType == RESOLVED_TYPE_DEL || resolvedType == RESOLVED_TYPE_DUP || resolvedType == RESOLVED_TYPE_INS)
        {
            return true;
        }

        return false;
    }

    public static boolean isSyntheticType(SvCluster cluster)
    {
        if(isSyntheticSimpleType(cluster.getResolvedType()))
            return true;

        // will be more when look for synthetic INVs, SGLs and translocations
        if(cluster.getResolvedType() == RESOLVED_TYPE_RECIPROCAL_INV || cluster.getResolvedType() == RESOLVED_TYPE_RECIPROCAL_TRANS)
            return cluster.getSvCount() > 2;

        return false;
    }

    public static boolean isFilteredResolvedType(final String resolvedType)
    {
        return resolvedType.equals(RESOLVED_TYPE_POLY_G_C) || resolvedType.equals(RESOLVED_TYPE_DUP_BE);
    }

    public static void setClusterResolvedState(
            SvCluster cluster, boolean isFinal, long delDupLongThreshold, int proximityThreshold)
    {
        if(cluster.getResolvedType() != RESOLVED_TYPE_NONE)
            return;

        if(cluster.hasLinkingLineElements())
        {
            // skip further classification for now
            cluster.setResolved(true, RESOLVED_TYPE_LINE);
            return;
        }

        if (smallSimpleSVs(cluster))
        {
            if(cluster.getSvCount() == 2 && (cluster.getTypeCount(DEL) + cluster.getTypeCount(DUP) == 2) && cluster.getTypeCount(DEL) != 2)
            {
                markSyntheticDelDups(cluster, delDupLongThreshold);

                if(cluster.getResolvedType() != RESOLVED_TYPE_NONE)
                    return;
            }

            boolean hasLongSVs = false;
            for(SvVarData var : cluster.getSVs())
            {
                if(var.length() >= delDupLongThreshold)
                {
                    hasLongSVs = true;
                    break;
                }
            }

            if(cluster.getSvCount() == 1)
            {
                StructuralVariantType type = cluster.getSV(0).type();

                if(type == DEL)
                    cluster.setResolved(!hasLongSVs, RESOLVED_TYPE_DEL);
                else if(type == DUP)
                    cluster.setResolved(!hasLongSVs, RESOLVED_TYPE_DUP);
                else if(type == INS)
                    cluster.setResolved(!hasLongSVs, RESOLVED_TYPE_INS);
            }
            else
            {
                cluster.setResolved(!hasLongSVs, RESOLVED_TYPE_SIMPLE_GRP);
            }

            return;
        }

        markSyntheticTypes(cluster, delDupLongThreshold, proximityThreshold);

        if(cluster.getResolvedType() != RESOLVED_TYPE_NONE)
            return;

        if(isFinal)
        {
            if(cluster.getSvCount() == 1)
            {
                StructuralVariantType type = cluster.getSV(0).type();

                if(type == BND)
                    cluster.setResolved(false, RESOLVED_TYPE_UNBALANCED_TRANS);
                else if(type == INV)
                    cluster.setResolved(false, RESOLVED_TYPE_INV);
                else if(type == SGL && cluster.getSV(0).isNoneSegment())
                    cluster.setResolved(false, RESOLVED_TYPE_INF);
                else if(type == SGL)
                    cluster.setResolved(false, RESOLVED_TYPE_SGL);

                return;
            }
            else if(cluster.getSvCount() == 2)
            {
                cluster.setResolved(false, RESOLVED_TYPE_PAIR_OTHER);
                return;
            }

            cluster.setResolved(false, RESOLVED_TYPE_COMPLEX);
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

    public static void markSyntheticTypes(SvCluster cluster, long delDupLongThreshold, int proximityThreshold)
    {
        if(cluster.getTypeCount(SGL) == 0)
        {
            markSyntheticReciprocalInversion(cluster, proximityThreshold);

            if (cluster.getResolvedType() != RESOLVED_TYPE_NONE)
                return;

            markSyntheticReciprocalTranslocation(cluster, proximityThreshold);

            if (cluster.getResolvedType() != RESOLVED_TYPE_NONE)
                return;

            markSyntheticDelDups(cluster, delDupLongThreshold);
        }
        else if (cluster.getSvCount() == 2 && cluster.isConsistent() && cluster.getTypeCount(SGL) == 2)
        {
            final SvVarData sgl1 = cluster.getSV(0);
            final SvVarData sgl2 = cluster.getSV(1);
            String resolvedType = markSinglePairResolvedType(sgl1, sgl2);

            if(resolvedType != RESOLVED_TYPE_NONE)
            {
                long length = abs(sgl1.position(true) - sgl2.position(true));
                cluster.setResolved(true, resolvedType);
                cluster.setSyntheticData(length, 0);
            }
        }
    }

    public static void markSyntheticDelDups(SvCluster cluster, long delDupLongThreshold)
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

        String resolvedType = "";

        if(longTICount == 0)
        {
            resolvedType = faceAway ? RESOLVED_TYPE_DEL : RESOLVED_TYPE_DUP;
        }
        else if(longTICount == 1)
        {
            // skip DEL-DUP constituents due to likely being false positives
            if(cluster.getSvCount() == 2 && cluster.getTypeCount(DEL) + cluster.getTypeCount(DUP) == 2)
                return;

            if(!faceAway && syntheticLength > SHORT_TI_LENGTH)
            {
                resolvedType = RESOLVED_TYPE_RECIPROCAL_DUP_PAIR;
            }
            else if(faceAway)
            {
                resolvedType = RESOLVED_TYPE_RECIPROCAL_DUP_DEL;
            }
        }

        if(resolvedType.isEmpty())
            return;

        LOGGER.debug("cluster({}) chain(links=({} len={} tiLen(longest={} avg={}) synLen({}) marked as {}",
                cluster.id(), chain.getLinkCount(), totalChainLength, longestTILength, avgTiLength, syntheticLength, resolvedType);

        boolean setResolved = (syntheticLength < delDupLongThreshold) && (longTICount == 0);
        cluster.setResolved(setResolved, resolvedType);
        cluster.setSyntheticData(syntheticLength, longestTILength);

        if(resolvedType == RESOLVED_TYPE_RECIPROCAL_DUP_PAIR && !longestPair.chromosome().equals(startBreakend.chromosome()))
        {
            LOGGER.debug("cluster({}) spliting reciprocal DUP pair into 2 chains between chrs({} & {})",
                    cluster.id(), longestPair.chromosome(), startBreakend.chromosome());

            // remove the long chain link on the assumption that these breakends aren't joined
            cluster.getChains().clear();

            SvChain newChain = new SvChain(0);

            longestTILength = 0;

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
                    if(newChain != null && newChain.getLinkCount() > 0)
                    {
                        cluster.addChain(newChain, false);
                        newChain = null;
                    }
                }
            }

            if(newChain != null && newChain.getLinkCount() > 0)
                cluster.addChain(newChain, false);

            cluster.setResolved(true, resolvedType);
            cluster.setSyntheticData(0, longestTILength);
        }
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

        String resolvedType = RESOLVED_TYPE_RECIPROCAL_INV;

        LOGGER.debug("cluster({}) chain(links=({} longestTI={}) synLen({}) marked as {}",
                cluster.id(), totalLinks, tiPair.length(), syntheticLength, resolvedType);

        cluster.setResolved(true, resolvedType);
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

        String resolvedType = RESOLVED_TYPE_RECIPROCAL_TRANS ;

        LOGGER.debug("cluster({}) chain(links=({} longestTI={}) marked as {}",
                cluster.id(), totalLinks, longestTILength, resolvedType);

        cluster.setResolved(true, resolvedType);
        cluster.setSyntheticData(0, longestTILength);
    }

}
