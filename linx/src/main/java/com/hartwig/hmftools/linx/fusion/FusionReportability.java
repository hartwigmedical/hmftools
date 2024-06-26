package com.hartwig.hmftools.linx.fusion;

import static java.lang.Math.min;
import static java.lang.Math.round;

import static com.hartwig.hmftools.common.fusion.KnownFusionType.PROMISCUOUS_ENHANCER_TARGET;
import static com.hartwig.hmftools.common.gene.TranscriptProteinData.BIOTYPE_NONSENSE_MED_DECAY;
import static com.hartwig.hmftools.common.gene.TranscriptProteinData.BIOTYPE_PROTEIN_CODING;
import static com.hartwig.hmftools.common.fusion.KnownFusionType.EXON_DEL_DUP;
import static com.hartwig.hmftools.common.fusion.KnownFusionType.IG_KNOWN_PAIR;
import static com.hartwig.hmftools.common.fusion.KnownFusionType.IG_PROMISCUOUS;
import static com.hartwig.hmftools.common.fusion.KnownFusionType.KNOWN_PAIR;
import static com.hartwig.hmftools.common.fusion.KnownFusionType.NONE;
import static com.hartwig.hmftools.common.fusion.KnownFusionType.PROMISCUOUS_3;
import static com.hartwig.hmftools.common.fusion.KnownFusionType.PROMISCUOUS_5;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.SGL;
import static com.hartwig.hmftools.linx.fusion.FusionConstants.FUSION_MAX_CHAIN_LINKS;
import static com.hartwig.hmftools.linx.fusion.FusionConstants.MAX_UPSTREAM_DISTANCE_IG_KNOWN;
import static com.hartwig.hmftools.linx.fusion.FusionConstants.MAX_UPSTREAM_DISTANCE_KNOWN;
import static com.hartwig.hmftools.linx.fusion.FusionConstants.MAX_UPSTREAM_DISTANCE_OTHER;
import static com.hartwig.hmftools.linx.fusion.FusionConstants.PROTEINS_REQUIRED_KEPT;
import static com.hartwig.hmftools.linx.fusion.FusionConstants.PROTEINS_REQUIRED_LOST;
import static com.hartwig.hmftools.linx.fusion.FusionConstants.SHORT_UNPHASED_DISTANCE_KNOWN;
import static com.hartwig.hmftools.common.linx.FusionReportableReason.CHAIN_LINKS;
import static com.hartwig.hmftools.common.linx.FusionReportableReason.CHAIN_TERMINATED;
import static com.hartwig.hmftools.common.linx.FusionReportableReason.EXON_SKIPPING;
import static com.hartwig.hmftools.common.linx.FusionReportableReason.INVALID_TRAVERSAL;
import static com.hartwig.hmftools.common.linx.FusionReportableReason.NOT_KNOWN;
import static com.hartwig.hmftools.common.linx.FusionReportableReason.NEG_SPLICE_ACC_DISTANCE;
import static com.hartwig.hmftools.common.linx.FusionReportableReason.NMD;
import static com.hartwig.hmftools.common.linx.FusionReportableReason.NON_DISRUPTIVE_CHAIN;
import static com.hartwig.hmftools.common.linx.FusionReportableReason.SGL_NOT_KNOWN;
import static com.hartwig.hmftools.common.linx.FusionReportableReason.UNPHASED_5P_UTR;
import static com.hartwig.hmftools.common.linx.FusionReportableReason.UNPHASED_NOT_KNOWN;
import static com.hartwig.hmftools.common.linx.FusionReportableReason.UNPHASED_SHORT;
import static com.hartwig.hmftools.common.linx.FusionReportableReason.PRE_GENE_DISTANCE;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.fusion.KnownFusionType;
import com.hartwig.hmftools.common.linx.FusionReportableReason;
import com.hartwig.hmftools.linx.gene.BreakendTransData;

public class FusionReportability
{
    public static boolean isReportable(final GeneFusion fusion)
    {
        return determineReportability(fusion).isEmpty();
    }

    public static List<FusionReportableReason> determineReportability(final GeneFusion fusion)
    {
        // first check whether a fusion is known or not - a key requirement of it being potentially reportable
        if(fusion.knownType() == NONE || fusion.knownType() == IG_PROMISCUOUS)
            return Lists.newArrayList(NOT_KNOWN);

        List<FusionReportableReason> nonReportableReasons = Lists.newArrayList();

        final BreakendTransData upTrans = fusion.upstreamTrans();
        final BreakendTransData downTrans = fusion.downstreamTrans();

        if(!fusion.phaseMatched())
        {
            if(fusion.knownType() == EXON_DEL_DUP && fusion.isExonic())
            {
                // currently allowed for specific exon-to-exon fusions
            }
            else
            {
                if(fusion.knownType() != KNOWN_PAIR)
                    nonReportableReasons.add(UNPHASED_NOT_KNOWN);

                if(!upTrans.nonCoding() && !upTrans.preCoding())
                    nonReportableReasons.add(UNPHASED_5P_UTR);

                if(fusion.sameChromosome() && fusion.distance() <= SHORT_UNPHASED_DISTANCE_KNOWN)
                    nonReportableReasons.add(UNPHASED_SHORT);
            }
        }

        if(upTrans.gene().type() == SGL || downTrans.gene().type() == SGL)
        {
            if(fusion.knownType() != KNOWN_PAIR && fusion.knownType() != IG_KNOWN_PAIR)
                nonReportableReasons.add(SGL_NOT_KNOWN);
        }

        // set limits on how far upstream the breakend can be - adjusted for whether the fusions is known or not
        int maxUpstreamDistance = getMaxUpstreamDistance(fusion);

        if(upTrans.getDistanceUpstream() > maxUpstreamDistance || downTrans.getDistanceUpstream() > maxUpstreamDistance)
            nonReportableReasons.add(PRE_GENE_DISTANCE);

        if(downTrans.bioType().equals(BIOTYPE_NONSENSE_MED_DECAY))
            nonReportableReasons.add(NMD);

        if(!fusion.isIG() && downTrans.hasNegativePrevSpliceAcceptorDistance())
            nonReportableReasons.add(NEG_SPLICE_ACC_DISTANCE);

        if(!permittedExonSkipping(fusion))
            nonReportableReasons.add(EXON_SKIPPING);

        if(fusion.isTerminated() && !allowSuspectChains(fusion.knownType()))
            nonReportableReasons.add(CHAIN_TERMINATED);

        if(fusion.nonDisruptiveChain())
            nonReportableReasons.add(NON_DISRUPTIVE_CHAIN);

        if(!fusion.validChainTraversal() && !allowSuspectChains(fusion.knownType()))
            nonReportableReasons.add(INVALID_TRAVERSAL);

        if(fusion.getChainLinks() > FUSION_MAX_CHAIN_LINKS)
            nonReportableReasons.add(CHAIN_LINKS);

        return nonReportableReasons;
    }

    private static int getMaxUpstreamDistance(final GeneFusion fusion)
    {
        if(fusion.knownType() == IG_KNOWN_PAIR)
            return MAX_UPSTREAM_DISTANCE_IG_KNOWN;
        else if(fusion.knownType() == KNOWN_PAIR  || fusion.isHighImpactPromiscuous())
            return MAX_UPSTREAM_DISTANCE_KNOWN;
        else
            return MAX_UPSTREAM_DISTANCE_OTHER;
    }

    public static boolean allowSuspectChains(final KnownFusionType type)
    {
        return (type == KNOWN_PAIR || type == EXON_DEL_DUP || type == IG_KNOWN_PAIR);
    }

    public static GeneFusion findTopPriorityFusion(final List<GeneFusion> fusions)
    {
        GeneFusion topFusion = null;

        // form a score by allocating 0/1 or length value to each power of 10 descending
        double highestScore = 0;

        for(final GeneFusion fusion : fusions)
        {
            double fusionPriorityScore = calcFusionPriority(fusion);
            fusion.setPriority(fusionPriorityScore);

            if(fusionPriorityScore > highestScore)
            {
                topFusion = fusion;
                highestScore = fusionPriorityScore;
            }
        }

        return topFusion;
    }

    private static double calcFusionPriority(final GeneFusion fusion)
    {
        final BreakendTransData upTrans = fusion.upstreamTrans();
        final BreakendTransData downTrans = fusion.downstreamTrans();

            /* prioritisation rules:
                - Reportable
                - Known pair or DEL-DUP
                - inframe
                - chain not terminated for known fusions
                - 3’ partner biotype is protein_coding
                - No exons skipped
                - Best 3' partner by canonical, not NMD then coding bases (or exon count if not coding)
                - Best 5' partner by canonical, protein-coding then coding bases
            */

        double fusionPriorityScore = 0;
        double factor = 1000000;

        // Reportable
        if(fusion.reportable())
            fusionPriorityScore += factor;

        factor /= 10;

        // 0. Known pair
        if(isHighPriorityType(fusion.knownType()))
            fusionPriorityScore += factor;

        factor /= 10;

        // 1. Phase matched
        if(fusion.phaseMatched())
            fusionPriorityScore += factor;

        factor /= 10;

        // 2. Chain not terminated (only applicable for chained & known fusions)
        if(!fusion.isTerminated() && fusion.validChainTraversal())
            fusionPriorityScore += factor;

        factor /= 10;

        // 3' protein coding
        if(downTrans.bioType().equals(BIOTYPE_PROTEIN_CODING))
            fusionPriorityScore += factor;

        factor /= 10;

        // 4. Not skipping exons
        if(fusion.getExonsSkipped(true) == 0 && fusion.getExonsSkipped(false) == 0)
            fusionPriorityScore += factor;

        factor /= 10;

        // 5. Best 3' partner
        if(downTrans.isCanonical())
            fusionPriorityScore += factor;

        factor /= 10;

        if(!downTrans.bioType().equals(BIOTYPE_NONSENSE_MED_DECAY))
            fusionPriorityScore += factor;

        factor /= 100;

        int length = downTrans.isCoding() ? downTrans.CodingBases : downTrans.exonCount();

        // will be a range between 1-99 * current factor
        length = min(round(length/10), 99);
        fusionPriorityScore += length * factor;

        factor /= 10;

        // 6. Best 5' partner
        if(upTrans.isCanonical())
            fusionPriorityScore += factor;

        factor /= 10;

        if(upTrans.bioType().equals(BIOTYPE_PROTEIN_CODING))
            fusionPriorityScore += factor;

        factor /= 100;

        length = upTrans.isCoding() ? upTrans.CodingBases : upTrans.exonCount();
        length = min(round(length/10), 99);
        fusionPriorityScore += length * factor;

        factor /= 10;

        fusionPriorityScore += upTrans.gene().jcn() * factor;

        return fusionPriorityScore;
    }

    public static boolean isHighPriorityType(final KnownFusionType type)
    {
        return (type == KNOWN_PAIR || type == IG_KNOWN_PAIR || type == EXON_DEL_DUP || type == PROMISCUOUS_ENHANCER_TARGET);
    }

    public static boolean checkProteinDomains(final KnownFusionType type)
    {
        return (type == PROMISCUOUS_5 || type == PROMISCUOUS_3 || type == KNOWN_PAIR);
    }

    public static boolean validProteinDomains(final GeneFusion fusion)
    {
        final BreakendTransData downTrans = fusion.downstreamTrans();

        long requiredKeptButLost = PROTEINS_REQUIRED_KEPT.stream().filter(f -> downTrans.getProteinFeaturesLost().contains(f)).count();
        long requiredLostButKept = PROTEINS_REQUIRED_LOST.stream().filter(f -> downTrans.getProteinFeaturesKept().contains(f)).count();

        return requiredKeptButLost == 0 && requiredLostButKept == 0;
    }

    private static boolean permittedExonSkipping(final GeneFusion fusion)
    {
        if(fusion.knownType() == KNOWN_PAIR || fusion.knownType() == IG_KNOWN_PAIR || fusion.knownType() == IG_PROMISCUOUS || fusion.knownType() == EXON_DEL_DUP)
            return true;

        if(fusion.isHighImpactPromiscuous() || fusion.knownExons())
            return true;

        // otherwise not allowed
        return fusion.getExonsSkipped(true) == 0 && fusion.getExonsSkipped(false) == 0;
    }
}
