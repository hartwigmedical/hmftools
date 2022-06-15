package com.hartwig.hmftools.linx.fusion;

import static java.lang.Math.min;
import static java.lang.Math.round;

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
import static com.hartwig.hmftools.linx.fusion.FusionConstants.MAX_UPSTREAM_DISTANCE_KNOWN;
import static com.hartwig.hmftools.linx.fusion.FusionConstants.MAX_UPSTREAM_DISTANCE_OTHER;
import static com.hartwig.hmftools.linx.fusion.FusionConstants.SHORT_UNPHASED_DISTANCE_KNOWN;
import static com.hartwig.hmftools.linx.fusion.ReportableReason.CHAIN_LINKS;
import static com.hartwig.hmftools.linx.fusion.ReportableReason.CHAIN_TERMINATED;
import static com.hartwig.hmftools.linx.fusion.ReportableReason.EXON_SKIPPING;
import static com.hartwig.hmftools.linx.fusion.ReportableReason.INVALID_TRAVERSAL;
import static com.hartwig.hmftools.linx.fusion.ReportableReason.KNOWN_TYPE;
import static com.hartwig.hmftools.linx.fusion.ReportableReason.NEG_SPLICE_ACC_DISTANCE;
import static com.hartwig.hmftools.linx.fusion.ReportableReason.NMD;
import static com.hartwig.hmftools.linx.fusion.ReportableReason.NON_DISRUPTIVE_CHAIN;
import static com.hartwig.hmftools.linx.fusion.ReportableReason.OK;
import static com.hartwig.hmftools.linx.fusion.ReportableReason.SGL_NOT_KNOWN;
import static com.hartwig.hmftools.linx.fusion.ReportableReason.UNPHASED_5P_UTR;
import static com.hartwig.hmftools.linx.fusion.ReportableReason.UNPHASED_NOT_KNOWN;
import static com.hartwig.hmftools.linx.fusion.ReportableReason.UNPHASED_SHORT;
import static com.hartwig.hmftools.linx.fusion.ReportableReason.PRE_GENE_DISTANCE;

import java.util.List;

import com.hartwig.hmftools.common.fusion.KnownFusionType;
import com.hartwig.hmftools.linx.gene.BreakendTransData;

import org.apache.commons.compress.utils.Lists;

public class FusionReportability
{
    private static final List<String> mProteinsRequiredKept = Lists.newArrayList();
    private static final List<String> mProteinsRequiredLost = Lists.newArrayList();

    public static ReportableReason determineReportability(final GeneFusion fusion)
    {
        // first check whether a fusion is known or not - a key requirement of it being potentially reportable
        if (fusion.knownType() == NONE || fusion.knownType() == IG_PROMISCUOUS)
            return KNOWN_TYPE;

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
                    return UNPHASED_NOT_KNOWN;

                if(!upTrans.nonCoding() && !upTrans.preCoding())
                    return UNPHASED_5P_UTR;

                if(fusion.sameChromosome() && fusion.distance() <= SHORT_UNPHASED_DISTANCE_KNOWN)
                    return UNPHASED_SHORT;
            }
        }

        if(upTrans.gene().type() == SGL || downTrans.gene().type() == SGL)
        {
            if(fusion.knownType() != KNOWN_PAIR && fusion.knownType() != IG_KNOWN_PAIR)
                return SGL_NOT_KNOWN;
        }

        // set limits on how far upstream the breakend can be - adjusted for whether the fusions is known or not
        int maxUpstreamDistance = getMaxUpstreamDistance(fusion);

        if(upTrans.getDistanceUpstream() > maxUpstreamDistance || downTrans.getDistanceUpstream() > maxUpstreamDistance)
            return PRE_GENE_DISTANCE;

        if(downTrans.bioType().equals(BIOTYPE_NONSENSE_MED_DECAY))
            return NMD;

        if(!fusion.isIG() && downTrans.hasNegativePrevSpliceAcceptorDistance())
            return NEG_SPLICE_ACC_DISTANCE;

        if(!permittedExonSkipping(fusion))
            return EXON_SKIPPING;

        if(fusion.isTerminated() && !allowSuspectChains(fusion.knownType()))
            return CHAIN_TERMINATED;

        if(fusion.nonDisruptiveChain())
            return NON_DISRUPTIVE_CHAIN;

        if(!fusion.validChainTraversal() && !allowSuspectChains(fusion.knownType()))
            return INVALID_TRAVERSAL;

        if(fusion.getChainLinks() > FUSION_MAX_CHAIN_LINKS)
            return CHAIN_LINKS;

        return OK;
    }

    private static int getMaxUpstreamDistance(final GeneFusion fusion)
    {
        if(fusion.knownType() == KNOWN_PAIR || fusion.knownType() == IG_KNOWN_PAIR || fusion.isHighImpactPromiscuous())
            return MAX_UPSTREAM_DISTANCE_KNOWN;

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
        if(fusion.knownType() == KNOWN_PAIR || fusion.knownType() == EXON_DEL_DUP)
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

        return fusionPriorityScore;
    }

    public static boolean checkProteinDomains(final KnownFusionType type)
    {
        return (type == PROMISCUOUS_5 || type == PROMISCUOUS_3 || type == KNOWN_PAIR);
    }

    public static boolean validProteinDomains(final GeneFusion fusion)
    {
        final BreakendTransData downTrans = fusion.downstreamTrans();

        long requiredKeptButLost =
                mProteinsRequiredKept.stream().filter(f -> downTrans.getProteinFeaturesLost().contains(f)).count();
        long requiredLostButKept =
                mProteinsRequiredLost.stream().filter(f -> downTrans.getProteinFeaturesKept().contains(f)).count();

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

    public static void populateRequiredProteins()
    {
        mProteinsRequiredLost.add("Raf-like Ras-binding");

        mProteinsRequiredKept.add("Ets domain");
        mProteinsRequiredKept.add("Protein kinase domain");
        mProteinsRequiredKept.add("Epidermal growth factor-like domain");
        mProteinsRequiredKept.add("Ankyrin repeat-containing domain");
        mProteinsRequiredKept.add("Basic-leucine zipper domain");
        mProteinsRequiredKept.add("High mobility group box domain");
        mProteinsRequiredKept.add("Bromodomain");
    }

}
