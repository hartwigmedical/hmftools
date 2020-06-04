package com.hartwig.hmftools.linx.fusion;

import static java.lang.Math.min;
import static java.lang.Math.round;

import static com.hartwig.hmftools.common.ensemblcache.TranscriptProteinData.BIOTYPE_NONSENSE_MED_DECAY;
import static com.hartwig.hmftools.common.ensemblcache.TranscriptProteinData.BIOTYPE_PROTEIN_CODING;
import static com.hartwig.hmftools.common.fusion.KnownFusionType.EXON_DEL_DUP;
import static com.hartwig.hmftools.common.fusion.KnownFusionType.IG_KNOWN_PAIR;
import static com.hartwig.hmftools.common.fusion.KnownFusionType.IG_PROMISCUOUS;
import static com.hartwig.hmftools.common.fusion.KnownFusionType.KNOWN_PAIR;
import static com.hartwig.hmftools.common.fusion.KnownFusionType.NONE;
import static com.hartwig.hmftools.common.fusion.KnownFusionType.PROMISCUOUS_3;
import static com.hartwig.hmftools.common.fusion.KnownFusionType.PROMISCUOUS_5;
import static com.hartwig.hmftools.linx.fusion.FusionConstants.MAX_UPSTREAM_DISTANCE_KNOWN;
import static com.hartwig.hmftools.linx.fusion.FusionConstants.MAX_UPSTREAM_DISTANCE_OTHER;

import java.util.List;

import com.hartwig.hmftools.common.fusion.KnownFusionType;
import com.hartwig.hmftools.common.fusion.Transcript;

public class FusionReportability
{
    public static boolean couldBeReportable(GeneFusion fusion)
    {
        if(!fusion.phaseMatched() || fusion.neoEpitopeOnly())
            return false;

        // first check whether a fusion is known or not - a key requirement of it being potentially reportable
        if (fusion.knownType() == NONE)
            return false;

        // set limits on how far upstream the breakend can be - adjusted for whether the fusions is known or not
        int maxUpstreamDistance = fusion.knownType() == KNOWN_PAIR ?
                MAX_UPSTREAM_DISTANCE_KNOWN : MAX_UPSTREAM_DISTANCE_OTHER;

        final Transcript upTrans = fusion.upstreamTrans();
        final Transcript downTrans = fusion.downstreamTrans();

        if(upTrans.getDistanceUpstream() > maxUpstreamDistance || downTrans.getDistanceUpstream() > maxUpstreamDistance)
            return false;

        if(downTrans.bioType().equals(BIOTYPE_NONSENSE_MED_DECAY))
            return false;

        if(downTrans.hasNegativePrevSpliceAcceptorDistance())
            return false;

        if(!allowExonSkipping(fusion.knownType()))
        {
            if(fusion.getExonsSkipped(true) > 0 || fusion.getExonsSkipped(false) > 0)
                return false;
        }

        return true;
    }

    public static GeneFusion determineReportableFusion(final List<GeneFusion> fusions, boolean requireReportable)
    {
        GeneFusion reportableFusion = null;

        // form a score by allocating 0/1 or length value to each power of 10 descending
        double highestScore = 0;

        for(final GeneFusion fusion : fusions)
        {
            if(requireReportable && !couldBeReportable(fusion))
                continue;

            double fusionPriorityScore = calcFusionPriority(fusion);
            fusion.setPriority(fusionPriorityScore);

            if(fusionPriorityScore > highestScore)
            {
                reportableFusion = fusion;
                highestScore = fusionPriorityScore;
            }
        }

        return reportableFusion;
    }

    private static double calcFusionPriority(final GeneFusion fusion)
    {
        // first check whether a fusion is known or not - a key requirement of it being potentially reportable
        final Transcript upTrans = fusion.upstreamTrans();
        final Transcript downTrans = fusion.downstreamTrans();

            /* prioritisation rules:
            1. inframe
            2. chain not terminated for known fusions
            3. 3’ partner biotype is protein_coding
            4. No exons skipped
            5. Best 3' partner by canonical, not NMD then coding bases (or exon count if not coding)
            6. Best 5' partner by canonical, protein-coding then coding bases
            */

        double fusionPriorityScore = 0;
        double factor = 1000000;

        // 1. Phase matched
        if(fusion.phaseMatched())
            fusionPriorityScore += factor;

        factor /= 10;

        // 2. Chain not terminated (only applicable for chained & known fusions)
        if(!fusion.isTerminated())
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

        int length = downTrans.isCoding() ? downTrans.calcCodingBases() : downTrans.ExonMax;

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

        length = upTrans.isCoding() ? upTrans.calcCodingBases() : upTrans.ExonMax;
        length = min(round(length/10), 99);
        fusionPriorityScore += length * factor;

        return fusionPriorityScore;
    }

    public static boolean checkProteinDomains(final KnownFusionType type)
    {
        return (type == PROMISCUOUS_5 || type == PROMISCUOUS_3);
    }

    public static boolean allowExonSkipping(final KnownFusionType type)
    {
        return (type == KNOWN_PAIR || type == EXON_DEL_DUP || type == IG_KNOWN_PAIR || type == IG_PROMISCUOUS);
    }

}
