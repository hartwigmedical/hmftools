package com.hartwig.hmftools.lilac.evidence;

import static com.hartwig.hmftools.lilac.LilacConfig.LL_LOGGER;

import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.lilac.seq.HlaSequenceLoci;

public class PhasedEvidenceValidation
{
    // TODO - rename
    public static void validateExpected(final String gene, final List<PhasedEvidence> evidence, final List<HlaSequenceLoci> candidates)
    {
        List<HlaSequenceLoci> expectedSequences = candidates.stream().filter(x -> x.getAllele().Gene.equals(gene)).collect(Collectors.toList());

        for (HlaSequenceLoci sequence : expectedSequences)
        {
            for (PhasedEvidence phasedEvidence : evidence)
            {
                if (!sequence.consistentWith(phasedEvidence))
                {
                    LL_LOGGER.warn("Expected allele {} filtered by {}", sequence.getAllele(), phasedEvidence);
                }
            }
        }
    }

    // TODO - rename
    public static void validateAgainstFinalCandidates(
            final String gene, final List<PhasedEvidence> evidence, final List<HlaSequenceLoci> candidates)
    {
        for(PhasedEvidence inconsistentEvidence2 : unmatchedEvidence(evidence, candidates))
        {
            LL_LOGGER.warn("HLA-{} phased evidence not found in candidates: {}", gene, inconsistentEvidence2);
        }
    }

    private static List<PhasedEvidence> unmatchedEvidence(final List<PhasedEvidence> evidence, final List<HlaSequenceLoci> candidates)
    {
        return evidence.stream().map(x -> x.inconsistentEvidence(candidates))
                .filter(x -> x.getEvidence().isEmpty()).collect(Collectors.toList());
    }
}
