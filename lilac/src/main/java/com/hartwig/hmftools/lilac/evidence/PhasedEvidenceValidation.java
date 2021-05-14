package com.hartwig.hmftools.lilac.evidence;

import static com.hartwig.hmftools.lilac.LilacConfig.LL_LOGGER;

import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.lilac.seq.HlaSequenceLoci;

public class PhasedEvidenceValidation
{
    public static void logInconsistentEvidence(final String gene, final List<PhasedEvidence> evidence, final List<HlaSequenceLoci> candidates)
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
}
