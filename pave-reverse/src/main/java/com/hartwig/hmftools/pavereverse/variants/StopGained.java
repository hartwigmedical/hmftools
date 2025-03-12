package com.hartwig.hmftools.pavereverse.variants;

import java.util.Set;

import com.google.common.base.Preconditions;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptAminoAcids;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.pavereverse.aa.AminoAcid;
import com.hartwig.hmftools.pavereverse.aa.AminoAcidRange;
import com.hartwig.hmftools.pavereverse.aa.AminoAcidSequence;
import com.hartwig.hmftools.pavereverse.base.CodonWithinExons;

public class StopGained extends SingleCodonVariant
{
    private static final AminoAcid STOP = new AminoAcid("X");
    public final AminoAcid FirstChangedAminoAcid;

    public StopGained(GeneData gene,
            TranscriptData transcript,
            TranscriptAminoAcids aminoAcidSequence,
            AminoAcidRange refRange)
    {
        super(gene, transcript, aminoAcidSequence, refRange.startPosition());
        Preconditions.checkArgument(refRange.length() == 1);
        FirstChangedAminoAcid = refRange.aminoAcidAtStart();
    }

    @Override
    AminoAcidSequence variantSequence()
    {
        return completeReferenceAminoAcidSequence().subsequenceUpToInclusive(positionOfFirstAlteredCodon() - 1);
    }

    @Override
    Set<CodonChange> possibleVariants(CodonWithinExons codon)
    {
        return codon.possibleVariantsGivingStop();
    }

    @Override
    boolean isConsistentWithThisVariant(AminoAcidSequence candidate)
    {
        // The sequences need to match up to the change point and then stop.
        if(doesNotMatchVariantUpToLastAminoAcid(candidate))
        {
            return false;
        }
        return candidate.get(variantSequence().length()).equals(STOP);
    }
}
