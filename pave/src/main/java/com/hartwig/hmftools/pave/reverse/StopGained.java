package com.hartwig.hmftools.pave.reverse;

import java.util.Set;

import com.google.common.base.Preconditions;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptAminoAcids;
import com.hartwig.hmftools.common.gene.TranscriptData;

import org.jetbrains.annotations.NotNull;

class StopGained extends SingleCodonVariant
{
    private static final AminoAcid STOP = new AminoAcid("X");
    @NotNull
    final AminoAcid mFirstChangedAminoAcid;

    StopGained(@NotNull final GeneData gene,
            @NotNull final TranscriptData transcript,
            @NotNull final TranscriptAminoAcids aminoAcidSequence,
            @NotNull final AminoAcidRange refRange)
    {
        super(gene, transcript, aminoAcidSequence, refRange.startPosition());
        Preconditions.checkArgument(refRange.length() == 1);
        mFirstChangedAminoAcid = refRange.aminoAcidAtStart();
    }

    @Override
    AminoAcidSequence variantSequence()
    {
        return completeReferenceAminoAcidSequence().subsequenceUpToInclusive(positionOfFirstAlteredCodon() - 1);
    }

    @NotNull
    @Override
    Set<CodonChange> possibleVariants(@NotNull final CodonWithinExons codon)
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
