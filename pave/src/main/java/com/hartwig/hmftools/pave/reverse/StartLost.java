package com.hartwig.hmftools.pave.reverse;

import static com.hartwig.hmftools.pave.reverse.AminoAcid.START;

import java.util.HashSet;
import java.util.Set;

import com.google.common.base.Preconditions;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptAminoAcids;
import com.hartwig.hmftools.common.gene.TranscriptData;

import org.jetbrains.annotations.NotNull;

class StartLost extends SingleCodonVariant
{
    StartLost(@NotNull final GeneData gene,
            @NotNull final TranscriptData transcript,
            @NotNull final TranscriptAminoAcids aminoAcidSequence,
            @NotNull final AminoAcidRange refRange)
    {
        super(gene, transcript, aminoAcidSequence, refRange.startPosition());
        Preconditions.checkArgument(refRange.length() == 1);
        Preconditions.checkArgument(refRange.startPosition() == 1);
        Preconditions.checkArgument(refRange.aminoAcidAtStart().equals(START));
    }

    @Override
    AminoAcidSequence variantSequence()
    {
        return AminoAcidSequence.empty();
    }

    @NotNull
    @Override
    Set<CodonChange> possibleVariants(@NotNull final CodonWithinExons codon)
    {
        Set<CodonChange> changes = new HashSet<>();
        changes.add(new CodonChange("ATG", "ATA"));
        changes.add(new CodonChange("ATG", "ATC"));
        changes.add(new CodonChange("ATG", "ATT"));
        changes.add(new CodonChange("ATG", "AAG"));
        changes.add(new CodonChange("ATG", "ACG"));
        changes.add(new CodonChange("ATG", "AGG"));
        changes.add(new CodonChange("ATG", "CTG"));
        changes.add(new CodonChange("ATG", "GTG"));
        changes.add(new CodonChange("ATG", "TTG"));
        return changes;
    }

    @Override
    boolean isConsistentWithThisVariant(AminoAcidSequence candidate)
    {
        return candidate.length() == 0;
    }

    @Override
    AminoAcidSequence convertBaseSequence(String bases)
    {
        Preconditions.checkArgument(bases.length() > 3);
        Preconditions.checkArgument(!bases.startsWith("ATG"));
        return AminoAcidSequence.empty();
    }
}
