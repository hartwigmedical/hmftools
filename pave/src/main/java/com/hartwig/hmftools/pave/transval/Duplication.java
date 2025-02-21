package com.hartwig.hmftools.pave.transval;

import java.util.Set;

import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptAminoAcids;
import com.hartwig.hmftools.common.gene.TranscriptData;

import org.jetbrains.annotations.NotNull;

class Duplication extends ProteinVariant
{
    public Duplication(@NotNull final GeneData gene,
            @NotNull final TranscriptData transcript,
            @NotNull final TranscriptAminoAcids aminoAcidSequence,
            @NotNull final AminoAcidRange refRange)
    {
        super(gene, transcript, aminoAcidSequence, refRange.startPosition(), refRange.length());
    }

    @NotNull
    @Override
    Set<ChangeResult> applyChange(ChangeContext changeContext)
    {
        int changeStart = changeContext.StartPositionInExon;
        int changeEnd = changeContext.FinishPositionInExon + 1;
        PaddedExon exon = changeContext.mExon;
        String bases = exon.baseSequenceWithDuplicationApplied(changeStart, changeEnd, changeContext.IsPositiveStrand);
        AminoAcidSequence acids = AminoAcidSequence.fromNucleotides(bases);
        String duplicated = changeContext.refBases();
        String refBase = duplicated.substring(0, 1);
        return Set.of( new ChangeResult(acids, bases, changeContext.positionOfChangeStartInStrand() - 1, refBase, duplicated));
    }

    @Override
    AminoAcidSequence variantSequence()
    {
        int startOfDuplicatedSection = positionOfFirstAlteredCodon() - 1;
        int endOfDuplicatedSection = startOfDuplicatedSection + this.mRefLength;
        return completeReferenceAminoAcidSequence().duplicateRange(startOfDuplicatedSection, endOfDuplicatedSection);
    }
}
