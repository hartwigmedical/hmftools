package com.hartwig.hmftools.pave.reverse;

import java.util.Set;

import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptAminoAcids;
import com.hartwig.hmftools.common.gene.TranscriptData;

import org.jetbrains.annotations.NotNull;

class Deletion extends ProteinVariant
{
    public Deletion(@NotNull final GeneData gene,
            @NotNull final TranscriptData transcript,
            @NotNull final TranscriptAminoAcids aminoAcidSequence,
            @NotNull final AminoAcidRange refRange)
    {
        super(gene, transcript, aminoAcidSequence, refRange.startPosition(), refRange.length());
    }

    @NotNull
    @Override
    Set<ChangeResult> applyChange(ChangeContext context)
    {
        int start = context.StartPositionInExon;
        int end = context.FinishPositionInExon + 1;
        String bases = context.mExon.baseSequenceWithFramePreservingDeletionApplied(start, end, context.IsPositiveStrand);
        AminoAcidSequence resultSequence = AminoAcidSequence.fromNucleotides(bases);
        String deleted = context.refBases();
        String altBases = deleted.substring(0,  1);
        return Set.of(new ChangeResult(resultSequence, bases,context.positionOfChangeStartInStrand(), deleted, altBases));
    }

    @Override
    AminoAcidSequence variantSequence()
    {
        int startOfDeletedSection = positionOfFirstAlteredCodon() - 1;
        int endOfDeletedSection = startOfDeletedSection + this.mRefLength;
        return completeReferenceAminoAcidSequence().deleteRange(startOfDeletedSection, endOfDeletedSection);
    }
}
