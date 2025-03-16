package com.hartwig.hmftools.pavereverse.protein;

import java.util.Set;

import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptAminoAcids;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.pavereverse.aa.AminoAcidRange;
import com.hartwig.hmftools.pavereverse.aa.AminoAcidSequence;
import com.hartwig.hmftools.pavereverse.base.ChangeContext;

public class Deletion extends ProteinVariant
{
    public Deletion(GeneData gene,
            TranscriptData transcript,
            TranscriptAminoAcids aminoAcidSequence,
            AminoAcidRange refRange)
    {
        super(gene, transcript, aminoAcidSequence, refRange.startPosition(), refRange.length());
    }

    @Override
    Set<ChangeResult> applyChange(ChangeContext context)
    {
        int start = context.StartPositionInExon;
        int end = context.FinishPositionInExon + 1;
        String bases = context.Exon.baseSequenceWithFramePreservingDeletionApplied(start, end, context.IsPositiveStrand);
        AminoAcidSequence resultSequence = AminoAcidSequence.fromNucleotides(bases);
        String deleted = context.refBases();
        String altBases = deleted.substring(0, 1);
        return Set.of(new ChangeResult(resultSequence, bases, context.positionOfChangeStartInStrand(), deleted, altBases));
    }

    @Override
    AminoAcidSequence variantSequence()
    {
        int startOfDeletedSection = positionOfFirstAlteredCodon() - 1;
        int endOfDeletedSection = startOfDeletedSection + this.RefLength;
        return completeReferenceAminoAcidSequence().deleteRange(startOfDeletedSection, endOfDeletedSection);
    }
}
