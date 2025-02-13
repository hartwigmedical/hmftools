package com.hartwig.hmftools.pave.transval;

import com.hartwig.hmftools.common.codon.Nucleotides;
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

    @Override
    ChangeResult applyChange(ChangeContext context)
    {
        String residualBases = context.exonBasesAfterDeletion();
        AminoAcidSequence resultSequence = AminoAcidSequence.fromNucleotides(residualBases);
        String deleted = context.refBases();
        String altBases = deleted.substring(0,  1);

        return new ChangeResult(resultSequence, residualBases,context.positionOfChangeStartInStrand() - 1, deleted, altBases);
    }

//    @Override
//    TransvalHotspot convertToHotspot(final ChangeContext changeContext)
//    {
//        String deleted = changeContext.refBases();
//        String altBases = deleted.substring(0,  1);
//        return new TransvalHotspot(deleted, altBases, mGene.Chromosome, changeContext.positionOfChangeStartInStrand() - 1);
//    }
}
