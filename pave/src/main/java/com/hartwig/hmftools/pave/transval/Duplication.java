package com.hartwig.hmftools.pave.transval;

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

    @Override
    ChangeResult applyChange(ChangeContext changeContext)
    {
        int changeStart = changeContext.StartPositionInExon;
        int changeEnd = changeContext.FinishPositionInExon + 1;
        PaddedExon exon = changeContext.ContainingExon;
        String bases = exon.baseSequenceWithDuplicationApplied(changeStart, changeEnd, changeContext.IsPositiveStrand);
        AminoAcidSequence resultSequence = AminoAcidSequence.fromNucleotides(bases);
        String duplicated = changeContext.refBases();
        String refBase = duplicated.substring(0,  1);
        return new ChangeResult(resultSequence, bases, changeContext.positionOfChangeStartInStrand() - 1, refBase, duplicated);
    }

    @Override
    AminoAcidSequence variantSequence()
    {
        String rawAAs = this.mAminoAcidSequence.AminoAcids;
        int startOfDuplicatedSection = positionOfFirstAlteredCodon() - 1;
        String left = rawAAs.substring(0, startOfDuplicatedSection);
        int endOfDuplicatedSection = startOfDuplicatedSection + this.mRefLength;
        String toDuplicate = rawAAs.substring(startOfDuplicatedSection, endOfDuplicatedSection);
        String right = rawAAs.substring(endOfDuplicatedSection);

        String duplicatedAAs = left + toDuplicate + toDuplicate + right;
        return AminoAcidSequence.parse(duplicatedAAs);
    }

//    @Override
//    TransvalHotspot convertToHotspot(final ChangeContext changeContext)
//    {
//        String duplicated = changeContext.refBases();
//        String refBase = duplicated.substring(0,  1);
//        return new TransvalHotspot(refBase, duplicated, mGene.Chromosome, changeContext.positionOfChangeStartInStrand() - 1);
//    }
}
