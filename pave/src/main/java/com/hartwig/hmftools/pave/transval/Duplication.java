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
        return changeContext.applyDuplication();
    }

    @Override
    TransvalHotspot convertToHotspot(final ChangeContext changeContext)
    {
        String duplicated = changeContext.refBases();
        String refBase = duplicated.substring(0,  1);
        return new TransvalHotspot(refBase, duplicated, Gene.Chromosome, changeContext.positionOfChangeStartInStrand() - 1);
    }
}
