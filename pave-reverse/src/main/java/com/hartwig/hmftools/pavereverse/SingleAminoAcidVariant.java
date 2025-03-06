package com.hartwig.hmftools.pavereverse;

import java.util.Set;

import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptAminoAcids;
import com.hartwig.hmftools.common.gene.TranscriptData;

import org.jetbrains.annotations.NotNull;

class SingleAminoAcidVariant extends SingleCodonVariant
{
    @NotNull
    private final AminoAcidSpecification mAlt;

    public SingleAminoAcidVariant(
            @NotNull final GeneData gene,
            @NotNull final TranscriptData transcript,
            @NotNull final TranscriptAminoAcids aminoAcidSequence,
            @NotNull final AminoAcidSpecification alt)
    {
        super(gene, transcript, aminoAcidSequence, alt.mPosition);
        this.mAlt = alt;
    }

    @Override
    AminoAcidSequence variantSequence()
    {
        return completeReferenceAminoAcidSequence().replace(positionOfFirstAlteredCodon(), mAlt.value());
    }

    @NotNull
    @Override
    Set<CodonChange> possibleVariants(@NotNull final CodonWithinExons codon)
    {
        return codon.possibleVariantsGiving(mAlt.value());
    }

    public String altValue()
    {
        return mAlt.symbol();
    }
}
