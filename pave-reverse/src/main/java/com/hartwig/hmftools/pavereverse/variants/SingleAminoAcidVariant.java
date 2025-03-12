package com.hartwig.hmftools.pavereverse.variants;

import java.util.Set;

import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptAminoAcids;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.pavereverse.aa.AminoAcidSequence;
import com.hartwig.hmftools.pavereverse.aa.AminoAcidSpecification;
import com.hartwig.hmftools.pavereverse.base.CodonWithinExons;

public class SingleAminoAcidVariant extends SingleCodonVariant
{
    private final AminoAcidSpecification mAlt;

    public SingleAminoAcidVariant(
            GeneData gene,
            TranscriptData transcript,
            TranscriptAminoAcids aminoAcidSequence,
            AminoAcidSpecification alt)
    {
        super(gene, transcript, aminoAcidSequence, alt.Position);
        this.mAlt = alt;
    }

    @Override
    AminoAcidSequence variantSequence()
    {
        return completeReferenceAminoAcidSequence().replace(positionOfFirstAlteredCodon(), mAlt.value());
    }

    @Override
    Set<CodonChange> possibleVariants(CodonWithinExons codon)
    {
        return codon.possibleVariantsGiving(mAlt.value());
    }

    public String altValue()
    {
        return mAlt.symbol();
    }
}
