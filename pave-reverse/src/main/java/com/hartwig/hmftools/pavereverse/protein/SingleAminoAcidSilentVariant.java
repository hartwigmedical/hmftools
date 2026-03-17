package com.hartwig.hmftools.pavereverse.protein;

import java.util.Set;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptAminoAcids;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.pavereverse.aminoacids.AminoAcidSequence;
import com.hartwig.hmftools.pavereverse.base.CodonWithinExons;

public class SingleAminoAcidSilentVariant extends SingleCodonVariant
{
    public SingleAminoAcidSilentVariant(
            GeneData gene,
            TranscriptData transcript,
            TranscriptAminoAcids aminoAcidSequence,
            int position)
    {
        super(gene, transcript, aminoAcidSequence, position);
    }

    @Override
    AminoAcidSequence variantSequence()
    {
        return completeReferenceAminoAcidSequence();
    }

    @Override
    Set<CodonChange> possibleVariants(CodonWithinExons codon)
    {
        return codon.possibleVariantsGiving(variantSequence().get(positionOfFirstAlteredCodon() - 1))
                .stream()
                .filter(c -> c.editDistance() > 0)
                .collect(Collectors.toSet());
    }

    public String altValue()
    {
        return variantSequence().symbolAt(positionOfFirstAlteredCodon() - 1);
    }
}
