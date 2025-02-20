package com.hartwig.hmftools.pave.transval;

import java.util.HashSet;
import java.util.Set;

import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptAminoAcids;
import com.hartwig.hmftools.common.gene.TranscriptData;

import org.apache.commons.lang3.tuple.Pair;
import org.jetbrains.annotations.NotNull;

class SingleAminoAcidVariant extends ProteinVariant
{
    @NotNull
    private final AminoAcidSpecification Alt;

    public SingleAminoAcidVariant(
            @NotNull final GeneData gene,
            @NotNull final TranscriptData transcript,
            @NotNull final TranscriptAminoAcids aminoAcidSequence,
            @NotNull final AminoAcidSpecification alt)
    {
        super(gene, transcript, aminoAcidSequence, alt.mPosition, 1);
        this.Alt = alt;
    }

    @Override
    AminoAcidSequence variantSequence()
    {
        return referenceAminoAcidSequence().replace(positionOfFirstAlteredCodon(), Alt.value());
    }

    @Override
    int numberOfLeftShiftsToTry(ChangeContext changeContext)
    {
        return 0;
    }

    @NotNull
    @Override
    Set<ChangeResult> applyChange(ChangeContext context)
    {
        CodonWithinExons codon = context.codonForProteinChange(positionOfFirstAlteredCodon());
        Set<CodonChange> possibleVariants = codon.possibleVariantsGiving(Alt.value());
        Set<ChangeResult> result = new HashSet<>();
        possibleVariants.forEach(codonChange -> {
            CodonChange forwardStrandChange = context.IsPositiveStrand ? codonChange : codonChange.reverseComplement();
            CodonWithinExons forwardStrandCodon = context.IsPositiveStrand ? codon : codon.reverseComplement();
            Pair<String, String> refAlt = forwardStrandChange.differenceStrings();
            String forwardStrandReferenceBases = refAlt.getLeft();
            String forwardStrandReplacementBases = refAlt.getRight();
            int locationOnStrandOfChange = forwardStrandCodon.forwardStrandLocationOfChange(forwardStrandChange.positionOfFirstDifference());
            String exonBasesAfterChange = context.exonBasesWithReplacementAppliedAtStrandLocation(locationOnStrandOfChange, forwardStrandReplacementBases);
            AminoAcidSequence acids = AminoAcidSequence.fromNucleotides(exonBasesAfterChange);
            result.add(new ChangeResult(acids, exonBasesAfterChange, locationOnStrandOfChange, forwardStrandReferenceBases, forwardStrandReplacementBases));
        });
        return result;
    }

    public String altValue()
    {
        return Alt.symbol();
    }
}
