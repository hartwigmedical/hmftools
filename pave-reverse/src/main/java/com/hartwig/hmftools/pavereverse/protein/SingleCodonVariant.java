package com.hartwig.hmftools.pavereverse.protein;

import java.util.HashSet;
import java.util.Set;

import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptAminoAcids;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.pavereverse.aa.AminoAcidSequence;
import com.hartwig.hmftools.pavereverse.base.ChangeContext;
import com.hartwig.hmftools.pavereverse.base.CodonWithinExons;

import org.apache.commons.lang3.tuple.Pair;

abstract class SingleCodonVariant extends ProteinVariant
{
    SingleCodonVariant(
            GeneData gene,
            TranscriptData transcript,
            TranscriptAminoAcids aminoAcidSequence,
            final int positionOfFirstAlteredCodon)
    {
        super(gene, transcript, aminoAcidSequence, positionOfFirstAlteredCodon, 1);
    }

    @Override
    int numberOfLeftShiftsToTry(ChangeContext changeContext)
    {
        return 0;
    }

    abstract Set<CodonChange> possibleVariants(CodonWithinExons codon);

    @Override
    public Set<ChangeResult> applyChange(ChangeContext context)
    {
        CodonWithinExons codon = context.codonForProteinChange(positionOfFirstAlteredCodon());
        Set<CodonChange> possibleVariants = possibleVariants(codon);
        Set<ChangeResult> result = new HashSet<>();
        possibleVariants.forEach(codonChange ->
        {
            CodonChange forwardStrandChange = context.IsPositiveStrand ? codonChange : codonChange.reverseComplement();
            CodonWithinExons forwardStrandCodon = context.IsPositiveStrand ? codon : codon.reverseComplement();
            Pair<String, String> refAlt = forwardStrandChange.differenceStrings();
            String forwardStrandReferenceBases = refAlt.getLeft();
            String forwardStrandReplacementBases = refAlt.getRight();
            int locationOnStrandOfChange =
                    forwardStrandCodon.forwardStrandLocationOfChange(forwardStrandChange.positionOfFirstDifference());
            String exonBasesAfterChange =
                    context.exonBasesWithReplacementAppliedAtStrandLocation(locationOnStrandOfChange, forwardStrandReferenceBases, forwardStrandReplacementBases);
            AminoAcidSequence acids = convertBaseSequence(exonBasesAfterChange);
            result.add(new ChangeResult(acids, exonBasesAfterChange, locationOnStrandOfChange, forwardStrandReferenceBases, forwardStrandReplacementBases));
        });
        return result;
    }

    AminoAcidSequence convertBaseSequence(String bases)
    {
        return AminoAcidSequence.fromNucleotides(bases);
    }
}
