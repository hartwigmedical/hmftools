package com.hartwig.hmftools.pave.transval;

import static com.hartwig.hmftools.common.codon.Nucleotides.reverseComplementBases;

import java.util.ArrayList;
import java.util.List;
import java.util.SortedSet;

import com.google.common.base.Preconditions;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptAminoAcids;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;

import org.apache.commons.lang3.tuple.Pair;
import org.jetbrains.annotations.NotNull;

public class SingleAminoAcidVariant extends ProteinVariant
{
    @NotNull
    private final String Alt;
    @NotNull
    private final CodonRegions RegionsDefiningCodon;

    public SingleAminoAcidVariant(
            @NotNull final GeneData gene,
            @NotNull final TranscriptData transcript,
            @NotNull final TranscriptAminoAcids aminoAcidSequence,
            final int position,
            @NotNull final String variant)
    {
        super(gene, transcript, aminoAcidSequence, position);
        Preconditions.checkArgument(Checks.isValidAminoAcidName(variant));
        this.Alt = variant;
        int codonPosition = 3 * (positionOfFirstAlteredCodon() - 1);
        RegionsDefiningCodon = exonsForCodonPosition(codonPosition);
    }

    @Override
    public TransvalSnvMnv calculateVariant(RefGenomeInterface refGenome)
    {
        String referenceCodon = referenceCodon(refGenome);
        CodonVariantsCalculator change = new CodonVariantsCalculator(referenceAminoAcids(), altValue());
        SortedSet<CodonVariant> codonChanges = change.possibleCodonsForVariant(referenceCodon);
        if(codonChanges.isEmpty())
        {
            return null;
        }
        CodonVariant codonVariant = codonChanges.first();
        int positionInCodonOfChange = codonVariant.positionOfFirstDifference();
        int positionOfChangeInChromosome = RegionsDefiningCodon.translateCodonPosition(positionInCodonOfChange);
        Pair<String,String> nucleotideDifferences = codonVariant.differenceStrings();
        if(!Gene.forwardStrand())
        {
            nucleotideDifferences = Pair.of(reverseComplementBases(nucleotideDifferences.getLeft()), reverseComplementBases(nucleotideDifferences.getRight()));
        }
        List<String> alternateNucleotides = new ArrayList<>();
        codonChanges.forEach(cv -> alternateNucleotides.add(cv.AlternateCodon));
        return new TransvalSnvMnv(
                Transcript.TransName,
                Gene.Chromosome,
                positionOfChangeInChromosome,
                !codonIsInSingleExon(),
                nucleotideDifferences.getLeft(),
                nucleotideDifferences.getRight(),
                referenceCodon,
                alternateNucleotides
        );
    }

    public String altValue()
    {
        return Alt;
    }

    public String referenceCodon(RefGenomeInterface refGenomeSource)
    {
        return RegionsDefiningCodon.retrieveCodon(refGenomeSource);
    }

    public boolean codonIsInSingleExon()
    {
        return RegionsDefiningCodon.codonIsInSingleExon();
    }

    @Override
    int changedReferenceSequenceLength()
    {
        return 1;
    }
}
