package com.hartwig.hmftools.pave.transval;

import static com.hartwig.hmftools.common.codon.Nucleotides.reverseComplementBases;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.SortedSet;

import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptAminoAcids;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;

import org.apache.commons.lang3.tuple.Pair;
import org.jetbrains.annotations.NotNull;

class SingleAminoAcidVariant extends ProteinVariant
{
    @NotNull
    private final AminoAcidSpecification Alt;
    @NotNull
    private final CodonRegions RegionsDefiningCodon;

    public SingleAminoAcidVariant(
            @NotNull final GeneData gene,
            @NotNull final TranscriptData transcript,
            @NotNull final TranscriptAminoAcids aminoAcidSequence,
            @NotNull final AminoAcidSpecification alt)
    {
        super(gene, transcript, aminoAcidSequence, alt.position, 1);
        this.Alt = alt;
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
        if(!mGene.forwardStrand())
        {
            nucleotideDifferences = Pair.of(reverseComplementBases(nucleotideDifferences.getLeft()), reverseComplementBases(nucleotideDifferences.getRight()));
        }
        List<String> alternateNucleotides = new ArrayList<>();
        codonChanges.forEach(cv -> alternateNucleotides.add(cv.AlternateCodon));
        Set<TransvalHotspot> hotspots = new HashSet<>();
        codonChanges.forEach(cv -> {
            Pair<String,String> refAlt = cv.differenceStrings();
            if (mGene.forwardStrand())
            {
                int position = RegionsDefiningCodon.translateCodonPosition(cv.positionOfFirstDifference());
                hotspots.add(new TransvalHotspot(refAlt.getLeft(), refAlt.getRight(), mGene.Chromosome, position));
            } else {
                int position = RegionsDefiningCodon.translateCodonPosition(cv.positionOfFirstDifference()) - refAlt.getLeft().length() + 1;
                hotspots.add(new TransvalHotspot(
                        reverseComplementBases(refAlt.getLeft()),
                        reverseComplementBases(refAlt.getRight()),
                        mGene.Chromosome, position));
            }
        });
        return new TransvalSnvMnv(
                mTranscript,
                mGene.Chromosome,
                positionOfChangeInChromosome,
                !codonIsInSingleExon(),
                nucleotideDifferences.getLeft(),
                hotspots,
                nucleotideDifferences.getRight(),
                referenceCodon,
                alternateNucleotides
        );
    }

    public String altValue()
    {
        return Alt.symbol();
    }

    public String referenceCodon(RefGenomeInterface refGenomeSource)
    {
        return RegionsDefiningCodon.retrieveCodon(refGenomeSource);
    }

    public boolean codonIsInSingleExon()
    {
        return RegionsDefiningCodon.codonIsInSingleExon();
    }
}
