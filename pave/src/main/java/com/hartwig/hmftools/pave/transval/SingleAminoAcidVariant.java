package com.hartwig.hmftools.pave.transval;

import static com.hartwig.hmftools.common.codon.Nucleotides.reverseComplementBases;

import java.util.HashSet;
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
        super(gene, transcript, aminoAcidSequence, alt.mPosition, 1);
        this.Alt = alt;
        int codonPosition = 3 * (positionOfFirstAlteredCodon() - 1);
        RegionsDefiningCodon = exonsForCodonPosition(codonPosition);
    }

    @Override
    AminoAcidSequence variantSequence()
    {
        return null;
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
            int locationOnStrandOfChange = codon.strandLocationOfChange(codonChange.positionOfFirstDifference());
            Pair<String, String> refAlt = codonChange.differenceStrings();
            String exonBasesAfterChange = context.exonBasesWithReplacementApplied(locationOnStrandOfChange, refAlt.getRight());
            AminoAcidSequence acids = AminoAcidSequence.fromNucleotides(exonBasesAfterChange);
            result.add( new ChangeResult(acids, exonBasesAfterChange, locationOnStrandOfChange, refAlt.getLeft(), refAlt.getRight()));
        });
        return result;
    }

    public TransvalVariant calculateVariant2(RefGenomeInterface refGenome)
    {
        String referenceCodon = referenceCodon(refGenome);
        CodonVariantsCalculator change = new CodonVariantsCalculator(referenceAminoAcids(), altValue());
        SortedSet<CodonChange> codonChanges = change.possibleCodonsForVariant(referenceCodon);
        if(codonChanges.isEmpty())
        {
            return null;
        }
        Set<TransvalHotspot> hotspots = new HashSet<>();
        codonChanges.forEach(cv ->
        {
            Pair<String, String> refAlt = cv.differenceStrings();
            if(mGene.forwardStrand())
            {
                int position = RegionsDefiningCodon.translateCodonPosition(cv.positionOfFirstDifference());
                hotspots.add(new TransvalHotspot(refAlt.getLeft(), refAlt.getRight(), mGene.Chromosome, position));
            }
            else
            {
                int position = RegionsDefiningCodon.translateCodonPosition(cv.positionOfFirstDifference()) - refAlt.getLeft().length() + 1;
                hotspots.add(new TransvalHotspot(
                        reverseComplementBases(refAlt.getLeft()),
                        reverseComplementBases(refAlt.getRight()),
                        mGene.Chromosome, position));
            }
        });
        return new TransvalVariant(
                mTranscript,
                mGene.Chromosome,
                !codonIsInSingleExon(),
                hotspots
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
