package com.hartwig.hmftools.pavereverse.variants;

import java.util.HashSet;
import java.util.Set;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.base.Preconditions;
import com.hartwig.hmftools.common.codon.Nucleotides;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptAminoAcids;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.pavereverse.aa.AminoAcidRange;
import com.hartwig.hmftools.pavereverse.aa.AminoAcidSequence;
import com.hartwig.hmftools.pavereverse.base.ChangeContext;

public class Insertion extends ProteinVariant
{
    public final AminoAcidSequence InsertedSequence;

    public Insertion(GeneData gene,
            TranscriptData transcript,
            TranscriptAminoAcids aminoAcidSequence,
            AminoAcidRange refRange,
            AminoAcidSequence insertedSequence)
    {
        super(gene, transcript, aminoAcidSequence, refRange.startPosition(), refRange.length());
        Preconditions.checkArgument(insertedSequence.length() > 0);
        Preconditions.checkArgument(refRange.length() == 2);
        InsertedSequence = insertedSequence;
    }

    @VisibleForTesting
    Set<String> possibleInsertedNucleotideSequences()
    {
        final NucleotidesCalculator nucleotidesCalculator = new NucleotidesCalculator(InsertedSequence, "", "");
        if(InsertedSequence.length() > 1)
        {
            return Set.of(nucleotidesCalculator.anyBaseSequence());
        }
        return nucleotidesCalculator.allPossibleBaseSequences();
    }

    @Override
    Set<ChangeResult> applyChange(ChangeContext context)
    {
        int insertionPosition = context.StartPositionInExon;
        int locationOfChange = context.insertionPoint();
        String ref = context.IsPositiveStrand
                ? context.forwardStrandBaseAndLeftNeighbour().getLeft()
                : context.forwardStrandBaseAndLeftNeighbour().getRight();
        Set<String> baseOptions = possibleInsertedNucleotideSequences();
        Set<ChangeResult> result = new HashSet<>();
        baseOptions.forEach(bases ->
        {
            String basesToInsert = bases;
            if(!context.IsPositiveStrand)
            {
                basesToInsert = Nucleotides.reverseComplementBases(basesToInsert);
            }
            String withBasesInserted =
                    context.Exon.baseSequenceWithInsertionApplied(insertionPosition, basesToInsert, context.IsPositiveStrand);
            AminoAcidSequence acids = AminoAcidSequence.fromNucleotides(withBasesInserted);
            String alt = ref + basesToInsert;
            result.add(new ChangeResult(acids, withBasesInserted, locationOfChange, ref, alt));
        });

        return result;
    }

    @Override
    public int positionOfFirstAlteredCodon()
    {
        // The insertion happens just after the first AA in the ref seq.
        return super.positionOfFirstAlteredCodon() + 1;
    }

    @Override
    int numberOfLeftShiftsToTry(ChangeContext changeContext)
    {
        return super.numberOfLeftShiftsToTry(changeContext) - 1;
    }

    @Override
    AminoAcidSequence variantSequence()
    {
        int pointOfInsertion = super.positionOfFirstAlteredCodon();
        return completeReferenceAminoAcidSequence().insert(pointOfInsertion, InsertedSequence);
    }
}
