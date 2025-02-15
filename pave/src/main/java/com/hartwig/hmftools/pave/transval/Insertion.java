package com.hartwig.hmftools.pave.transval;

import java.util.HashSet;
import java.util.Set;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.base.Preconditions;
import com.hartwig.hmftools.common.codon.Nucleotides;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptAminoAcids;
import com.hartwig.hmftools.common.gene.TranscriptData;

import org.jetbrains.annotations.NotNull;

class Insertion extends ProteinVariant
{
    @NotNull
    final AminoAcidSequence mInsertedSequence;

    public Insertion(@NotNull final GeneData gene,
            @NotNull final TranscriptData transcript,
            @NotNull final TranscriptAminoAcids aminoAcidSequence,
            @NotNull final AminoAcidRange refRange,
            @NotNull final AminoAcidSequence insertedSequence)
    {
        super(gene, transcript, aminoAcidSequence, refRange.startPosition(), refRange.length());
        Preconditions.checkArgument(insertedSequence.length() > 0);
        Preconditions.checkArgument( refRange.length() == 2);
        mInsertedSequence = insertedSequence;
    }

    @VisibleForTesting
    Set<String> possibleInsertedNucleotideSequences()
    {
        final NucleotidesCalculator nucleotidesCalculator = new NucleotidesCalculator(mInsertedSequence, "", "");
        if(mInsertedSequence.length() > 1)
        {
            return Set.of(nucleotidesCalculator.anyBaseSequence());
        }
        return nucleotidesCalculator.allPossibleBaseSequences();
    }

    @NotNull
    @Override
    Set<ChangeResult> applyChange(ChangeContext context)
    {
        int insertionPosition = context.StartPositionInExon;
        int locationOfChange = context.insertionPoint();
        String ref = context.mExon.baseAt(insertionPosition, context.IsPositiveStrand);
        Set<String> baseOptions = possibleInsertedNucleotideSequences();
        Set<ChangeResult> result = new HashSet<>();
        baseOptions.forEach(bases -> {
            String basesToInsert = bases;
            if(!context.IsPositiveStrand)
            {
                basesToInsert = Nucleotides.reverseComplementBases(basesToInsert);
            }
            String withBasesInserted = context.mExon.baseSequenceWithInsertionApplied(insertionPosition, basesToInsert, context.IsPositiveStrand);
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
    AminoAcidSequence variantSequence()
    {
        String rawAAs = this.mAminoAcidSequence.AminoAcids;
        int pointOfInsertion = positionOfFirstAlteredCodon() - 1;
        String left = rawAAs.substring(0, pointOfInsertion);
        String right = rawAAs.substring(pointOfInsertion);
        String duplicatedAAs = left + mInsertedSequence.sequence() + right;
        return AminoAcidSequence.parse(duplicatedAAs);
    }
}
