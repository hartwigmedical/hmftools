package com.hartwig.hmftools.pave.transval;

import java.util.Set;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.base.Preconditions;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptAminoAcids;
import com.hartwig.hmftools.common.gene.TranscriptData;

import org.jetbrains.annotations.NotNull;

class Insertion extends ProteinVariant
{
    @NotNull
    private final AminoAcidSequence mInsertedSequence;

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
        Set<String> allPossibilities = new NucleotidesCalculator(mInsertedSequence, "", "").possibilities();
        if(mInsertedSequence.length() > 1)
        {
            return Set.of(allPossibilities.iterator().next());
        }
        return allPossibilities;
    }

    @Override
    ChangeResult applyChange(ChangeContext changeContext)
    {
        return changeContext.applyDuplication();
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

    @Override
    TransvalHotspot convertToHotspot(final ChangeContext changeContext)
    {
        String duplicated = changeContext.refBases();
        String refBase = duplicated.substring(0,  1);
        return new TransvalHotspot(refBase, duplicated, mGene.Chromosome, changeContext.positionOfChangeStartInStrand() - 1);
    }
}
