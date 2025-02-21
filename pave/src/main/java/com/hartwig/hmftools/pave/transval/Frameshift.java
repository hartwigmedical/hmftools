package com.hartwig.hmftools.pave.transval;

import java.util.Set;

import com.google.common.base.Preconditions;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptAminoAcids;
import com.hartwig.hmftools.common.gene.TranscriptData;

import org.apache.commons.lang3.tuple.Pair;
import org.jetbrains.annotations.NotNull;

class Frameshift extends ProteinVariant
{
    @NotNull
    final AminoAcid mFirstChangedAminoAcid;
    public Frameshift(@NotNull final GeneData gene,
            @NotNull final TranscriptData transcript,
            @NotNull final TranscriptAminoAcids aminoAcidSequence,
            @NotNull final AminoAcidRange refRange)
    {
        super(gene, transcript, aminoAcidSequence, refRange.startPosition(), refRange.length());
        Preconditions.checkArgument(refRange.length() == 1);
        mFirstChangedAminoAcid = refRange.aminoAcidAtStart();
    }

    @NotNull
    @Override
    Set<ChangeResult> applyChange(ChangeContext context)
    {
        Pair<String,String> baseToLeftAndBaseDeleted = context.mExon.exonBaseAndImmediatePriorStrandBase(context.StartPositionInExon);
        String newBases = context.mExon.baseSequenceWithSingleBaseRemoved(context.StartPositionInExon);
        AminoAcidSequence newAminoAcids = AminoAcidSequence.fromNucleotides(newBases);
        String alt = baseToLeftAndBaseDeleted.getLeft();
        String ref = alt + baseToLeftAndBaseDeleted.getRight();
        return Set.of(new ChangeResult(newAminoAcids, newBases, context.insertionPoint(), ref, alt));
    }

    @Override
    AminoAcidSequence variantSequence()
    {
        return completeReferenceAminoAcidSequence().subsequenceUpToInclusive(positionOfFirstAlteredCodon() - 1);
    }

    @Override
    boolean isConsistentWithThisVariant(AminoAcidSequence candidate)
    {
        // The sequences need to match up to the change point and then differ
        // at least at that position.
        AminoAcidSequence variant = variantSequence();
        if(candidate.length() < variant.length() + 1)
        {
            return false;
        }
        for (int i = 0; i < variant.length(); i++)
        {
            if(!candidate.get(i).equals(variant.get(i)))
            {
                return false;
            }
        }
        return !candidate.get(variant.length()).equals(mFirstChangedAminoAcid);
    }
}
