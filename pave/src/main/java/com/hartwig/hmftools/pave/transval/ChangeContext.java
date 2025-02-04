package com.hartwig.hmftools.pave.transval;

import java.util.Objects;

import org.jetbrains.annotations.NotNull;

public class ChangeContext
{
    @NotNull
    final ExtendedExon containingExon;
    final int startPositionInExon;
    final int finishPositionInExon;

    public ChangeContext(@NotNull final ExtendedExon containingExon, final int startPositionInExon, final int finishPositionInExon)
    {
        this.containingExon = containingExon;
        this.startPositionInExon = startPositionInExon;
        this.finishPositionInExon = finishPositionInExon;
    }

    int positionOfChangeStartInStrand()
    {
        return containingExon.toStrandCoordinates(startPositionInExon);
    }

    AminoAcidSequence applyDeletion()
    {
        String exonBasesAfterDeletion = containingExon.baseSequenceWithDeletionApplied(startPositionInExon, finishPositionInExon);
        return AminoAcidSequence.fromNucleotides(exonBasesAfterDeletion);
    }

    public String affectedBases()
    {
        return containingExon.basesBetween(startPositionInExon, finishPositionInExon);
    }

    @Override
    public String toString()
    {
        return "ChangeContext{" +
                "containingExon=" + containingExon +
                ", startPositionInExon=" + startPositionInExon +
                ", finishPositionInExon=" + finishPositionInExon +
                '}';
    }

    @Override
    public boolean equals(final Object o)
    {
        if(o == null || getClass() != o.getClass())
        {
            return false;
        }
        final ChangeContext that = (ChangeContext) o;
        return startPositionInExon == that.startPositionInExon && finishPositionInExon == that.finishPositionInExon
                && Objects.equals(containingExon, that.containingExon);
    }

    @Override
    public int hashCode()
    {
        return Objects.hash(containingExon, startPositionInExon, finishPositionInExon);
    }
}
