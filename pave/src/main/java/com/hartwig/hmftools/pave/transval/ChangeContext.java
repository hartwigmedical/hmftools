package com.hartwig.hmftools.pave.transval;

import java.util.Objects;

import com.google.common.base.Preconditions;

import org.jetbrains.annotations.NotNull;

public class ChangeContext
{
    @NotNull
    final PaddedExon mExon;
    final int StartPositionInExon;
    final int FinishPositionInExon;
    final boolean IsPositiveStrand;
    private final int AminoAcidNumberOfFirstAminoAcid;

    public ChangeContext(@NotNull final PaddedExon containingExon, final int startPositionInExon,
            final int finishPositionInExon,
            final boolean isPositiveStrand,
            int aminoAcidNumberOfFirstAminoAcidStartingInExon)
    {
        this.mExon = containingExon;
        this.StartPositionInExon = startPositionInExon;
        this.FinishPositionInExon = finishPositionInExon;
        this.IsPositiveStrand = isPositiveStrand;
        AminoAcidNumberOfFirstAminoAcid = aminoAcidNumberOfFirstAminoAcidStartingInExon;
    }

    int positionOfChangeStartInStrand()
    {
        if(IsPositiveStrand)
        {
            return mExon.toStrandCoordinates(StartPositionInExon, IsPositiveStrand);
        }
        return mExon.toStrandCoordinates(FinishPositionInExon + 1, IsPositiveStrand);
    }

    int insertionPoint()
    {
        if(IsPositiveStrand)
        {
            return mExon.toStrandCoordinates(StartPositionInExon, IsPositiveStrand) - 1;
        }
        return mExon.toStrandCoordinates(StartPositionInExon, IsPositiveStrand) - 1;
    }

    public SplitCodonSequence basesForProteinChange(int firstAminoAcid, int numberOfAminoAcidsChanged, boolean isPositiveStrand)
    {
        Preconditions.checkArgument(firstAminoAcid >= 0);
        Preconditions.checkArgument(numberOfAminoAcidsChanged >= 0);
        int codonNumber = firstAminoAcid - AminoAcidNumberOfFirstAminoAcid + 1;
        return mExon.getSplitSequenceForCodons(codonNumber, numberOfAminoAcidsChanged, isPositiveStrand);
    }

    public String refBases()
    {
        if(IsPositiveStrand)
        {
            return mExon.baseImmediatelyBefore(StartPositionInExon) + mExon.basesBetween(StartPositionInExon, FinishPositionInExon);
        }
        int end = mExon.inExonLength() - StartPositionInExon - 1;
        int start = mExon.inExonLength() - FinishPositionInExon - 1;
        return mExon.baseImmediatelyBefore(start) + mExon.basesBetween(start, end);
    }

    @Override
    public String toString()
    {
        return "ChangeContext{" +
                "containingExon=" + mExon +
                ", startPositionInExon=" + StartPositionInExon +
                ", finishPositionInExon=" + FinishPositionInExon +
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
        return StartPositionInExon == that.StartPositionInExon && FinishPositionInExon == that.FinishPositionInExon
                && Objects.equals(mExon, that.mExon);
    }

    @Override
    public int hashCode()
    {
        return Objects.hash(mExon, StartPositionInExon, FinishPositionInExon);
    }
}
