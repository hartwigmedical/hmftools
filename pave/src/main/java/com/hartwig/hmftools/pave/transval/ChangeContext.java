package com.hartwig.hmftools.pave.transval;

import java.util.Objects;

import com.google.common.base.Preconditions;

import org.jetbrains.annotations.NotNull;

public class ChangeContext
{
    @NotNull
    final PaddedExon ContainingExon;
    final int StartPositionInExon;
    final int FinishPositionInExon;
    final boolean IsPositiveStrand;
    private final int AminoAcidNumberOfFirstAminoAcid;

    public ChangeContext(@NotNull final PaddedExon containingExon, final int startPositionInExon,
            final int finishPositionInExon,
            final boolean isPositiveStrand,
            int aminoAcidNumberOfFirstAminoAcidStartingInExon)
    {
        this.ContainingExon = containingExon;
        this.StartPositionInExon = startPositionInExon;
        this.FinishPositionInExon = finishPositionInExon;
        this.IsPositiveStrand = isPositiveStrand;
        AminoAcidNumberOfFirstAminoAcid = aminoAcidNumberOfFirstAminoAcidStartingInExon;
    }

    int positionOfChangeStartInStrand()
    {
        if(IsPositiveStrand)
        {
            return ContainingExon.toStrandCoordinates(StartPositionInExon, IsPositiveStrand);
        }
        return ContainingExon.toStrandCoordinates(FinishPositionInExon + 1, IsPositiveStrand);
    }

    String affectedBases()
    {
        return ContainingExon.basesBetween(StartPositionInExon, FinishPositionInExon);
    }

    String baseImmediatelyBeforeChange()
    {
        return ContainingExon.baseImmediatelyBefore(StartPositionInExon);
    }

    ChangeResult applyDuplication()
    {
        String bases = ContainingExon.baseSequenceWithDuplicationApplied(StartPositionInExon, FinishPositionInExon + 1, IsPositiveStrand);
        return new ChangeResult(AminoAcidSequence.fromNucleotides(bases), bases);
    }

    AminoAcidSequence applyDeletion()
    {
        return AminoAcidSequence.fromNucleotides(exonBasesAfterDeletion());
    }

    String exonBasesAfterDeletion()
    {
        return ContainingExon.baseSequenceWithDeletionApplied(StartPositionInExon, FinishPositionInExon, IsPositiveStrand);
    }

    public SplitCodonSequence basesForProteinChange(int firstAminoAcid, int numberOfAminoAcidsChanged, boolean isPositiveStrand)
    {
        Preconditions.checkArgument(firstAminoAcid >= 0);
        Preconditions.checkArgument(numberOfAminoAcidsChanged >= 0);
        int codonNumber = firstAminoAcid - AminoAcidNumberOfFirstAminoAcid + 1;
        return ContainingExon.getSplitSequenceForCodons(codonNumber, numberOfAminoAcidsChanged, isPositiveStrand);
    }

    public String refBases()
    {
        if(IsPositiveStrand)
        {
            return ContainingExon.basesBetween(StartPositionInExon - 1, FinishPositionInExon);
        }
        int excisionEnd = ContainingExon.inExonLength() - StartPositionInExon - 1;
        int excisionStart = ContainingExon.inExonLength() - FinishPositionInExon - 1;
        return ContainingExon.basesBetween(excisionStart - 1, excisionEnd);
    }

    @Override
    public String toString()
    {
        return "ChangeContext{" +
                "containingExon=" + ContainingExon +
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
                && Objects.equals(ContainingExon, that.ContainingExon);
    }

    @Override
    public int hashCode()
    {
        return Objects.hash(ContainingExon, StartPositionInExon, FinishPositionInExon);
    }
}
