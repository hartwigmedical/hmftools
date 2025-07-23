package com.hartwig.hmftools.pavereverse.base;

import java.util.Objects;

import com.google.common.base.Preconditions;
import com.hartwig.hmftools.common.codon.Nucleotides;

import org.apache.commons.lang3.tuple.Pair;

public class ChangeContext
{
    public final PaddedExon Exon;
    public final int StartPositionInExon;
    public final int FinishPositionInExon;
    public final boolean IsPositiveStrand;
    private final int AminoAcidNumberOfFirstAminoAcid;
    private ChangeContext mCompanionContext = null;

    public ChangeContext(PaddedExon containingExon, int startPositionInExon,
            int finishPositionInExon,
            boolean isPositiveStrand,
            int aminoAcidNumberOfFirstAminoAcidStartingInExon)
    {
        Exon = containingExon;
        StartPositionInExon = startPositionInExon;
        FinishPositionInExon = finishPositionInExon;
        IsPositiveStrand = isPositiveStrand;
        AminoAcidNumberOfFirstAminoAcid = aminoAcidNumberOfFirstAminoAcidStartingInExon;
    }

    public ChangeContext companionContext()
    {
        return mCompanionContext;
    }

    public void setCompanionContext(ChangeContext companionContext)
    {
        Preconditions.checkArgument(companionContext.Exon.Index == Exon.Index - 1);
        Preconditions.checkArgument(companionContext.IsPositiveStrand == IsPositiveStrand);
        mCompanionContext = companionContext;
    }

    public Pair<String, String> forwardStrandBaseAndLeftNeighbour()
    {
        return Exon.forwardStrandBaseAndLeftNeighbour(StartPositionInExon, IsPositiveStrand);
    }

    public ChangeContext shiftLeft(int i)
    {
        int start = IsPositiveStrand ? StartPositionInExon - i : StartPositionInExon + i;
        int end = IsPositiveStrand ? FinishPositionInExon - i : FinishPositionInExon + i;
        return new ChangeContext(Exon, start, end, IsPositiveStrand, AminoAcidNumberOfFirstAminoAcid);
    }

    public int positionOfChangeStartInStrand()
    {
        if(IsPositiveStrand)
        {
            return positionBeforeStrandPositionFor(StartPositionInExon);
        }
        return positionBeforeStrandPositionFor(FinishPositionInExon + 1);
    }

    public int insertionPoint()
    {
        if(IsPositiveStrand)
        {
            return positionBeforeStrandPositionFor(StartPositionInExon);
        }
        return positionBeforeStrandPositionFor(StartPositionInExon);
    }

    public SplitCodonSequence basesForProteinChange(int firstAminoAcid, int numberOfAminoAcidsChanged)
    {
        Preconditions.checkArgument(firstAminoAcid >= 0);
        Preconditions.checkArgument(numberOfAminoAcidsChanged >= 0);
        int codonNumber = firstAminoAcid - AminoAcidNumberOfFirstAminoAcid + 1;
        return Exon.getSplitSequenceForCodons(codonNumber, numberOfAminoAcidsChanged, IsPositiveStrand);
    }

    public CodonWithinExons codonForProteinChange(int position)
    {
        Preconditions.checkArgument(position >= 0);
        int codonNumber = position - AminoAcidNumberOfFirstAminoAcid + 1;
        return Exon.getCodon(codonNumber, IsPositiveStrand);
    }

    public String exonBasesWithReplacementAppliedAtStrandLocation(int strandLocationOfChange, String currentBases, String alternateBases)
    {
        String result = Exon.baseSequenceWithBasesReplacedAtStrandLocation(strandLocationOfChange, currentBases, alternateBases);
        return IsPositiveStrand ? result : Nucleotides.reverseComplementBases(result);
    }

    public String refBases()
    {
        if(IsPositiveStrand)
        {
            return Exon.baseImmediatelyBefore(StartPositionInExon) + Exon.basesBetween(StartPositionInExon, FinishPositionInExon);
        }
        int end = Exon.inExonLength() - StartPositionInExon - 1;
        int start = Exon.inExonLength() - FinishPositionInExon - 1;
        return Exon.baseImmediatelyBefore(start) + Exon.basesBetween(start, end);
    }

    @Override
    public String toString()
    {
        return "ChangeContext{" +
                "containingExon=" + Exon +
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
                && Objects.equals(Exon, that.Exon);
    }

    @Override
    public int hashCode()
    {
        return Objects.hash(Exon, StartPositionInExon, FinishPositionInExon);
    }

    private int positionBeforeStrandPositionFor(int position)
    {
        return Exon.toStrandCoordinates(position, IsPositiveStrand) - 1;
    }
}
