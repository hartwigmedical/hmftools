package com.hartwig.hmftools.pave.transval;

import java.util.Objects;

import com.hartwig.hmftools.common.codon.Nucleotides;

import org.jetbrains.annotations.NotNull;

public class ChangeContext
{
    @NotNull
    final ExtendedExon containingExon;
    final int startPositionInExon;
    final int finishPositionInExon;
    final boolean IsPositiveStrand;

    public ChangeContext(@NotNull final ExtendedExon containingExon, final int startPositionInExon,
            final int finishPositionInExon,
            final boolean isPositiveStrand)
    {
        this.containingExon = containingExon;
        this.startPositionInExon = startPositionInExon;
        this.finishPositionInExon = finishPositionInExon;
        this.IsPositiveStrand = isPositiveStrand;
    }

    int positionOfChangeStartInStrand()
    {
        if (IsPositiveStrand) {
            return containingExon.toStrandCoordinates(startPositionInExon, IsPositiveStrand);
        }
        return containingExon.toStrandCoordinates(finishPositionInExon + 1, IsPositiveStrand);
    }

    AminoAcidSequence applyDeletion()
    {
        String exonBasesAfterDeletion = containingExon.baseSequenceWithDeletionApplied(startPositionInExon, finishPositionInExon, IsPositiveStrand);
        return AminoAcidSequence.fromNucleotides(exonBasesAfterDeletion);
    }

    public TransvalHotspot hotspot(String chromosome)
    {
//        if (IsPositiveStrand)
//        {
            return new TransvalHotspot(affectedBases(), "", chromosome, positionOfChangeStartInStrand());
//        }
//        int location = positionOfChangeStartInStrand();
//        return new TransvalHotspot(affectedBases(), "", chromosome, location);
    }

    public String affectedBases()
    {
        if (IsPositiveStrand)
        {
            return containingExon.basesBetween(startPositionInExon, finishPositionInExon);
        }
        int actualEnd = containingExon.inExonLength() - startPositionInExon - 1;
        int actualStart = containingExon.inExonLength() - finishPositionInExon - 1;
        final String positiveStrandBases = containingExon.basesBetween(actualStart, actualEnd);
//        return Nucleotides.reverseComplementBases(positiveStrandBases);
        return positiveStrandBases;
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
