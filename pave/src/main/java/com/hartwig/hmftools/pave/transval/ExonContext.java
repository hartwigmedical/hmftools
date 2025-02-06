package com.hartwig.hmftools.pave.transval;

import java.util.Objects;

import com.google.common.base.Preconditions;
import com.hartwig.hmftools.common.codon.Nucleotides;

import org.jetbrains.annotations.NotNull;

public class ExonContext
{
    @NotNull private final String basesOfFirstCodonInPreviousExon;
    @NotNull private final String basesOfLastCodontInFollowingExon;
    @NotNull private final String exonBases;
    private final int exonStart;
    @NotNull private final String intronicPrefix;

    public ExonContext(@NotNull final String basesOfFirstCodonInPreviousExon,
            @NotNull final String basesOfLastCodontInFollowingExon,
            @NotNull final String exonBases,
            final int exonStart,
            @NotNull final String intronicPrefix)
    {
        this.intronicPrefix = intronicPrefix;
        Preconditions.checkArgument(basesOfFirstCodonInPreviousExon.length() < 3);
        Preconditions.checkArgument(basesOfLastCodontInFollowingExon.length() < 3);
        this.basesOfFirstCodonInPreviousExon = basesOfFirstCodonInPreviousExon;
        this.basesOfLastCodontInFollowingExon = basesOfLastCodontInFollowingExon;
        this.exonBases = exonBases;
        this.exonStart = exonStart;
        Preconditions.checkArgument(totalLength() % 3 == 0);
    }

    public String baseSequenceWithDeletionApplied(int start, int end, boolean positiveStrand)
    {
        Preconditions.checkArgument(start >= 0);
        Preconditions.checkArgument(start <= end);
        Preconditions.checkArgument((end - start) % 3 == 2);
        if(!positiveStrand) {
            int exonEnd = exonBases.length() - 1;
            String left = exonBases.substring(0, exonBases.length() - end - 1);
            String right = exonBases.substring(exonEnd - start + 1, exonEnd + 1);
            String complete = basesOfLastCodontInFollowingExon + left + right + basesOfFirstCodonInPreviousExon;
            String result = Nucleotides.reverseComplementBases(complete);
            return result;
//            return Nucleotides.reverseComplementBases(right) + Nucleotides.reverseComplementBases(left);
        }
        String left = exonBases.substring(0, start);
        String right = exonBases.substring(end + 1);

        return basesOfFirstCodonInPreviousExon + left + right + basesOfLastCodontInFollowingExon;
    }

    public String basesBetween(int start, int end)
    {
        Preconditions.checkArgument(start >= 0);
        Preconditions.checkArgument(start <= end);
        Preconditions.checkArgument((end - start) % 3 == 2);
        return exonBases.substring(start, end + 1);
    }

    public int toStrandCoordinates(int positionInExon, boolean positiveStrand)
    {
        Preconditions.checkArgument(positionInExon >= 0);
        Preconditions.checkArgument(positionInExon < inExonLength());
        if(positiveStrand) {
            return positionInExon + exonStart;
        }
        return exonStart + exonBases.length() - positionInExon;
    }

    public int inExonLength()
    {
        return exonBases.length();
    }

    public int totalLength()
    {
        return basesOfFirstCodonInPreviousExon.length() + basesOfLastCodontInFollowingExon.length() + inExonLength();
    }

    @Override
    public boolean equals(final Object o)
    {
        if(o == null || getClass() != o.getClass())
        {
            return false;
        }
        final ExonContext that = (ExonContext) o;
        return Objects.equals(basesOfFirstCodonInPreviousExon, that.basesOfFirstCodonInPreviousExon)
                && Objects.equals(basesOfLastCodontInFollowingExon, that.basesOfLastCodontInFollowingExon) && Objects.equals(exonBases, that.exonBases);
    }

    @Override
    public int hashCode()
    {
        return Objects.hash(basesOfFirstCodonInPreviousExon, basesOfLastCodontInFollowingExon, exonBases);
    }

    @Override
    public String toString()
    {
        return "ExtendedExon{" +
                "prefixFromPreviousExon='" + basesOfFirstCodonInPreviousExon + '\'' +
                ", suffixFromNextExon='" + basesOfLastCodontInFollowingExon + '\'' +
                ", exonBases='" + exonBases + '\'' +
                '}';
    }
}
