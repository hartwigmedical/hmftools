package com.hartwig.hmftools.pave.transval;

import java.util.Objects;

import com.google.common.base.Preconditions;

import org.jetbrains.annotations.NotNull;

public class ExtendedExon
{
    @NotNull private final String prefixFromPreviousExon;
    @NotNull private final String suffixFromNextExon;
    @NotNull private final String exonBases;
    private final int exonStart;

    public ExtendedExon(@NotNull final String prefixFromPreviousExon,
            @NotNull final String suffixFromNextExon,
            @NotNull final String exonBases,
            final int exonStart)
    {
        Preconditions.checkArgument(prefixFromPreviousExon.length() < 3);
        Preconditions.checkArgument(suffixFromNextExon.length() < 3);
        this.prefixFromPreviousExon = prefixFromPreviousExon;
        this.suffixFromNextExon = suffixFromNextExon;
        this.exonBases = exonBases;
        this.exonStart = exonStart;
        Preconditions.checkArgument(totalLength() % 3 == 0);
    }

    public String baseSequenceWithDeletionApplied(int start, int end)
    {
        Preconditions.checkArgument(start >= 0);
        Preconditions.checkArgument(start <= end);
        Preconditions.checkArgument((end - start) % 3 == 2);
        String left = exonBases.substring(0, start);
        String right = exonBases.substring(end + 1);
        return left + right;
    }

    public String basesBetween(int start, int end)
    {
        Preconditions.checkArgument(start >= 0);
        Preconditions.checkArgument(start <= end);
        Preconditions.checkArgument((end - start) % 3 == 2);
        return exonBases.substring(start, end + 1);
    }

    public String completeBaseSequence()
    {
        return prefixFromPreviousExon + exonBases + suffixFromNextExon;
    }

    public int toStrandCoordinates(int positionInExon)
    {
        Preconditions.checkArgument(positionInExon >= 0);
        Preconditions.checkArgument(positionInExon < inExonLength());
        return positionInExon + exonStart;
    }

    public int inExonLength()
    {
        return exonBases.length();
    }

    public int totalLength()
    {
        return prefixFromPreviousExon.length() + suffixFromNextExon.length() + inExonLength();
    }

    @Override
    public boolean equals(final Object o)
    {
        if(o == null || getClass() != o.getClass())
        {
            return false;
        }
        final ExtendedExon that = (ExtendedExon) o;
        return Objects.equals(prefixFromPreviousExon, that.prefixFromPreviousExon)
                && Objects.equals(suffixFromNextExon, that.suffixFromNextExon) && Objects.equals(exonBases, that.exonBases);
    }

    @Override
    public int hashCode()
    {
        return Objects.hash(prefixFromPreviousExon, suffixFromNextExon, exonBases);
    }

    @Override
    public String toString()
    {
        return "ExtendedExon{" +
                "prefixFromPreviousExon='" + prefixFromPreviousExon + '\'' +
                ", suffixFromNextExon='" + suffixFromNextExon + '\'' +
                ", exonBases='" + exonBases + '\'' +
                '}';
    }
}
