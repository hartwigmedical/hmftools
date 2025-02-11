package com.hartwig.hmftools.pave.transval;

import java.util.Objects;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.base.Preconditions;
import com.hartwig.hmftools.common.codon.Nucleotides;

import org.jetbrains.annotations.NotNull;

public class PaddedExon
{
    @NotNull private final String BasesOfFirstCodonInPreviousExon;
    @NotNull private final String BasesOfLastCodonInFollowingExon;
    @NotNull private final String exonBases;
    private final int exonStart;
    @NotNull private final String IntronicPrefix;

    public PaddedExon(@NotNull final String basesOfFirstCodonInPreviousExon,
            @NotNull final String basesOfLastCodontInFollowingExon,
            @NotNull final String exonBases,
            final int exonStart,
            @NotNull final String intronicPrefix)
    {
        this.IntronicPrefix = intronicPrefix;
        Preconditions.checkArgument(basesOfFirstCodonInPreviousExon.length() < 3);
        Preconditions.checkArgument(basesOfLastCodontInFollowingExon.length() < 3);
        this.BasesOfFirstCodonInPreviousExon = basesOfFirstCodonInPreviousExon;
        this.BasesOfLastCodonInFollowingExon = basesOfLastCodontInFollowingExon;
        this.exonBases = exonBases;
        this.exonStart = exonStart;
        Preconditions.checkArgument(totalLength() % 3 == 0);
    }

    public String baseSequenceWithDeletionApplied(int start, int end, boolean positiveStrand)
    {
        Preconditions.checkArgument(start >= 0);
        Preconditions.checkArgument(start <= end);
        Preconditions.checkArgument((end - start + 1) % 3 == 0);
        if(!positiveStrand)
        {
            int exonEnd = exonBases.length() - 1;
            String left = exonBases.substring(0, exonBases.length() - end - 1);
            String right = exonBases.substring(exonEnd - start + 1, exonEnd + 1);
            String complete = BasesOfLastCodonInFollowingExon + left + right + BasesOfFirstCodonInPreviousExon;
            return Nucleotides.reverseComplementBases(complete);
        }
        String left = exonBases.substring(0, start);
        String right = exonBases.substring(end + 1);

        return BasesOfFirstCodonInPreviousExon + left + right + BasesOfLastCodonInFollowingExon;
    }

    public String basesBetween(int start, int end)
    {
        Preconditions.checkArgument(start >= 0);
        Preconditions.checkArgument(start <= end);
        return exonBases.substring(start, end + 1);
    }

    public SplitCodonSequence getSplitSequenceForCodons(int startCodon, int count, final boolean isPositiveStrand)
    {
        Preconditions.checkArgument(startCodon >= 0);
        if(startCodon == 0)
        {
            Preconditions.checkArgument(!BasesOfFirstCodonInPreviousExon.isBlank());
        }
        Preconditions.checkArgument(count > 0);
        if(!isPositiveStrand)
        {
            if(numberOfBasesInPreviousExon() + numberOfBasesInFollowingExon() > 0)
            {
                throw new IllegalArgumentException("Not yet implemented");
            }
            int deletionEnd = exonBases.length() - 3 * (startCodon - 1);
            int deletionStart = deletionEnd - 3 * count;
            String fragmentInExon = exonBases.substring(deletionStart, deletionEnd);
            String reversed = Nucleotides.reverseComplementBases(fragmentInExon);
            return new SplitCodonSequence(reversed, "", exonStart + deletionEnd - 1);
        }
        int start = codonLocationInExonBody(startCodon, isPositiveStrand);
        int stop = start + count * 3;

        if(start < 0) {
            return new SplitCodonSequence(BasesOfFirstCodonInPreviousExon, exonBases.substring(0, stop), exonStart);
        }
        if(stop >= exonBases.length()) {
            return new SplitCodonSequence(exonBases.substring(start), BasesOfLastCodonInFollowingExon, start + exonStart);
        }
        String fragmentInExon = exonBases.substring(start, stop);
        return new SplitCodonSequence(fragmentInExon, "", start + exonStart);
    }

    public int codonLocationInExonBody(int codonIndex, final boolean isPositiveStrand)
    {
        Preconditions.checkArgument(codonIndex >= 0);
        int offset = (3 - numberOfBasesInPreviousExon()) % 3;

        if(!isPositiveStrand)
        {
            return exonBases.length() - 3 * (codonIndex + 1);
//            int offset = (3 - numberOfBasesInFollowingExon()) %3;
        }
        return 3 * (codonIndex - 1) + offset;
    }

    @VisibleForTesting
    public int numberOfBasesInPreviousExon()
    {
        return BasesOfFirstCodonInPreviousExon.length();
    }

    @VisibleForTesting
    public int numberOfBasesInFollowingExon()
    {
        return BasesOfLastCodonInFollowingExon.length();
    }

    public int toStrandCoordinates(int positionInExon, boolean positiveStrand)
    {
        Preconditions.checkArgument(positionInExon >= 0);
        Preconditions.checkArgument(positionInExon <= inExonLength());
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
        return BasesOfFirstCodonInPreviousExon.length() + BasesOfLastCodonInFollowingExon.length() + inExonLength();
    }

    @Override
    public boolean equals(final Object o)
    {
        if(o == null || getClass() != o.getClass())
        {
            return false;
        }
        final PaddedExon that = (PaddedExon) o;
        return Objects.equals(BasesOfFirstCodonInPreviousExon, that.BasesOfFirstCodonInPreviousExon)
                && Objects.equals(BasesOfLastCodonInFollowingExon, that.BasesOfLastCodonInFollowingExon) && Objects.equals(exonBases, that.exonBases);
    }

    @Override
    public int hashCode()
    {
        return Objects.hash(BasesOfFirstCodonInPreviousExon, BasesOfLastCodonInFollowingExon, exonBases);
    }

    @Override
    public String toString()
    {
        return "ExtendedExon{" +
                "prefixFromPreviousExon='" + BasesOfFirstCodonInPreviousExon + '\'' +
                ", suffixFromNextExon='" + BasesOfLastCodonInFollowingExon + '\'' +
                ", exonBases='" + exonBases + '\'' +
                '}';
    }
}
