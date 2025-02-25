package com.hartwig.hmftools.pave.transval;

import static com.hartwig.hmftools.common.utils.Strings.last;

import java.util.Objects;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.base.Preconditions;
import com.hartwig.hmftools.common.codon.Nucleotides;

import org.apache.commons.lang3.tuple.Pair;
import org.jetbrains.annotations.NotNull;

public class PaddedExon
{
    public final int mIndex;
    @NotNull
    private final String BasesOfFirstCodonInPreviousExon;
    @NotNull
    private final String BasesOfLastCodonInFollowingExon;
    @NotNull
    private final String exonBases;
    private final int exonStart;
    @NotNull
    private final String IntronicPrefix;
    @NotNull
    private final String IntronicSuffix;

    public PaddedExon(
            final int mIndex, @NotNull final String basesOfFirstCodonInPreviousExon,
            @NotNull final String basesOfLastCodontInFollowingExon,
            @NotNull final String exonBases,
            final int exonStart,
            @NotNull final String intronicPrefix,
            @NotNull final String intronicSuffix)
    {
        Preconditions.checkArgument(basesOfFirstCodonInPreviousExon.length() < 3);
        Preconditions.checkArgument(basesOfLastCodontInFollowingExon.length() < 3);
        this.mIndex = mIndex;
        this.IntronicPrefix = intronicPrefix;
        this.IntronicSuffix = intronicSuffix;
        this.BasesOfFirstCodonInPreviousExon = basesOfFirstCodonInPreviousExon;
        this.BasesOfLastCodonInFollowingExon = basesOfLastCodontInFollowingExon;
        this.exonBases = exonBases;
        this.exonStart = exonStart;
        Preconditions.checkArgument(totalLength() % 3 == 0);
    }

    public String baseSequenceWithFramePreservingDeletionApplied(int start, int end, boolean positiveStrand)
    {
        Preconditions.checkArgument(start >= 0);
        Preconditions.checkArgument(start <= end);
        Preconditions.checkArgument((end - start) % 3 == 0);
        if(!positiveStrand)
        {
            String right = exonBases.substring(exonBases.length() - start);
            String left = exonBases.substring(0, exonBases.length() - end);
            String complete = BasesOfFirstCodonInPreviousExon + left + right + BasesOfLastCodonInFollowingExon;
            return Nucleotides.reverseComplementBases(complete);
        }
        String left = exonBases.substring(0, start);
        String right = exonBases.substring(end);

        return BasesOfFirstCodonInPreviousExon + left + right + BasesOfLastCodonInFollowingExon;
    }

    public String baseSequenceWithDuplicationApplied(final int start, final int end, final boolean positiveStrand)
    {
        Preconditions.checkArgument(start >= 0);
        Preconditions.checkArgument(start < end);
        Preconditions.checkArgument((end - start) % 3 == 0);
        if(!positiveStrand)
        {
            int exonEnd = exonBases.length();
            int s = exonEnd - end;
            int e = exonEnd - start;
            String r = exonBases.substring(e);
            String l = exonBases.substring(0, s);
            String m = exonBases.substring(s, e);
            Preconditions.checkArgument((l + m + r).equals(exonBases));
            String complete = BasesOfFirstCodonInPreviousExon + l + m + m + r + BasesOfLastCodonInFollowingExon;
            return Nucleotides.reverseComplementBases(complete);
        }
        String left = exonBases.substring(0, start);
        String middle = exonBases.substring(start, end);
        String right = exonBases.substring(end);
        Preconditions.checkArgument((left + middle + right).equals(exonBases));
        return BasesOfFirstCodonInPreviousExon + left + middle + middle + right + BasesOfLastCodonInFollowingExon;
    }

    public String baseSequenceWithInsertionApplied(final int position, final String basesToInsert, final boolean positiveStrand)
    {
        Preconditions.checkArgument(position >= 0);
        Preconditions.checkArgument(position <= exonBases.length());
        if(!positiveStrand)
        {
            int s = exonBases.length() - position;
            String l = exonBases.substring(0, s);
            String r = exonBases.substring(s);
            Preconditions.checkArgument((l + r).equals(exonBases));
            String complete = BasesOfFirstCodonInPreviousExon + l + basesToInsert + r + BasesOfLastCodonInFollowingExon;
            return Nucleotides.reverseComplementBases(complete);
        }
        String left = exonBases.substring(0, position);
        String right = exonBases.substring(position);
        Preconditions.checkArgument((left + right).equals(exonBases));
        return BasesOfFirstCodonInPreviousExon + left + basesToInsert + right + BasesOfLastCodonInFollowingExon;
    }

    public Pair<String, String> forwardStrandBaseAndLeftNeighbour(final int position, final boolean isPositiveStrand)
    {
        Preconditions.checkArgument(position >= 0);
        Preconditions.checkArgument(position < exonBases.length());
        if(!isPositiveStrand)
        {
            int p = exonBases.length() - position - 1;
            String right = exonBases.substring(p, p + 1);
            String left = p == 0 ? last(IntronicPrefix) : exonBases.substring(p - 1, p);
            return Pair.of(left, right);
        }
        String left = position == 0 ? last(IntronicPrefix) : exonBases.substring(position - 1, position);
        String right = exonBases.substring(position, position + 1);
        return Pair.of(left, right);
    }

    public String baseSequenceWithSingleBaseRemoved(final int position, final boolean isPositiveStrand)
    {
        Preconditions.checkArgument(position >= 0);
        Preconditions.checkArgument(position < exonBases.length());
        if(!isPositiveStrand)
        {
            String right = exonBases.substring(exonBases.length() - position);
            String left = exonBases.substring(0, exonBases.length() - position - 1);
            String leftPadding = last(IntronicPrefix);
            String complete = BasesOfFirstCodonInPreviousExon + leftPadding + left + right + BasesOfLastCodonInFollowingExon;
            return Nucleotides.reverseComplementBases(complete);
        }
        String left = exonBases.substring(0, position);
        String right = position == exonBases.length() - 1 ? "" : exonBases.substring(position + 1);
        // To return a whole number of codons we pad from the suffix.
        // This will allow a sequence of amino acids to be generated from
        // the result, even though it might correspond to lost junction.
        String rightPadding = IntronicSuffix.substring(0, 1);
        return BasesOfFirstCodonInPreviousExon + left + right + rightPadding + BasesOfLastCodonInFollowingExon;
    }

    public String baseSequenceWithBasesReplaced(final int position, final String replacementBases, final boolean positiveStrand)
    {
        Preconditions.checkArgument(position >= 0);
        if(!positiveStrand)
        {
            int s = exonBases.length() - position - replacementBases.length();
            String l = exonBases.substring(0, s);
            String r = exonBases.substring(s + replacementBases.length());
            String complete = BasesOfFirstCodonInPreviousExon + l + replacementBases + r + BasesOfLastCodonInFollowingExon;
            return Nucleotides.reverseComplementBases(complete);
        }
        String left = exonBases.substring(0, position);
        String right = exonBases.substring(position + replacementBases.length());
        return BasesOfFirstCodonInPreviousExon + left + replacementBases + right + BasesOfLastCodonInFollowingExon;
    }

    public String baseSequenceWithBasesReplacedAtStrandLocation(final int location, final String forwardStrandReplacementBases,
            final boolean positiveStrand)
    {
        Preconditions.checkArgument(location >= 0);
        if(!positiveStrand)
        {
            int leftEnd = location - exonStart;
            String l = exonBases.substring(0, leftEnd);
            String r = exonBases.substring(leftEnd + forwardStrandReplacementBases.length());
            String complete = BasesOfFirstCodonInPreviousExon + l + forwardStrandReplacementBases + r + BasesOfLastCodonInFollowingExon;
            return Nucleotides.reverseComplementBases(complete);
        }
        int leftEnd = location - exonStart;

        String left = exonBases.substring(0, leftEnd);
        String right = exonBases.substring(leftEnd + forwardStrandReplacementBases.length());
        return BasesOfFirstCodonInPreviousExon + left + forwardStrandReplacementBases + right + BasesOfLastCodonInFollowingExon;
    }

    public String basesBetween(int start, int end)
    {
        Preconditions.checkArgument(start >= 0);
        Preconditions.checkArgument(start <= end);
        return exonBases.substring(start, end + 1);
    }

    public String baseAt(int start, boolean positiveStrand)
    {
        Preconditions.checkArgument(start >= 0);
        if(!positiveStrand)
        {
            int s = exonBases.length() - start - 1;
            return exonBases.substring(s, s + 1);
        }
        return exonBases.substring(start, start + 1);
    }

    public String baseImmediatelyBefore(int position)
    {
        Preconditions.checkArgument(position >= 0);
        Preconditions.checkArgument(position <= inExonLength());
        if(position == 0)
        {
            return last(IntronicPrefix);
        }
        return exonBases.substring(position - 1, position);
    }

    CodonWithinExons getCodon(int codonIndex, boolean isPositiveStrand)
    {
        Preconditions.checkArgument(codonIndex >= 0);
        if(!isPositiveStrand)
        {
            int start = codonLocationInExonBody(codonIndex, false);
            int stop = start + 3;
            if(start < 0)
            {
                String forwardBasesInExon = exonBases.substring(0, stop);
                String fragmentInExon = Nucleotides.reverseComplementBases(forwardBasesInExon);
                BaseSequence bodySequence = new BaseSequence(exonStart, fragmentInExon, false);
                String suffix = Nucleotides.reverseComplementBases(BasesOfFirstCodonInPreviousExon);
                return CodonWithinExons.factory("", bodySequence, suffix);
            }
            if(stop >= exonBases.length())
            {
                String forwardBasesInExon = exonBases.substring(start);
                String fragmentInExon = Nucleotides.reverseComplementBases(forwardBasesInExon);
                BaseSequence bodySequence = new BaseSequence(start + exonStart, fragmentInExon, false);
                String prefix = Nucleotides.reverseComplementBases(BasesOfLastCodonInFollowingExon);
                return CodonWithinExons.factory(prefix, bodySequence, "");
            }
            String forwardBases = exonBases.substring(start, start + 3);
            String fragmentInExon = Nucleotides.reverseComplementBases(forwardBases);
            BaseSequence bodySequence = new BaseSequence(start + exonStart, fragmentInExon, false);
            return CodonWithinExons.factory("", bodySequence, "");
        }
        int start = codonLocationInExonBody(codonIndex, true);
        int stop = start + 3;
        if(start < 0)
        {
            String fragmentInExon = exonBases.substring(0, stop);
            BaseSequence bodySequence = new BaseSequence(exonStart, fragmentInExon, true);
            return CodonWithinExons.factory(BasesOfFirstCodonInPreviousExon, bodySequence, "");
        }
        if(stop >= exonBases.length())
        {
            String fragmentInExon = exonBases.substring(start);
            BaseSequence bodySequence = new BaseSequence(start + exonStart, fragmentInExon, true);
            return CodonWithinExons.factory("", bodySequence, BasesOfLastCodonInFollowingExon);
        }
        String fragmentInExon = exonBases.substring(start, stop);
        BaseSequence bodySequence = new BaseSequence(start + exonStart, fragmentInExon, true);
        return CodonWithinExons.factory("", bodySequence, "");
    }

    SplitCodonSequence getSplitSequenceForCodons(int startCodon, int count, final boolean isPositiveStrand)
    {
        Preconditions.checkArgument(startCodon >= 0);
        if(startCodon == 0)
        {
            Preconditions.checkArgument(!BasesOfFirstCodonInPreviousExon.isBlank());
        }
        Preconditions.checkArgument(count > 0);
        if(!isPositiveStrand)
        {
            int stop = codonLocationInExonBody(startCodon, false) + 3; //end of deletion is the codon end.
            int start = stop - (3 * count);
            if(start < 0)
            {
                String forwardBasesInExon = exonBases.substring(0, stop);
                String reversed = Nucleotides.reverseComplementBases(forwardBasesInExon);
                String suffix = Nucleotides.reverseComplementBases(BasesOfFirstCodonInPreviousExon);
                return new SplitCodonSequence(reversed, suffix, exonStart);
            }
            if(stop >= exonBases.length())
            {
                String forwardBasesInExon = exonBases.substring(start);
                String reversed = Nucleotides.reverseComplementBases(forwardBasesInExon);
                String prefix = Nucleotides.reverseComplementBases(BasesOfLastCodonInFollowingExon);
                return new SplitCodonSequence(prefix, reversed, exonStart + start);
            }
            String forwardBases = exonBases.substring(start, stop);
            String reversed = Nucleotides.reverseComplementBases(forwardBases);
            return new SplitCodonSequence(reversed, "", exonStart + start);
        }
        int start = codonLocationInExonBody(startCodon, true);
        int stop = start + count * 3;

        if(start < 0)
        {
            return new SplitCodonSequence(BasesOfFirstCodonInPreviousExon, exonBases.substring(0, stop), exonStart);
        }
        if(stop >= exonBases.length())
        {
            return new SplitCodonSequence(exonBases.substring(start), BasesOfLastCodonInFollowingExon, start + exonStart);
        }
        String fragmentInExon = exonBases.substring(start, stop);
        return new SplitCodonSequence(fragmentInExon, "", start + exonStart);
    }

    public int codonLocationInExonBody(int codonIndex, final boolean isPositiveStrand)
    {
        Preconditions.checkArgument(codonIndex >= 0);

        if(!isPositiveStrand)
        {
            int b = exonBases.length();
            int numberUsedInExon0 = (3 - numberOfBasesInFollowingExon()) % 3;
            int startOfExon0 = b - numberUsedInExon0;
            return startOfExon0 - 3 * codonIndex;
        }
        int offset = (3 - numberOfBasesInPreviousExon()) % 3;
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
        if(positiveStrand)
        {
            return positionInExon + exonStart;
        }
        return exonStart + exonBases.length() - positionInExon;
    }

    public int fromStrandCoordinates(int positionInStrand, boolean positiveStrand)
    {
        Preconditions.checkArgument(positionInStrand >= 0);
        if(positiveStrand)
        {
            return positionInStrand - exonStart;
        }
        return exonStart + exonBases.length() - positionInStrand;
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
                && Objects.equals(BasesOfLastCodonInFollowingExon, that.BasesOfLastCodonInFollowingExon)
                && Objects.equals(exonBases, that.exonBases);
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
