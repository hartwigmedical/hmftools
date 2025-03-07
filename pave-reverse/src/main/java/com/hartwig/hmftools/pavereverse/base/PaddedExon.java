package com.hartwig.hmftools.pavereverse.base;

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
    private final String mBasesOfFirstCodonInPreviousExon;
    @NotNull
    private final String mBasesOfLastCodonInFollowingExon;
    @NotNull
    private final String mExonBases;
    private final int mExonStart;
    @NotNull
    private final String mIntronicPrefix;
    @NotNull
    private final String mIntronicSuffix;

    public PaddedExon(
            final int mIndex, @NotNull final String basesOfFirstCodonInPreviousExon,
            @NotNull final String basesOfLastCodonInFollowingExon,
            @NotNull final String exonBases,
            final int exonStart,
            @NotNull final String intronicPrefix,
            @NotNull final String intronicSuffix)
    {
        Preconditions.checkArgument(basesOfFirstCodonInPreviousExon.length() < 3);
        Preconditions.checkArgument(basesOfLastCodonInFollowingExon.length() < 3);
        this.mIndex = mIndex;
        this.mIntronicPrefix = intronicPrefix;
        this.mIntronicSuffix = intronicSuffix;
        this.mBasesOfFirstCodonInPreviousExon = basesOfFirstCodonInPreviousExon;
        this.mBasesOfLastCodonInFollowingExon = basesOfLastCodonInFollowingExon;
        this.mExonBases = exonBases;
        this.mExonStart = exonStart;
        Preconditions.checkArgument(totalLength() % 3 == 0);
    }

    public String baseSequenceWithFramePreservingDeletionApplied(int start, int end, boolean positiveStrand)
    {
        Preconditions.checkArgument(start >= 0);
        Preconditions.checkArgument(start <= end);
        Preconditions.checkArgument((end - start) % 3 == 0);
        if(!positiveStrand)
        {
            String right = mExonBases.substring(mExonBases.length() - start);
            String left = mExonBases.substring(0, mExonBases.length() - end);
            String complete = mBasesOfFirstCodonInPreviousExon + left + right + mBasesOfLastCodonInFollowingExon;
            return Nucleotides.reverseComplementBases(complete);
        }
        String left = mExonBases.substring(0, start);
        String right = mExonBases.substring(end);

        return mBasesOfFirstCodonInPreviousExon + left + right + mBasesOfLastCodonInFollowingExon;
    }

    public String baseSequenceWithDuplicationApplied(final int start, final int end, final boolean positiveStrand)
    {
        Preconditions.checkArgument(start >= 0);
        Preconditions.checkArgument(start < end);
        Preconditions.checkArgument((end - start) % 3 == 0);
        if(!positiveStrand)
        {
            int exonEnd = mExonBases.length();
            int s = exonEnd - end;
            int e = exonEnd - start;
            String r = mExonBases.substring(e);
            String l = mExonBases.substring(0, s);
            String m = mExonBases.substring(s, e);
            Preconditions.checkArgument((l + m + r).equals(mExonBases));
            String complete = mBasesOfFirstCodonInPreviousExon + l + m + m + r + mBasesOfLastCodonInFollowingExon;
            return Nucleotides.reverseComplementBases(complete);
        }
        String left = mExonBases.substring(0, start);
        String middle = mExonBases.substring(start, end);
        String right = mExonBases.substring(end);
        Preconditions.checkArgument((left + middle + right).equals(mExonBases));
        return mBasesOfFirstCodonInPreviousExon + left + middle + middle + right + mBasesOfLastCodonInFollowingExon;
    }

    public String baseSequenceWithInsertionApplied(final int position, final String basesToInsert, final boolean positiveStrand)
    {
        Preconditions.checkArgument(position >= 0);
        Preconditions.checkArgument(position <= mExonBases.length());
        if(!positiveStrand)
        {
            int s = mExonBases.length() - position;
            String l = mExonBases.substring(0, s);
            String r = mExonBases.substring(s);
            Preconditions.checkArgument((l + r).equals(mExonBases));
            String complete = mBasesOfFirstCodonInPreviousExon + l + basesToInsert + r + mBasesOfLastCodonInFollowingExon;
            return Nucleotides.reverseComplementBases(complete);
        }
        String left = mExonBases.substring(0, position);
        String right = mExonBases.substring(position);
        Preconditions.checkArgument((left + right).equals(mExonBases));
        return mBasesOfFirstCodonInPreviousExon + left + basesToInsert + right + mBasesOfLastCodonInFollowingExon;
    }

    public Pair<String, String> forwardStrandBaseAndLeftNeighbour(final int position, final boolean isPositiveStrand)
    {
        Preconditions.checkArgument(position >= 0);
        Preconditions.checkArgument(position < mExonBases.length());
        if(!isPositiveStrand)
        {
            int p = mExonBases.length() - position - 1;
            String right = mExonBases.substring(p, p + 1);
            String left = p == 0 ? last(mIntronicPrefix) : mExonBases.substring(p - 1, p);
            return Pair.of(left, right);
        }
        String left = position == 0 ? last(mIntronicPrefix) : mExonBases.substring(position - 1, position);
        String right = mExonBases.substring(position, position + 1);
        return Pair.of(left, right);
    }

    public String baseSequenceWithSingleBaseRemoved(final int position, final boolean isPositiveStrand)
    {
        Preconditions.checkArgument(position >= 0);
        Preconditions.checkArgument(position < mExonBases.length());
        if(!isPositiveStrand)
        {
            String right = mExonBases.substring(mExonBases.length() - position);
            String left = mExonBases.substring(0, mExonBases.length() - position - 1);
            String leftPadding = last(mIntronicPrefix);
            String complete = mBasesOfFirstCodonInPreviousExon + leftPadding + left + right + mBasesOfLastCodonInFollowingExon;
            return Nucleotides.reverseComplementBases(complete);
        }
        String left = mExonBases.substring(0, position);
        String right = position == mExonBases.length() - 1 ? "" : mExonBases.substring(position + 1);
        // To return a whole number of codons we pad from the suffix.
        // This will allow a sequence of amino acids to be generated from
        // the result, even though it might correspond to lost junction.
        String rightPadding = mIntronicSuffix.substring(0, 1);
        return mBasesOfFirstCodonInPreviousExon + left + right + rightPadding + mBasesOfLastCodonInFollowingExon;
    }

    public String baseSequenceWithBasesReplacedAtStrandLocation(final int location, final String current, final String replacementBases)
    {
        // Prefix|-------l      r----|Suffix
        int leftEnd = location - mExonStart;
        int rightEnd = leftEnd + current.length();
        String existing = mExonBases.substring(leftEnd, rightEnd);
        Preconditions.checkArgument(existing.equals(current));

        String left = mExonBases.substring(0, leftEnd);
        String right = mExonBases.substring(rightEnd);
        return mBasesOfFirstCodonInPreviousExon + left + replacementBases + right + mBasesOfLastCodonInFollowingExon;
    }

    public String basesBetween(int start, int end)
    {
        Preconditions.checkArgument(start >= 0);
        Preconditions.checkArgument(start <= end);
        return mExonBases.substring(start, end + 1);
    }

    public String baseAt(int start, boolean positiveStrand)
    {
        Preconditions.checkArgument(start >= 0);
        if(!positiveStrand)
        {
            int s = mExonBases.length() - start - 1;
            return mExonBases.substring(s, s + 1);
        }
        return mExonBases.substring(start, start + 1);
    }

    public String baseImmediatelyBefore(int position)
    {
        Preconditions.checkArgument(position >= 0);
        Preconditions.checkArgument(position <= inExonLength());
        if(position == 0)
        {
            return last(mIntronicPrefix);
        }
        return mExonBases.substring(position - 1, position);
    }

    public CodonWithinExons getCodon(int codonIndex, boolean isPositiveStrand)
    {
        Preconditions.checkArgument(codonIndex >= 0);
        if(!isPositiveStrand)
        {
            int start = codonLocationInExonBody(codonIndex, false);
            int stop = start + 3;
            if(start < 0)
            {
                String forwardBasesInExon = mExonBases.substring(0, stop);
                String fragmentInExon = Nucleotides.reverseComplementBases(forwardBasesInExon);
                BaseSequence bodySequence = new BaseSequence(mExonStart, fragmentInExon, false);
                String suffix = Nucleotides.reverseComplementBases(mBasesOfFirstCodonInPreviousExon);
                return CodonWithinExons.factory("", bodySequence, suffix);
            }
            if(stop >= mExonBases.length())
            {
                String forwardBasesInExon = mExonBases.substring(start);
                String fragmentInExon = Nucleotides.reverseComplementBases(forwardBasesInExon);
                BaseSequence bodySequence = new BaseSequence(start + mExonStart, fragmentInExon, false);
                String prefix = Nucleotides.reverseComplementBases(mBasesOfLastCodonInFollowingExon);
                return CodonWithinExons.factory(prefix, bodySequence, "");
            }
            String forwardBases = mExonBases.substring(start, start + 3);
            String fragmentInExon = Nucleotides.reverseComplementBases(forwardBases);
            BaseSequence bodySequence = new BaseSequence(start + mExonStart, fragmentInExon, false);
            return CodonWithinExons.factory("", bodySequence, "");
        }
        int start = codonLocationInExonBody(codonIndex, true);
        int stop = start + 3;
        if(start < 0)
        {
            String fragmentInExon = mExonBases.substring(0, stop);
            BaseSequence bodySequence = new BaseSequence(mExonStart, fragmentInExon, true);
            return CodonWithinExons.factory(mBasesOfFirstCodonInPreviousExon, bodySequence, "");
        }
        if(stop >= mExonBases.length())
        {
            String fragmentInExon = mExonBases.substring(start);
            BaseSequence bodySequence = new BaseSequence(start + mExonStart, fragmentInExon, true);
            return CodonWithinExons.factory("", bodySequence, mBasesOfLastCodonInFollowingExon);
        }
        String fragmentInExon = mExonBases.substring(start, stop);
        BaseSequence bodySequence = new BaseSequence(start + mExonStart, fragmentInExon, true);
        return CodonWithinExons.factory("", bodySequence, "");
    }

    public SplitCodonSequence getSplitSequenceForCodons(int startCodon, int count, final boolean isPositiveStrand)
    {
        Preconditions.checkArgument(startCodon >= 0);
        if(startCodon == 0)
        {
            Preconditions.checkArgument(!mBasesOfFirstCodonInPreviousExon.isBlank());
        }
        Preconditions.checkArgument(count > 0);
        if(!isPositiveStrand)
        {
            int stop = codonLocationInExonBody(startCodon, false) + 3; //end of deletion is the codon end.
            int start = stop - (3 * count);
            if(start < 0)
            {
                String forwardBasesInExon = mExonBases.substring(0, stop);
                String reversed = Nucleotides.reverseComplementBases(forwardBasesInExon);
                String suffix = Nucleotides.reverseComplementBases(mBasesOfFirstCodonInPreviousExon);
                return new SplitCodonSequence(reversed, suffix, mExonStart);
            }
            if(stop >= mExonBases.length())
            {
                String forwardBasesInExon = mExonBases.substring(start);
                String reversed = Nucleotides.reverseComplementBases(forwardBasesInExon);
                String prefix = Nucleotides.reverseComplementBases(mBasesOfLastCodonInFollowingExon);
                return new SplitCodonSequence(prefix, reversed, mExonStart + start);
            }
            String forwardBases = mExonBases.substring(start, stop);
            String reversed = Nucleotides.reverseComplementBases(forwardBases);
            return new SplitCodonSequence(reversed, "", mExonStart + start);
        }
        int start = codonLocationInExonBody(startCodon, true);
        int stop = start + count * 3;

        if(start < 0)
        {
            return new SplitCodonSequence(mBasesOfFirstCodonInPreviousExon, mExonBases.substring(0, stop), mExonStart);
        }
        if(stop >= mExonBases.length())
        {
            return new SplitCodonSequence(mExonBases.substring(start), mBasesOfLastCodonInFollowingExon, start + mExonStart);
        }
        String fragmentInExon = mExonBases.substring(start, stop);
        return new SplitCodonSequence(fragmentInExon, "", start + mExonStart);
    }

    public int codonLocationInExonBody(int codonIndex, final boolean isPositiveStrand)
    {
        Preconditions.checkArgument(codonIndex >= 0);

        if(!isPositiveStrand)
        {
            int b = mExonBases.length();
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
        return mBasesOfFirstCodonInPreviousExon.length();
    }

    @VisibleForTesting
    public int numberOfBasesInFollowingExon()
    {
        return mBasesOfLastCodonInFollowingExon.length();
    }

    public int toStrandCoordinates(int positionInExon, boolean positiveStrand)
    {
        Preconditions.checkArgument(positionInExon >= 0);
        Preconditions.checkArgument(positionInExon <= inExonLength());
        if(positiveStrand)
        {
            return positionInExon + mExonStart;
        }
        return mExonStart + mExonBases.length() - positionInExon;
    }

    public int fromStrandCoordinates(int positionInStrand, boolean positiveStrand)
    {
        Preconditions.checkArgument(positionInStrand >= 0);
        if(positiveStrand)
        {
            return positionInStrand - mExonStart;
        }
        return mExonStart + mExonBases.length() - positionInStrand;
    }

    public int inExonLength()
    {
        return mExonBases.length();
    }

    public int totalLength()
    {
        return mBasesOfFirstCodonInPreviousExon.length() + mBasesOfLastCodonInFollowingExon.length() + inExonLength();
    }

    @Override
    public boolean equals(final Object o)
    {
        if(o == null || getClass() != o.getClass())
        {
            return false;
        }
        final PaddedExon that = (PaddedExon) o;
        return Objects.equals(mBasesOfFirstCodonInPreviousExon, that.mBasesOfFirstCodonInPreviousExon)
                && Objects.equals(mBasesOfLastCodonInFollowingExon, that.mBasesOfLastCodonInFollowingExon)
                && Objects.equals(mExonBases, that.mExonBases);
    }

    @Override
    public int hashCode()
    {
        return Objects.hash(mBasesOfFirstCodonInPreviousExon, mBasesOfLastCodonInFollowingExon, mExonBases);
    }

    @Override
    public String toString()
    {
        return "PaddedExon{" +
                "prefixFromPreviousExon='" + mBasesOfFirstCodonInPreviousExon + '\'' +
                ", suffixFromNextExon='" + mBasesOfLastCodonInFollowingExon + '\'' +
                ", exonBases='" + mExonBases + '\'' +
                '}';
    }
}
