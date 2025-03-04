package com.hartwig.hmftools.pave.reverse;

import static com.hartwig.hmftools.pave.reverse.Checks.isNucleotideSequence;

import java.util.Objects;

import com.google.common.base.Preconditions;

import org.apache.commons.lang3.StringUtils;
import org.jetbrains.annotations.NotNull;

class DeletionInsertionChange
{
    @NotNull
    private final String Ref;

    @NotNull
    private final String Alt;

    @NotNull
    private final String Inserted;

    @NotNull
    private final String Deleted;

    private final int DeletionStart;

    static int lengthOfCommonPrefix(@NotNull String s, @NotNull String t)
    {
        int length = Math.min(s.length(), t.length());
        for (int i = 0; i < length; i++)
        {
            if (s.charAt(i) != t.charAt(i))
            {
                return i;
            }
        }
        return length;
    }

    static int lengthOfCommonSuffix(@NotNull String s, @NotNull String t)
    {
        return lengthOfCommonPrefix(StringUtils.reverse(s), StringUtils.reverse(t));
    }

    DeletionInsertionChange(@NotNull final String ref, @NotNull final String alt)
    {
        Preconditions.checkArgument(isNucleotideSequence(ref));
        Preconditions.checkArgument(isNucleotideSequence(alt));
        this.Ref = ref;
        this.Alt = alt;
        int commonSuffixLength = lengthOfCommonSuffix(ref, alt);
        String refWithoutCommonSuffix = ref.substring(0, ref.length() - commonSuffixLength);
        String altWithoutCommonSuffix = alt.substring(0, alt.length() - commonSuffixLength);
        int commonPrefixLength = lengthOfCommonPrefix(refWithoutCommonSuffix, altWithoutCommonSuffix);
        Deleted = refWithoutCommonSuffix.substring(commonPrefixLength);
        Inserted = altWithoutCommonSuffix.substring(commonPrefixLength);
        DeletionStart = commonPrefixLength;
    }

    @NotNull
    BaseSequenceChange toHotspot(@NotNull final ChangeLocation changeLocation)
    {
        int localPosition = changeLocation.Location + DeletionStart;
        return new BaseSequenceChange(Deleted, Inserted, changeLocation.mChromosome, localPosition);
    }

    @NotNull
    String deleted()
    {
        return Deleted;
    }

    @NotNull
    String inserted()
    {
        return Inserted;
    }

    int positionOfDeletion()
    {
        return DeletionStart;
    }

    @Override
    public boolean equals(final Object o)
    {
        if(o == null || getClass() != o.getClass())
        {
            return false;
        }
        final DeletionInsertionChange that = (DeletionInsertionChange) o;
        return Objects.equals(Ref, that.Ref) && Objects.equals(Alt, that.Alt);
    }

    @Override
    public int hashCode()
    {
        return Objects.hash(Ref, Alt);
    }

    @Override
    public String toString()
    {
        return "DeletionInsertionChange{" +
                "Ref='" + Ref + '\'' +
                ", Alt='" + Alt + '\'' +
                ", Inserted='" + Inserted + '\'' +
                ", Deleted='" + Deleted + '\'' +
                ", DeletionStart=" + DeletionStart +
                '}';
    }
}
