package com.hartwig.hmftools.pavereverse.protein;

import static com.hartwig.hmftools.pavereverse.util.Checks.isNucleotideSequence;

import java.util.Objects;

import com.google.common.base.Preconditions;

import org.apache.commons.lang3.StringUtils;

public class DeletionInsertionChange
{
    private final String mRef;
    private final String mAlt;
    private final String mInserted;
    private final String mDeleted;
    private final int DeletionStart;

    static int lengthOfCommonPrefix(String s, String t)
    {
        int length = Math.min(s.length(), t.length());
        for(int i = 0; i < length; i++)
        {
            if(s.charAt(i) != t.charAt(i))
            {
                return i;
            }
        }
        return length;
    }

    static int lengthOfCommonSuffix(String s, String t)
    {
        return lengthOfCommonPrefix(StringUtils.reverse(s), StringUtils.reverse(t));
    }

    DeletionInsertionChange(String ref, String alt)
    {
        Preconditions.checkArgument(isNucleotideSequence(ref));
        Preconditions.checkArgument(isNucleotideSequence(alt));
        mRef = ref;
        mAlt = alt;
        int commonSuffixLength = lengthOfCommonSuffix(ref, alt);
        String refWithoutCommonSuffix = ref.substring(0, ref.length() - commonSuffixLength);
        String altWithoutCommonSuffix = alt.substring(0, alt.length() - commonSuffixLength);
        int commonPrefixLength = lengthOfCommonPrefix(refWithoutCommonSuffix, altWithoutCommonSuffix);
        mDeleted = refWithoutCommonSuffix.substring(commonPrefixLength);
        mInserted = altWithoutCommonSuffix.substring(commonPrefixLength);
        DeletionStart = commonPrefixLength;
    }

    String deleted()
    {
        return mDeleted;
    }

    String inserted()
    {
        return mInserted;
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
        return Objects.equals(mRef, that.mRef) && Objects.equals(mAlt, that.mAlt);
    }

    @Override
    public int hashCode()
    {
        return Objects.hash(mRef, mAlt);
    }

    @Override
    public String toString()
    {
        return "DeletionInsertionChange{" +
                "Ref='" + mRef + '\'' +
                ", Alt='" + mAlt + '\'' +
                ", Inserted='" + mInserted + '\'' +
                ", Deleted='" + mDeleted + '\'' +
                ", DeletionStart=" + DeletionStart +
                '}';
    }
}
