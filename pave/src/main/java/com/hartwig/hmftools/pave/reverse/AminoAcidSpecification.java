package com.hartwig.hmftools.pave.reverse;

import static com.hartwig.hmftools.pave.reverse.Checks.matchPattern;

import java.util.List;
import java.util.Objects;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import com.google.common.base.Preconditions;
import com.hartwig.hmftools.common.gene.TranscriptAminoAcids;

import org.jetbrains.annotations.NotNull;

/**
 * A 1-based position in a gene transcript and an amino acid that is in that position,
 * or null if none is specified.
 */
class AminoAcidSpecification implements TranscriptFilter
{
    final int mPosition;
    private final AminoAcid mAminoAcid;

    @NotNull
    public static AminoAcidSpecification parse(String aaPos)
    {
        Pattern pattern = Pattern.compile("([A-Z](?:[a-z][a-z])?+)(\\d+)");
        final Matcher matcher = matchPattern(pattern, aaPos);

        return new AminoAcidSpecification(Integer.parseInt(matcher.group(2)), matcher.group(1));
    }

    public AminoAcidSpecification(final int position, final String aminoAcid)
    {
        this(position, aminoAcid != null ? new AminoAcid(aminoAcid) : null);
    }

    public AminoAcidSpecification(final int position, final AminoAcid aminoAcid)
    {
        Preconditions.checkArgument(position > 0);
        this.mPosition = position;
        this.mAminoAcid = aminoAcid;
    }

    public AminoAcid value()
    {
        return mAminoAcid;
    }

    @Override
    public boolean applies(final TranscriptAminoAcids aminoAcids)
    {
        if (mPosition > aminoAcids.AminoAcids.length())
        {
            return false;
        }
        if (mAminoAcid == null)
        {
            return true;
        }
        return aminoAcids.AminoAcids.substring(mPosition - 1, mPosition).equals(mAminoAcid.mSymbol);
    }

    @Override
    public boolean equals(final Object o)
    {
        if(o == null || getClass() != o.getClass())
        {
            return false;
        }
        final AminoAcidSpecification that = (AminoAcidSpecification) o;
        return mPosition == that.mPosition && Objects.equals(mAminoAcid, that.mAminoAcid);
    }

    @Override
    public int hashCode()
    {
        return Objects.hash(mPosition, mAminoAcid);
    }

    @Override
    public String toString()
    {
        return "[" + mPosition + "," + symbol() + ']';
    }

    @NotNull
    public String symbol()
    {
        return mAminoAcid == null ? "?" : mAminoAcid.mSymbol;
    }
}
