package com.hartwig.hmftools.pavereverse.aminoacids;

import static com.hartwig.hmftools.pavereverse.util.Checks.matchPattern;

import java.util.Objects;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import com.google.common.base.Preconditions;
import com.hartwig.hmftools.common.gene.TranscriptAminoAcids;
import com.hartwig.hmftools.pavereverse.TranscriptFilter;

/**
 * A 1-based position in a gene transcript and an amino acid that is in that position,
 * or null if none is specified.
 */
public class AminoAcidSpecification implements TranscriptFilter
{
    public final int Position;
    private final AminoAcid Acid;

    public static AminoAcidSpecification parse(String aaPos)
    {
        Pattern pattern = Pattern.compile("([A-Z](?:[a-z][a-z])?+)(\\d+)");
        Matcher matcher = matchPattern(pattern, aaPos);

        return new AminoAcidSpecification(Integer.parseInt(matcher.group(2)), matcher.group(1));
    }

    public AminoAcidSpecification(int position, String aminoAcid)
    {
        this(position, aminoAcid != null ? new AminoAcid(aminoAcid) : null);
    }

    public AminoAcidSpecification(int position, AminoAcid aminoAcid)
    {
        Preconditions.checkArgument(position > 0);
        Position = position;
        Acid = aminoAcid;
    }

    public AminoAcid value()
    {
        return Acid;
    }

    @Override
    public boolean applies(TranscriptAminoAcids aminoAcids)
    {
        if(Position > aminoAcids.AminoAcids.length())
        {
            return false;
        }
        if(Acid == null)
        {
            return true;
        }
        return aminoAcids.AminoAcids.substring(Position - 1, Position).equals(Acid.Symbol);
    }

    @Override
    public boolean equals(final Object o)
    {
        if(o == null || getClass() != o.getClass())
        {
            return false;
        }
        final AminoAcidSpecification that = (AminoAcidSpecification) o;
        return Position == that.Position && Objects.equals(Acid, that.Acid);
    }

    @Override
    public int hashCode()
    {
        return Objects.hash(Position, Acid);
    }

    @Override
    public String toString()
    {
        return "[" + Position + "," + symbol() + ']';
    }

    public String symbol()
    {
        return Acid == null ? "?" : Acid.Symbol;
    }
}
