package com.hartwig.hmftools.pave.transval;

import java.util.Objects;

import com.google.common.base.Preconditions;
import com.hartwig.hmftools.common.gene.TranscriptAminoAcids;

import org.jetbrains.annotations.NotNull;

/**
 * A 1-based position in a gene transcript and an amino acid that is in that position,
 * or null if none is specified.
 */
class AminoAcidSpecification implements TranscriptFilter
{
    final int position;
    private final AminoAcid aminoAcid;

    public AminoAcidSpecification(final int position, final String aminoAcid)
    {
        this(position, aminoAcid != null ? new AminoAcid(aminoAcid) : null);
    }

    public AminoAcidSpecification(final int position, final AminoAcid aminoAcid)
    {
        Preconditions.checkArgument(position > 0);
        this.position = position;
        this.aminoAcid = aminoAcid;
    }

    @Override
    public boolean applies(final TranscriptAminoAcids aminoAcids)
    {
        if (position > aminoAcids.AminoAcids.length())
        {
            return false;
        }
        if (aminoAcid == null)
        {
            return true;
        }
        return aminoAcids.AminoAcids.substring(position - 1, position).equals(aminoAcid.symbol);
    }

    @Override
    public boolean equals(final Object o)
    {
        if(o == null || getClass() != o.getClass())
        {
            return false;
        }
        final AminoAcidSpecification that = (AminoAcidSpecification) o;
        return position == that.position && Objects.equals(aminoAcid, that.aminoAcid);
    }

    @Override
    public int hashCode()
    {
        return Objects.hash(position, aminoAcid);
    }

    @Override
    public String toString()
    {
        return "[" + position + "," + symbol() + ']';
    }

    @NotNull
    public String symbol()
    {
        return aminoAcid == null ? "?" : aminoAcid.symbol;
    }
}
