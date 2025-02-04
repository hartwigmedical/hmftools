package com.hartwig.hmftools.pave.transval;

import java.util.Objects;

import com.google.common.base.Preconditions;
import com.hartwig.hmftools.common.gene.TranscriptAminoAcids;

import org.jetbrains.annotations.NotNull;

/**
 * Specifies the start and stop of a range of amino acids in a transcript.
 */
class AminoAcidRange implements TranscriptFilter
{
    @NotNull
    private final AminoAcidSpecification first;
    @NotNull
    private final AminoAcidSpecification last;

    AminoAcidRange(@NotNull final AminoAcidSpecification first, @NotNull final AminoAcidSpecification last)
    {
        Preconditions.checkArgument(first.position <= last.position, "End position must not be before start position");
        this.first = first;
        this.last = last;
    }

    public int startPosition()
    {
        return first.position;
    }

    public int length()
    {
        return last.position - first.position + 1;
    }

    @Override
    public boolean applies(final TranscriptAminoAcids aminoAcids)
    {
        return first.applies(aminoAcids) && last.applies(aminoAcids);
    }

    @Override
    public boolean equals(final Object o)
    {
        if(o == null || getClass() != o.getClass())
        {
            return false;
        }
        final AminoAcidRange that = (AminoAcidRange) o;
        return Objects.equals(first, that.first) && Objects.equals(last, that.last);
    }

    @Override
    public int hashCode()
    {
        return Objects.hash(first, last);
    }

    @Override
    public String toString()
    {
        return first + "_" + last;
    }
}
