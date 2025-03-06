package com.hartwig.hmftools.pavereverse;

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
    private final AminoAcidSpecification mFirst;
    @NotNull
    private final AminoAcidSpecification mLast;

    AminoAcidRange(@NotNull final AminoAcidSpecification first, @NotNull final AminoAcidSpecification last)
    {
        Preconditions.checkArgument(first.mPosition <= last.mPosition, "End position must not be before start position");
        this.mFirst = first;
        this.mLast = last;
    }

    public int startPosition()
    {
        return mFirst.mPosition;
    }

    public AminoAcid aminoAcidAtStart()
    {
        return mFirst.value();
    }

    public int length()
    {
        return mLast.mPosition - mFirst.mPosition + 1;
    }

    @Override
    public boolean applies(final TranscriptAminoAcids aminoAcids)
    {
        return mFirst.applies(aminoAcids) && mLast.applies(aminoAcids);
    }

    @Override
    public boolean equals(final Object o)
    {
        if(o == null || getClass() != o.getClass())
        {
            return false;
        }
        final AminoAcidRange that = (AminoAcidRange) o;
        return Objects.equals(mFirst, that.mFirst) && Objects.equals(mLast, that.mLast);
    }

    @Override
    public int hashCode()
    {
        return Objects.hash(mFirst, mLast);
    }

    @Override
    public String toString()
    {
        return mFirst + "_" + mLast;
    }
}
