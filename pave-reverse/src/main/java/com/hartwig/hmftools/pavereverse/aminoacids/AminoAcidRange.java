package com.hartwig.hmftools.pavereverse.aminoacids;

import java.util.Objects;

import com.google.common.base.Preconditions;
import com.hartwig.hmftools.common.gene.TranscriptAminoAcids;
import com.hartwig.hmftools.pavereverse.TranscriptFilter;

/**
 * Specifies the start and stop of a range of amino acids in a transcript.
 */
public class AminoAcidRange implements TranscriptFilter
{
    private final AminoAcidSpecification mFirst;
    private final AminoAcidSpecification mLast;

    public AminoAcidRange(AminoAcidSpecification first, AminoAcidSpecification last)
    {
        Preconditions.checkArgument(first.Position <= last.Position, "End position must not be before start position");
        mFirst = first;
        mLast = last;
    }

    public int startPosition()
    {
        return mFirst.Position;
    }

    public AminoAcid aminoAcidAtStart()
    {
        return mFirst.value();
    }

    public int length()
    {
        return mLast.Position - mFirst.Position + 1;
    }

    @Override
    public boolean applies(TranscriptAminoAcids aminoAcids)
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
