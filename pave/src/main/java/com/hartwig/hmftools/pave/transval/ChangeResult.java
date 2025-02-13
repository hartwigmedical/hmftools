package com.hartwig.hmftools.pave.transval;

import java.util.Objects;

import com.google.common.base.Preconditions;

import org.jetbrains.annotations.NotNull;

class ChangeResult
{
    @NotNull
    public final AminoAcidSequence mAminoAcids;

    @NotNull
    public final String mBases;

    public ChangeResult(@NotNull final AminoAcidSequence mAminoAcids, @NotNull final String mBases)
    {
        this.mAminoAcids = mAminoAcids;
        this.mBases = mBases;
    }

    @Override
    public boolean equals(final Object o)
    {
        if(o == null || getClass() != o.getClass())
        {
            return false;
        }
        final ChangeResult that = (ChangeResult) o;
        return Objects.equals(mAminoAcids, that.mAminoAcids) && Objects.equals(mBases, that.mBases);
    }

    @Override
    public int hashCode()
    {
        return Objects.hash(mAminoAcids, mBases);
    }

    @Override
    public String toString()
    {
        return "ChangeResult{" +
                "mAminoAcids='" + mAminoAcids + '\'' +
                ", mBases='" + mBases + '\'' +
                '}';
    }
}
