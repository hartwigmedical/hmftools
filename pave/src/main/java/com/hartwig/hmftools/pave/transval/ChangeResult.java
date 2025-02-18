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

    public final int mLocation;

    @NotNull
    public final String mRefBases;

    @NotNull
    public final String altBases;

    public ChangeResult(@NotNull final AminoAcidSequence mAminoAcids, @NotNull final String mBases, final int mLocation,
            @NotNull final String mRefBases, @NotNull final String altBases)
    {
        this.mAminoAcids = mAminoAcids;
        this.mBases = mBases;
        this.mLocation = mLocation;
        this.mRefBases = mRefBases;
        this.altBases = altBases;
    }

    public TransvalHotspot toHotspot(String chromosome)
    {
        return new TransvalHotspot(mRefBases, altBases, chromosome, mLocation);
    }

    @Override
    public boolean equals(final Object o)
    {
        if(o == null || getClass() != o.getClass())
        {
            return false;
        }
        final ChangeResult that = (ChangeResult) o;
        return mLocation == that.mLocation && Objects.equals(mAminoAcids, that.mAminoAcids)
                && Objects.equals(mBases, that.mBases) && Objects.equals(mRefBases, that.mRefBases)
                && Objects.equals(altBases, that.altBases);
    }

    @Override
    public int hashCode()
    {
        return Objects.hash(mAminoAcids, mBases, mLocation, mRefBases, altBases);
    }

    @Override
    public String toString()
    {
        return "ChangeResult{" +
                "mAminoAcids=" + mAminoAcids +
                ", mBases='" + mBases + '\'' +
                ", mLocation=" + mLocation +
                ", mRefBases='" + mRefBases + '\'' +
                ", altBases='" + altBases + '\'' +
                '}';
    }
}
