package com.hartwig.hmftools.pave.reverse;

import java.util.Objects;

import org.jetbrains.annotations.NotNull;

class ChangeResult
{
    @NotNull
    final AminoAcidSequence mAminoAcids;

    @NotNull
    final String mBases;

    final int mLocation;

    @NotNull
    final String mRefBases;

    @NotNull
    final String mAltBases;

    ChangeResult(@NotNull final AminoAcidSequence mAminoAcids, @NotNull final String mBases, final int mLocation,
            @NotNull final String mRefBases, @NotNull final String altBases)
    {
        this.mAminoAcids = mAminoAcids;
        this.mBases = mBases;
        this.mLocation = mLocation;
        this.mRefBases = mRefBases;
        this.mAltBases = altBases;
    }

    BaseSequenceChange asChange(String chromosome)
    {
        return new BaseSequenceChange(mRefBases, mAltBases, chromosome, mLocation);
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
                && Objects.equals(mAltBases, that.mAltBases);
    }

    @Override
    public int hashCode()
    {
        return Objects.hash(mAminoAcids, mBases, mLocation, mRefBases, mAltBases);
    }

    @Override
    public String toString()
    {
        return "ChangeResult{" +
                "mAminoAcids=" + mAminoAcids +
                ", mBases='" + mBases + '\'' +
                ", mLocation=" + mLocation +
                ", mRefBases='" + mRefBases + '\'' +
                ", altBases='" + mAltBases + '\'' +
                '}';
    }
}
