package com.hartwig.hmftools.lilac.utils;

import org.jetbrains.annotations.NotNull;

public record AminoAcid(int locus, @NotNull String acid)
{
    public Builder builder()
    {
        return new Builder();
    }

    public class Builder
    {
        private int mLocus = locus;
        private @NotNull String mAcid = acid;

        public Builder setLocus(int locus)
        {
            mLocus = locus;
            return this;
        }

        public Builder setAcid(@NotNull final String acid)
        {
            mAcid = acid;
            return this;
        }

        public AminoAcid build()
        {
            return new AminoAcid(mLocus, mAcid);
        }
    }
}
