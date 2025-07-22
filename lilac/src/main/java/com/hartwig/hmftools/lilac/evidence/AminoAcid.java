package com.hartwig.hmftools.lilac.evidence;

public record AminoAcid(int locus, String acid)
{
    public Builder builder()
    {
        return new Builder();
    }

    public class Builder
    {
        private int mLocus = locus;
        private String mAcid = acid;

        public Builder setLocus(int locus)
        {
            mLocus = locus;
            return this;
        }

        public Builder setAcid(final String acid)
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
