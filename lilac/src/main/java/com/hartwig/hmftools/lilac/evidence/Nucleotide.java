package com.hartwig.hmftools.lilac.evidence;

import java.util.Collection;
import java.util.List;

public record Nucleotide(int locus, byte qual, String bases)
{
    public static List<Byte> qualities(final Collection<Nucleotide> nucleotides)
    {
        return nucleotides.stream().map(Nucleotide::qual).toList();
    }

    public Builder builder()
    {
        return new Builder();
    }

    public class Builder
    {
        private int mLocus = locus;
        private byte mQual = qual;
        private String mBases = bases;

        public Builder setLocus(int locus)
        {
            mLocus = locus;
            return this;
        }

        public Builder setQual(byte qual)
        {
            mQual = qual;
            return this;
        }

        public Builder setBases(final String bases)
        {
            mBases = bases;
            return this;
        }

        public Nucleotide build()
        {
            return new Nucleotide(mLocus, mQual, mBases);
        }
    }
}