package com.hartwig.hmftools.lilac.evidence;

import static com.hartwig.hmftools.common.redux.BaseQualAdjustment.isLowBaseQual;

import java.util.Collection;
import java.util.List;

public record Nucleotide(int locus, boolean isLowQual, byte qual, String bases)
{
    public static final byte MISSING_BASE_QUAL = -1;

    public static Nucleotide create(int locus, byte qual, final String bases)
    {
        return new Nucleotide(locus, isLowBaseQual(qual), qual, bases);
    }

    public static Nucleotide create(int locus, boolean isLowQual, final String bases)
    {
        return new Nucleotide(locus, isLowQual, MISSING_BASE_QUAL, bases);
    }

    public static Nucleotide createHighQual(int locus, final String bases)
    {
        return new Nucleotide(locus, false, MISSING_BASE_QUAL, bases);
    }

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
        private boolean mIsLowQual = isLowQual;
        private byte mQual = qual;
        private String mBases = bases;

        public Builder setLocus(int locus)
        {
            mLocus = locus;
            return this;
        }

        public Builder setIsLowQual(boolean isLowQual)
        {
            mIsLowQual = isLowQual;
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
            return new Nucleotide(mLocus, mIsLowQual, mQual, mBases);
        }
    }
}
