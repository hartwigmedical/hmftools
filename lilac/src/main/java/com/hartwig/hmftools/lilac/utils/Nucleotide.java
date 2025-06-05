package com.hartwig.hmftools.lilac.utils;

import java.util.List;

import org.jetbrains.annotations.NotNull;

public record Nucleotide(int locus, byte qual, @NotNull String bases)
{
    // TODO: get around having this?
    public static List<Integer> loci(@NotNull final List<Nucleotide> nucleotides)
    {
        return nucleotides.stream().map(Nucleotide::locus).toList();
    }

    // TODO: get around having this?
    public static List<Byte> qualities(@NotNull final List<Nucleotide> nucleotides)
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
        private @NotNull String mBases = bases;

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

        public Builder setBases(@NotNull final String bases)
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