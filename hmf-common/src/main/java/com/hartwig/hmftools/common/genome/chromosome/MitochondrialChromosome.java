package com.hartwig.hmftools.common.genome.chromosome;

import org.jetbrains.annotations.NotNull;

public enum MitochondrialChromosome implements Chromosome
{
    MT {
        @Override
        public String contig()
        {
            return "MT"; // todo test
        }
    };

    public static final int MT_LENGTH = 16569;

    @Override
    public boolean isAutosome()
    {
        return false;
    }

    @Override
    public boolean isAllosome()
    {
        return false;
    }

    @NotNull
    public static MitochondrialChromosome fromString(@NotNull final String contig)
    {
        if(!contains(contig))
        {
            throw new IllegalArgumentException("Invalid mitochondrial contig " + contig);
        }

        return MT;
    }

    public static boolean contains(@NotNull final String contig)
    {
        return contig.equals("chrM") || contig.equals("MT");
    }
}
