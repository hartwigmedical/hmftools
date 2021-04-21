package com.hartwig.hmftools.common.genome.refgenome;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.SamReader;

public final class RefGenomeFunctions {

    private static final String CHR_PREFIX = "chr";

    private RefGenomeFunctions() {
    }

    public static boolean samReaderUsesChrInContigs(@NotNull SamReader samReader) {
        return samReader.getFileHeader()
                .getSequenceDictionary()
                .getSequences()
                .stream()
                .anyMatch(x -> x.getSequenceName().contains(CHR_PREFIX));
    }

    @NotNull
    public static String stripChromosome(@NotNull final String chromosome) {
        if (chromosome.startsWith(CHR_PREFIX)) {
            return chromosome.substring(CHR_PREFIX.length());
        }

        return chromosome;
    }

    @NotNull
    public static String enforceChromosome(@NotNull final String chromosome) {
        if (!chromosome.startsWith(CHR_PREFIX)) {
            return CHR_PREFIX + chromosome;
        }

        return chromosome;
    }
}
