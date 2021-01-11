package com.hartwig.hmftools.serve.extraction.util;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;

public class GenerateAltBase {

    @NotNull
    private final IndexedFastaSequenceFile fastaSequenceFile;

    public GenerateAltBase(@NotNull final IndexedFastaSequenceFile fastaSequenceFile) {
        this.fastaSequenceFile = fastaSequenceFile;
    }

    @NotNull
    public String createAltForRefBase(@NotNull String chromosome, long position) {
        String refBaseAtPosition = extractRefBaseAtGenomicPosition(chromosome, position);

        switch (refBaseAtPosition) {
            case "A":
                return "T";
            case "C":
                return "A";
            case "T":
                return "G";
            case "G":
                return "C";
            default:
                throw new IllegalArgumentException("Cannot generate alt for ref base '" + refBaseAtPosition + "'");
        }
    }

    @NotNull
    public String extractRefBaseAtGenomicPosition(@NotNull String chromosome, long position) {
        return fastaSequenceFile.getSubsequenceAt(chromosome, position, position).getBaseString();
    }
}
