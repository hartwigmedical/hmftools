package com.hartwig.hmftools.retentionchecker;

import org.jetbrains.annotations.NotNull;

class OriginalFastqHeaderNormalizer implements FastqHeaderNormalizer {

    @NotNull
    private final FastqHeaderNormalizer recreatedNormalizer = new RecreatedFastqHeaderNormalizer();

    @NotNull
    public String apply(@NotNull final String fastqHeader) {
        // KODU: This is a bit ugly, but if we recognize an original fastq file is actually recreated,
        // we try to recreated normalizer...
        final int spaceIndex = fastqHeader.indexOf(" ");
        return spaceIndex > 0 ? fastqHeader.substring(0, spaceIndex) : recreatedNormalizer.apply(fastqHeader);
    }
}
