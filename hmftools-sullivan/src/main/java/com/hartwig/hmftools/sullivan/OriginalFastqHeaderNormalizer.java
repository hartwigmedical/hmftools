package com.hartwig.hmftools.sullivan;

import org.jetbrains.annotations.NotNull;

class OriginalFastqHeaderNormalizer implements FastqHeaderNormalizer {

    FastqHeaderNormalizer recreatedNormalizer = new RecreatedFastqHeaderNormalizer();

    public String apply(@NotNull String fastqHeader) {
        // KODU: This is a bit ugly, but if we recognize an original fastq file is actually recreated,
        // we try to recreated normalizer...
        int spaceIndex = fastqHeader.indexOf(" ");
        return spaceIndex > 0 ? fastqHeader.substring(0, spaceIndex) : recreatedNormalizer.apply(fastqHeader);
    }
}
