package com.hartwig.hmftools.retentionchecker;

import org.jetbrains.annotations.NotNull;

class RecreatedFastqHeaderNormalizer implements FastqHeaderNormalizer {

    @NotNull
    public String apply(@NotNull final String fastqHeader) {
        return fastqHeader.substring(0, fastqHeader.length() - 2);
    }
}
