package com.hartwig.hmftools.retentionchecker;

import org.jetbrains.annotations.NotNull;

class RecreatedFastqHeaderNormalizer implements FastqHeaderNormalizer {

    public String apply(@NotNull String fastqHeader) {
        return fastqHeader.substring(0, fastqHeader.length() - 2);
    }
}
