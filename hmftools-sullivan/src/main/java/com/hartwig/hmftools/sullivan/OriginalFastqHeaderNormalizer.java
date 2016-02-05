package com.hartwig.hmftools.sullivan;

import org.jetbrains.annotations.NotNull;

class OriginalFastqHeaderNormalizer implements FastqHeaderNormalizer {

    public String apply(@NotNull String fastqHeader) {
        return fastqHeader.substring(0, fastqHeader.indexOf(" "));
    }
}
