package com.hartwig.hmftools.sullivan;

import org.jetbrains.annotations.NotNull;

class RecreatedFastqHeaderParser implements FastqHeaderParser {

    public String apply(@NotNull String fastqHeader) {
        return fastqHeader.substring(0, fastqHeader.length() - 2);
    }
}
