package com.hartwig.hmftools.sullivan;

import org.jetbrains.annotations.NotNull;

class OriginalFastqHeaderParser implements FastqHeaderParser {

    public String apply(@NotNull String fastqHeader) {
        return fastqHeader.substring(0, fastqHeader.indexOf(" "));
    }
}
