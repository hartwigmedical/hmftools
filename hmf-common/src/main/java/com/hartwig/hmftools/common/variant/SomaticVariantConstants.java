package com.hartwig.hmftools.common.variant;

import java.util.Arrays;
import java.util.List;

final class SomaticVariantConstants {

    static final String MUTECT = "mutect";
    static final String VARSCAN = "varscan";
    static final String STRELKA = "strelka";
    static final String FREEBAYES = "freebayes";

    static final List<String> ALL_CALLERS = Arrays.asList(MUTECT, VARSCAN, STRELKA, FREEBAYES);

    private SomaticVariantConstants() {
    }
}
