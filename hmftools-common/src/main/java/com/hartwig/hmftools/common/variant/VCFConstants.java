package com.hartwig.hmftools.common.variant;

import java.util.Arrays;
import java.util.List;

public final class VCFConstants {

    public static final String MUTECT = "mutect";
    public static final String VARSCAN = "varscan";
    public static final String STRELKA = "strelka";
    public static final String FREEBAYES = "freebayes";

    public static final List<String> ALL_CALLERS = Arrays.asList(MUTECT, VARSCAN, STRELKA, FREEBAYES);

    private VCFConstants() {
    }
}
