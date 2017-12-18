package com.hartwig.hmftools.bachelor;

import java.util.function.Predicate;

public class BachelorProgram {
    final Predicate<VariantModel> vcfProcessor;

    public BachelorProgram(final Predicate<VariantModel> vcfProcessor) {
        this.vcfProcessor = vcfProcessor;
    }
}
