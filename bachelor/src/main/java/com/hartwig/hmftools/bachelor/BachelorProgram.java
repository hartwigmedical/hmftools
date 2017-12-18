package com.hartwig.hmftools.bachelor;

import java.util.function.Predicate;

import com.hartwig.hmftools.common.gene.GeneCopyNumber;

public class BachelorProgram {
    final Predicate<VariantModel> vcfProcessor;
    final Predicate<GeneCopyNumber> copyNumberProcessor;

    public BachelorProgram(final Predicate<VariantModel> vcfProcessor, final Predicate<GeneCopyNumber> copyNumberProcessor) {
        this.vcfProcessor = vcfProcessor;
        this.copyNumberProcessor = copyNumberProcessor;
    }
}
