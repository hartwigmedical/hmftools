package com.hartwig.hmftools.bachelor;

import java.util.function.Predicate;

import com.hartwig.hmftools.common.gene.GeneCopyNumber;
import com.hartwig.hmftools.common.region.hmfslicer.HmfGenomeRegion;

public class BachelorProgram {
    final String name;
    final Predicate<VariantModel> vcfProcessor;
    final Predicate<GeneCopyNumber> copyNumberProcessor;
    final Predicate<HmfGenomeRegion> disruptionProcessor;

    BachelorProgram(final String name, final Predicate<VariantModel> vcfProcessor, final Predicate<GeneCopyNumber> copyNumberProcessor,
            final Predicate<HmfGenomeRegion> disruptionProcessor) {
        this.name = name;
        this.vcfProcessor = vcfProcessor;
        this.copyNumberProcessor = copyNumberProcessor;
        this.disruptionProcessor = disruptionProcessor;
    }
}
