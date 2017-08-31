package com.hartwig.hmftools.svannotation;

import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.region.GenomeRegion;
import com.hartwig.hmftools.common.variant.structural.StructuralVariant;

public class NullAnnotator implements VariantAnnotator {

    public static VariantAnnotator make() {
        return new NullAnnotator();
    }

    @Override
    public List<VariantAnnotation> annotateVariants(final List<StructuralVariant> variants) {
        return variants.stream().map(VariantAnnotation::new).collect(Collectors.toList());
    }

    @Override
    public VariantAnnotation annotateRegion(final GenomeRegion region) {
        return new VariantAnnotation(null);
    }
}
