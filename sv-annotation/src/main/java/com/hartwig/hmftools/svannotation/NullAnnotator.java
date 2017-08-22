package com.hartwig.hmftools.svannotation;

import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.variant.structural.StructuralVariant;

public class NullAnnotator implements StructuralVariantAnnotator {

    public static StructuralVariantAnnotator make() {
        return new NullAnnotator();
    }

    @Override
    public List<StructuralVariantAnnotation> annotateVariants(final List<StructuralVariant> variants) {
        return variants.stream().map(StructuralVariantAnnotation::new).collect(Collectors.toList());
    }
}
