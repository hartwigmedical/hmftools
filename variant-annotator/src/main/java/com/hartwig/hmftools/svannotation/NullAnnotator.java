package com.hartwig.hmftools.svannotation;

import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.variant.structural.StructuralVariant;
import com.hartwig.hmftools.svannotation.annotations.StructuralVariantAnnotation;

public class NullAnnotator implements VariantAnnotator {

    public static VariantAnnotator make() {
        return new NullAnnotator();
    }

    @Override
    public List<StructuralVariantAnnotation> annotateVariants(final List<StructuralVariant> variants) {
        return variants.stream().map(StructuralVariantAnnotation::new).collect(Collectors.toList());
    }
}
