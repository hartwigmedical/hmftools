package com.hartwig.hmftools.svannotation;

import java.util.List;

import com.hartwig.hmftools.common.region.GenomeRegion;
import com.hartwig.hmftools.common.variant.structural.StructuralVariant;

public interface VariantAnnotator {
    List<VariantAnnotation> annotateVariants(final List<StructuralVariant> variants);
    VariantAnnotation annotateRegion(final GenomeRegion region);
}
