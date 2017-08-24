package com.hartwig.hmftools.svannotation;

import java.util.List;

import com.hartwig.hmftools.common.region.GenomeRegion;
import com.hartwig.hmftools.common.variant.structural.StructuralVariant;

public interface StructuralVariantAnnotator {
    List<StructuralVariantAnnotation> annotateVariants(final List<StructuralVariant> variants);
    StructuralVariantAnnotation annotateRegion(final GenomeRegion region);
}
