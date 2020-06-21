package com.hartwig.hmftools.common.variant.tml;

import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.msi.MicrosatelliteIndels;
import com.hartwig.hmftools.common.variant.snpeff.SnpEffAnnotation;
import com.hartwig.hmftools.common.variant.snpeff.SnpEffAnnotationFactory;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.filter.PassingVariantFilter;
import htsjdk.variant.variantcontext.filter.VariantContextFilter;

public class TumorMutationalLoad  {

    private static final VariantContextFilter PASS = new PassingVariantFilter();

    private int load;
    private int burden;

    public int load() {
        return load;
    }

    public double burdenPerMb() {
        return burden / MicrosatelliteIndels.NUMBER_OF_MB_PER_GENOME;
    }

    public void accept(final SomaticVariant variant) {
        if (!variant.isFiltered()) {
            burden++;

            if (variant.worstCodingEffect() == CodingEffect.MISSENSE) {
                load++;
            }
        }
    }

    public void accept(final VariantContext context) {
        //TODO: Don't use this version until after NEAR_INDEL_PON incorporated into VCF

        if (PASS.test(context)) {
            burden++;

            final List<SnpEffAnnotation> allAnnotations = SnpEffAnnotationFactory.fromContext(context);
            final List<SnpEffAnnotation> transcriptAnnotations =
                    allAnnotations.stream().filter(SnpEffAnnotation::isTranscriptFeature).collect(Collectors.toList());

            if (!transcriptAnnotations.isEmpty()) {
                final SnpEffAnnotation worstAnnotation = transcriptAnnotations.get(0);
                final CodingEffect codingEffect = CodingEffect.effect(worstAnnotation.gene(), worstAnnotation.consequences());
                if (codingEffect == CodingEffect.MISSENSE) {
                    load++;
                }
            }
        }
    }
}
