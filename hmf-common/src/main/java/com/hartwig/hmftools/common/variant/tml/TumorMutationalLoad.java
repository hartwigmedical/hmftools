package com.hartwig.hmftools.common.variant.tml;

import java.util.function.Consumer;

import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.common.variant.msi.MicrosatelliteIndels;
import com.hartwig.hmftools.common.variant.snpeff.SnpEffSummary;
import com.hartwig.hmftools.common.variant.snpeff.SnpEffSummaryFactory;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.filter.PassingVariantFilter;
import htsjdk.variant.variantcontext.filter.VariantContextFilter;

public class TumorMutationalLoad implements Consumer<VariantContext> {

    private static final VariantContextFilter PASS = new PassingVariantFilter();

    private int load;
    private int burden;

    public int load() {
        return load;
    }

    public double burdenPerMb() {
        return burden / MicrosatelliteIndels.NUMBER_OF_MB_PER_GENOME;
    }

    @Override
    public void accept(final VariantContext context) {
        if (PASS.test(context)) {
            burden++;

            final SnpEffSummary snpEffSummary = SnpEffSummaryFactory.fromSnpEffEnrichment(context);
            if (snpEffSummary.worstCodingEffect().equals(CodingEffect.MISSENSE)) {
                load++;
            }
        }
    }
}
