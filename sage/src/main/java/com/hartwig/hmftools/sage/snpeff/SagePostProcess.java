package com.hartwig.hmftools.sage.snpeff;

import java.util.List;
import java.util.function.Consumer;

import com.hartwig.hmftools.common.genome.region.CanonicalTranscript;
import com.hartwig.hmftools.common.sage.SagePostProcessVCF;
import com.hartwig.hmftools.common.variant.snpeff.SnpEffSummary;
import com.hartwig.hmftools.common.variant.snpeff.SnpEffSummaryFactory;
import com.hartwig.hmftools.common.variant.snpeff.SnpEffSummarySerialiser;

import org.jetbrains.annotations.NotNull;

import htsjdk.variant.variantcontext.VariantContext;

public class SagePostProcess implements AutoCloseable, Consumer<VariantContext> {

    private final Consumer<VariantContext> consumer;
    private final SnpEffSummaryFactory snpEffSummaryFactory;

    public SagePostProcess(@NotNull final List<CanonicalTranscript> transcripts, @NotNull final Consumer<VariantContext> consumer) {
        this.consumer = consumer;
        this.snpEffSummaryFactory = new SnpEffSummaryFactory(transcripts);
    }

    @Override
    public void close() {
    }

    @Override
    public void accept(@NotNull final VariantContext context) {
        final SnpEffSummary snpEffSummary = snpEffSummaryFactory.fromAnnotations(context);
        if (!snpEffSummary.worstGene().isEmpty()) {
            context.getCommonInfo().putAttribute(SagePostProcessVCF.SNPEFF_WORST, SnpEffSummarySerialiser.worstDetails(snpEffSummary));
        }
        if (!snpEffSummary.canonicalGene().isEmpty()) {
            context.getCommonInfo()
                    .putAttribute(SagePostProcessVCF.SNPEFF_CANONICAL, SnpEffSummarySerialiser.canonicalDetails(snpEffSummary));
        }
        consumer.accept(context);
    }
}
