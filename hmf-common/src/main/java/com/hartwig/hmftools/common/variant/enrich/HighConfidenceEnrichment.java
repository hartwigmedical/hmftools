package com.hartwig.hmftools.common.variant.enrich;

import java.util.Optional;
import java.util.function.Consumer;

import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.genome.position.GenomePositions;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.genome.region.GenomeRegionSelector;
import com.hartwig.hmftools.common.genome.region.GenomeRegionSelectorFactory;

import org.jetbrains.annotations.NotNull;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

public class HighConfidenceEnrichment implements VariantContextEnrichment {

    public static final String HIGH_CONFIDENCE_FLAG = "HC";
    private static final String HIGH_CONFIDENCE_FLAG_DESCRIPTION = "High confidence region";

    private final Consumer<VariantContext> consumer;
    private final GenomeRegionSelector<GenomeRegion> highConfidenceSelector;

    HighConfidenceEnrichment(@NotNull final Multimap<String, GenomeRegion> highConfidenceRegions,
            @NotNull final Consumer<VariantContext> consumer) {
        this.highConfidenceSelector = GenomeRegionSelectorFactory.create(highConfidenceRegions);
        this.consumer = consumer;
    }

    @Override
    public void flush() {
        // Empty
    }

    @NotNull
    @Override
    public VCFHeader enrichHeader(@NotNull final VCFHeader template) {
        template.addMetaDataLine(new VCFInfoHeaderLine(HIGH_CONFIDENCE_FLAG,
                0,
                VCFHeaderLineType.Flag,
                HIGH_CONFIDENCE_FLAG_DESCRIPTION));

        return template;
    }

    @Override
    public void accept(@NotNull final VariantContext context) {
        Optional<GenomeRegion> region = highConfidenceSelector.select(GenomePositions.create(context.getContig(), context.getStart()));
        if (region.isPresent()) {
            context.getCommonInfo().putAttribute(HIGH_CONFIDENCE_FLAG, true);
        }

        consumer.accept(context);
    }
}
