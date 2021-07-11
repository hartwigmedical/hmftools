package com.hartwig.hmftools.common.variant.snpeff;

import java.util.List;
import java.util.Set;
import java.util.function.Consumer;

import com.hartwig.hmftools.common.genome.region.HmfTranscriptRegion;
import com.hartwig.hmftools.common.variant.enrich.VariantContextEnrichment;
import com.hartwig.hmftools.common.variant.snpeff.SnpEffSummary;
import com.hartwig.hmftools.common.variant.snpeff.SnpEffSummaryFactory;
import com.hartwig.hmftools.common.variant.snpeff.SnpEffSummarySerialiser;

import org.jetbrains.annotations.NotNull;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

public class SnpEffEnrichment implements VariantContextEnrichment
{

    public static final String SNPEFF_WORST = "SEW";
    public static final String SNPEFF_CANONICAL = "SEC";

    private final Consumer<VariantContext> consumer;
    private final SnpEffSummaryFactory snpEffSummaryFactory;

    public SnpEffEnrichment(@NotNull final Set<String> driverGenes, @NotNull final List<HmfTranscriptRegion> transcripts,
            @NotNull final Consumer<VariantContext> consumer) {
        this.consumer = consumer;
        this.snpEffSummaryFactory = new SnpEffSummaryFactory(driverGenes, transcripts);
    }

    @Override
    public void accept(@NotNull final VariantContext context) {
        final SnpEffSummary snpEffSummary = snpEffSummaryFactory.fromSnpEffAnnotations(context);
        if (!snpEffSummary.worstGene().isEmpty()) {
            context.getCommonInfo().putAttribute(SNPEFF_WORST, SnpEffSummarySerialiser.worstDetails(snpEffSummary), true);
        }
        if (!snpEffSummary.canonicalGene().isEmpty()) {
            context.getCommonInfo().putAttribute(SNPEFF_CANONICAL, SnpEffSummarySerialiser.canonicalDetails(snpEffSummary), true);
        }
        consumer.accept(context);
    }

    @Override
    public void flush() {
        // Do nothing
    }

    @NotNull
    @Override
    public VCFHeader enrichHeader(@NotNull final VCFHeader header) {
        header.addMetaDataLine(new VCFInfoHeaderLine(SNPEFF_WORST,
                5,
                VCFHeaderLineType.String,
                "SnpEff worst transcript summary [Gene, Transcript, Effect, CodingEffect, GenesAffected]"));
        header.addMetaDataLine(new VCFInfoHeaderLine(SNPEFF_CANONICAL,
                6,
                VCFHeaderLineType.String,
                "SnpEff canonical transcript summary [Gene, Transcript, Effect, CodingEffect, HgvsCodingImpact, HgvsProteinImpact]"));

        return header;
    }
}
