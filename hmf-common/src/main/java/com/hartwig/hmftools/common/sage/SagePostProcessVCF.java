package com.hartwig.hmftools.common.sage;

import java.util.List;
import java.util.function.Consumer;

import com.hartwig.hmftools.common.genome.region.CanonicalTranscript;
import com.hartwig.hmftools.common.variant.enrich.SnpEffEnrichment;

import org.jetbrains.annotations.NotNull;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFHeader;

public class SagePostProcessVCF implements AutoCloseable, Consumer<VariantContext> {
    private final VariantContextWriter writer;
    private final SnpEffEnrichment snpEffEnrichment;
    final List<CanonicalTranscript> transcripts;

    public SagePostProcessVCF(@NotNull final String outputVCF, final List<CanonicalTranscript> transcripts) {
        writer = new VariantContextWriterBuilder().setOutputFile(outputVCF).build();
        this.transcripts = transcripts;
        this.snpEffEnrichment = new SnpEffEnrichment(transcripts, writer::add);
    }

    public void writeHeader(@NotNull final VCFHeader header) {
        writer.writeHeader(snpEffEnrichment.enrichHeader(header));
    }

    @Override
    public void close() {
        writer.close();
    }

    @Override
    public void accept(final VariantContext variantContext) {
        snpEffEnrichment.accept(variantContext);
    }
}
