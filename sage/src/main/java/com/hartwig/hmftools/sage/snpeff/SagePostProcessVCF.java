package com.hartwig.hmftools.sage.snpeff;

import org.jetbrains.annotations.NotNull;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFHeader;

public class SagePostProcessVCF implements AutoCloseable {

    private final VariantContextWriter writer;

    public SagePostProcessVCF(@NotNull final String outputVCF) {
        writer = new VariantContextWriterBuilder().setOutputFile(outputVCF).build();
    }

    public void writeHeader(VCFHeader header) {
        writer.writeHeader(header);
    }

    public void writeVariant(@NotNull final VariantContext context) {
        writer.add(context);
    }

    @Override
    public void close() {
        writer.close();
    }
}
