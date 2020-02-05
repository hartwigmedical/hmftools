package com.hartwig.hmftools.sage.snpeff;

import org.jetbrains.annotations.NotNull;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

public class SagePostProcessVCF implements AutoCloseable {

    public static String PHASED_INFRAME_INDEL = "PII";

    private final VariantContextWriter writer;

    public SagePostProcessVCF(@NotNull final String outputVCF) {
        writer = new VariantContextWriterBuilder().setOutputFile(outputVCF).build();
    }

    public void writeHeader(@NotNull final VCFHeader header) {
        header.addMetaDataLine(new VCFInfoHeaderLine(PHASED_INFRAME_INDEL, 1, VCFHeaderLineType.Flag, "Phased inframe indel"));
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
