package com.hartwig.hmftools.common.sage;

import org.jetbrains.annotations.NotNull;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

public class SagePostProcessVCF implements AutoCloseable {

    public static String HMF_CANONICAL_GENE = "HCG";
    public static String HMF_CANONICAL_EFFECT = "HCE";
    public static String HMF_CANONICAL_IMPACT = "HCI";

    private final VariantContextWriter writer;

    public SagePostProcessVCF(@NotNull final String outputVCF) {
        writer = new VariantContextWriterBuilder().setOutputFile(outputVCF).build();
    }

    public void writeHeader(@NotNull final VCFHeader header) {
        header.addMetaDataLine(new VCFInfoHeaderLine(HMF_CANONICAL_GENE, 1, VCFHeaderLineType.String, "HMF canonical gene"));
        header.addMetaDataLine(new VCFInfoHeaderLine(HMF_CANONICAL_EFFECT, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "HMF canonical effect"));
        header.addMetaDataLine(new VCFInfoHeaderLine(HMF_CANONICAL_IMPACT, 1, VCFHeaderLineType.String, "HMF canonical impact"));
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
