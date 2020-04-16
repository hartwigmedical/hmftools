package com.hartwig.hmftools.common.sage;

import org.jetbrains.annotations.NotNull;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

public class SagePostProcessVCF implements AutoCloseable {

    public static String SNPEFF_WORST = "SEW";
    public static String SNPEFF_CANONICAL = "SEC";

    private final VariantContextWriter writer;

    public SagePostProcessVCF(@NotNull final String outputVCF) {
        writer = new VariantContextWriterBuilder().setOutputFile(outputVCF).build();
    }

    public void writeHeader(@NotNull final VCFHeader header) {
        header.addMetaDataLine(new VCFInfoHeaderLine(SNPEFF_WORST,
                5,
                VCFHeaderLineType.String,
                "SnpEff worst transcript summary [Gene, Transcript, Effect, CodingEffect, GenesAffected]"));
        header.addMetaDataLine(new VCFInfoHeaderLine(SNPEFF_CANONICAL,
                6,
                VCFHeaderLineType.String,
                "SnpEff canonical transcript summary [Gene, Transcript, Effect, CodingEffect, HgvsCodingImpact, HgvsProteinImpact]"));
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
