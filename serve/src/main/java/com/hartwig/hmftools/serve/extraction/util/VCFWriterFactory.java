package com.hartwig.hmftools.serve.extraction.util;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;

import org.jetbrains.annotations.NotNull;

import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

public final class VCFWriterFactory {

    private VCFWriterFactory() {
    }

    @NotNull
    public static VariantContextWriter generateVCFWriterWithInputAndSources(@NotNull String outputFile) {
        VariantContextWriter writer = new VariantContextWriterBuilder().setOutputFile(outputFile)
                .setOutputFileType(VariantContextWriterBuilder.OutputType.BLOCK_COMPRESSED_VCF)
                .modifyOption(Options.INDEX_ON_THE_FLY, false)
                .build();

        VCFHeader header = new VCFHeader(Sets.newHashSet(), Lists.newArrayList());
        header.addMetaDataLine(new VCFInfoHeaderLine("input", VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "input"));
        header.addMetaDataLine(new VCFInfoHeaderLine("sources", VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "sources"));

        writer.writeHeader(header);
        return writer;
    }
}
