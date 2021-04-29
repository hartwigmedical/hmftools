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

    public static final String INPUT_FIELD = "input";
    public static final String SOURCES_FIELD = "sources";

    private VCFWriterFactory() {
    }

    @NotNull
    public static VariantContextWriter openVCFWriterWithInputAndSources(@NotNull String outputFile, @NotNull String sources) {
        VariantContextWriter writer = new VariantContextWriterBuilder().setOutputFile(outputFile)
                .setOutputFileType(VariantContextWriterBuilder.OutputType.BLOCK_COMPRESSED_VCF)
                .modifyOption(Options.INDEX_ON_THE_FLY, false)
                .build();

        VCFHeader header = new VCFHeader(Sets.newHashSet(), Lists.newArrayList());
        header.addMetaDataLine(new VCFInfoHeaderLine(INPUT_FIELD,
                VCFHeaderLineCount.UNBOUNDED,
                VCFHeaderLineType.String,
                "Input based on which SERVE generated this record"));
        header.addMetaDataLine(new VCFInfoHeaderLine(SOURCES_FIELD,
                VCFHeaderLineCount.UNBOUNDED,
                VCFHeaderLineType.String,
                "SERVE sources that contained this record [" + sources + "]"));

        writer.writeHeader(header);
        return writer;
    }
}
