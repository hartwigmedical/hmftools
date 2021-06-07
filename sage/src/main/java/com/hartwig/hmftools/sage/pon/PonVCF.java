package com.hartwig.hmftools.sage.pon;

import java.util.List;

import org.jetbrains.annotations.NotNull;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

public class PonVCF implements AutoCloseable
{

    public static final String PON_COUNT = "PON_COUNT";
    public static final String PON_TOTAL = "PON_TOTAL";
    public static final String PON_MAX = "PON_MAX";

    private final VariantContextWriter writer;

    PonVCF(final String output, int sampleSize)
    {
        writer = new VariantContextWriterBuilder().setOutputFile(output)
                .modifyOption(Options.INDEX_ON_THE_FLY, false)
                .modifyOption(Options.USE_ASYNC_IO, false)
                .modifyOption(Options.DO_NOT_WRITE_GENOTYPES, true)
                .build();

        final VCFHeader header = new VCFHeader();
        header.addMetaDataLine(new VCFInfoHeaderLine(PON_COUNT, 1, VCFHeaderLineType.Integer, "how many samples had the variant"));
        header.addMetaDataLine(new VCFInfoHeaderLine(PON_TOTAL, 1, VCFHeaderLineType.Integer, "total depth"));
        header.addMetaDataLine(new VCFInfoHeaderLine(PON_MAX, 1, VCFHeaderLineType.Integer, "max depth"));
        header.addMetaDataLine(new VCFHeaderLine("PonInputSampleCount", String.valueOf(sampleSize)));
        writer.writeHeader(header);
    }

    public void write(@NotNull final List<VariantContext> contexts)
    {
        contexts.forEach(writer::add);
    }

    @Override
    public void close()
    {
        writer.close();
    }

}
