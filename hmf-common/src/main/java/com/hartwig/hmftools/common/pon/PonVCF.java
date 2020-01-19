package com.hartwig.hmftools.common.pon;

import java.util.List;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

public class PonVCF implements AutoCloseable {

    public static final String PON_COUNT = "PON_COUNT";

    private final VariantContextWriter writer;

    PonVCF(String output) {
        writer = new VariantContextWriterBuilder().setOutputFile(output)
                .modifyOption(Options.INDEX_ON_THE_FLY, false)
                .modifyOption(Options.USE_ASYNC_IO, false)
                .modifyOption(Options.DO_NOT_WRITE_GENOTYPES, true)
                .build();

    }

    public void write(int samples, List<VariantContext> contexts) {

        final VCFHeader header = new VCFHeader();
        header.addMetaDataLine(new VCFInfoHeaderLine(PON_COUNT, 1, VCFHeaderLineType.Integer, "how many samples had the variant"));
        header.addMetaDataLine(new VCFHeaderLine("PonInputSampleCount", String.valueOf(samples)));
        writer.writeHeader(header);

        contexts.forEach(writer::add);
    }

    @Override
    public void close() {
        writer.close();
    }

}
