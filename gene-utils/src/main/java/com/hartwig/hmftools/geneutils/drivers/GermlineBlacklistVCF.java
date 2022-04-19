package com.hartwig.hmftools.geneutils.drivers;

import java.util.List;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

public final class GermlineBlacklistVCF
{
    private static final String BLACKLIST_FLAG = "BLACKLIST";

    public static void write(final String outputFile, final List<VariantContext> variants)
    {
        VariantContextWriter writer = new VariantContextWriterBuilder().setOutputFile(outputFile)
                .modifyOption(Options.INDEX_ON_THE_FLY, true)
                .modifyOption(Options.USE_ASYNC_IO, false)
                .build();

        VCFHeader header = new VCFHeader();
        header.addMetaDataLine(new VCFInfoHeaderLine(BLACKLIST_FLAG, 1, VCFHeaderLineType.Flag, "Blacklisted location"));
        writer.writeHeader(header);

        variants.forEach(writer::add);
        writer.close();
    }
}
