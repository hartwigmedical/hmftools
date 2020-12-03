package com.hartwig.hmftools.common.drivercatalog.panel;

import java.util.List;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

class GermlineBlacklistVCF {

    public static void main(String[] args) {
        new GermlineBlacklistVCF().process("/Users/jon/hmf/resources/KnownBlacklist.germline.hg19.vcf.gz",  GermlineBlacklist.grch37Blacklist());
        new GermlineBlacklistVCF().process("/Users/jon/hmf/resources/KnownBlacklist.germline.hg38.vcf.gz",  GermlineBlacklist.grch38Blacklist());
    }


    static final String BLACKLIST_FLAG = "BLACKLIST";

    public void process(final String outputFile, List<VariantContext> variants) {

        // WRITE
        VariantContextWriter writer = new VariantContextWriterBuilder().setOutputFile(outputFile)
                .modifyOption(Options.INDEX_ON_THE_FLY, true)
                .modifyOption(Options.USE_ASYNC_IO, false)
                .build();

        VCFHeader writerheader = new VCFHeader();
        writerheader.addMetaDataLine(new VCFInfoHeaderLine(BLACKLIST_FLAG, 1, VCFHeaderLineType.Flag, "Blacklisted location"));
        writer.writeHeader(writerheader);

        variants.forEach(writer::add);
        writer.close();

    }

}
