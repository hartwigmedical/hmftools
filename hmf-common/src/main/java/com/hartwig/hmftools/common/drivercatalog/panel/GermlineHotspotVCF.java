package com.hartwig.hmftools.common.drivercatalog.panel;

import java.io.File;
import java.util.Set;

import com.google.common.collect.Sets;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;

class GermlineHotspotVCF {

    private final String contigPrefix;
    private final Set<String> germlineGenes;

    public GermlineHotspotVCF(final String contigPrefix, final Set<String> genes) {
        this.contigPrefix = contigPrefix;
        this.germlineGenes = genes;
    }

    public void process(final String inputFile, final String outputFile) {

        VariantContextWriter writer = new VariantContextWriterBuilder().setOutputFile(outputFile)
                .modifyOption(Options.INDEX_ON_THE_FLY, true)
                .modifyOption(Options.USE_ASYNC_IO, false)
                .build();

        VCFFileReader reader = new VCFFileReader(new File(inputFile), false);

        VCFHeader readerHeader = reader.getFileHeader();
        VCFHeader writerheader = new VCFHeader();
        writerheader.addMetaDataLine(readerHeader.getInfoHeaderLine("GENEINFO"));
        writerheader.addMetaDataLine(readerHeader.getInfoHeaderLine("CLNSIG"));
        writer.writeHeader(writerheader);

        for (VariantContext context : reader) {
            final Set<String> variantGenes = Sets.newHashSet();
            final String geneinfo = context.getAttributeAsString("GENEINFO", "UNKNOWN:0000");
            for (String singleGeneInfo : geneinfo.split("\\|")) {
                variantGenes.add(singleGeneInfo.split(":")[0]);
            }

            // Intersection of variant and approved genes
            variantGenes.retainAll(germlineGenes);

            String clinsig = context.getAttributeAsString("CLNSIG", "unknown");
            if (!variantGenes.isEmpty() && isPathogenicOrLikelyPathogenic(clinsig) && context.getAlleles().size() == 2) {
                VariantContextBuilder builder = new VariantContextBuilder("clinvar",
                        contigPrefix + context.getContig(),
                        context.getStart(),
                        context.getEnd(),
                        context.getAlleles()).attribute("GENEINFO", geneinfo).attribute("CLNSIG", clinsig);

                writer.add(builder.make());
            }
        }

        reader.close();
        writer.close();

    }

    private boolean isPathogenicOrLikelyPathogenic(String clnsig) {
        return clnsig.startsWith("Pathogenic") || clnsig.equals("Likely_pathogenic");
    }

}
