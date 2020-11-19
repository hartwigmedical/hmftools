package com.hartwig.hmftools.common.drivercatalog.panel;

import java.io.File;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Sets;

import org.apache.commons.compress.utils.Lists;
import org.apache.logging.log4j.util.Strings;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.VariantContextComparator;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

class GermlineHotspotVCF {

    static final String WHITELIST_FLAG = "WHITELIST";
    private final Set<String> germlineGenes;

    public GermlineHotspotVCF(final Set<String> genes) {
        this.germlineGenes = genes;
    }

    public void process(final String inputFile, final String outputFile) {

        // READ
        final VCFFileReader reader = new VCFFileReader(new File(inputFile), false);
        final VCFHeader readerHeader = reader.getFileHeader();
        final String assembly = readerHeader.getMetaDataLine("reference").getValue();
        if (!assembly.equals("GRCh37") && !assembly.equals("GRCh38")) {
            throw new IllegalArgumentException();
        }

        final String contigPrefix = assembly.equals("GRCh38") ? "chr" : Strings.EMPTY;
        final List<VariantContext> variants = Lists.newArrayList();

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

                variants.add(builder.make());
            }
        }
        reader.close();

        // ADD WHITE LIST (OUT OF ORDER)
        variants.addAll(assembly.equals("GRCh37") ? GermlineWhitelist.grch37Whitelist() : GermlineWhitelist.grch38Whitelist());

        // Get sorted contigs
        final List<String> contigs = variants.stream().map(VariantContext::getContig).distinct().collect(Collectors.toList());


        // WRITE
        VariantContextWriter writer = new VariantContextWriterBuilder().setOutputFile(outputFile)
                .modifyOption(Options.INDEX_ON_THE_FLY, true)
                .modifyOption(Options.USE_ASYNC_IO, false)
                .build();

        VCFHeader writerheader = new VCFHeader();
        writerheader.addMetaDataLine(readerHeader.getInfoHeaderLine("GENEINFO"));
        writerheader.addMetaDataLine(readerHeader.getInfoHeaderLine("CLNSIG"));
        writerheader.addMetaDataLine(new VCFInfoHeaderLine(WHITELIST_FLAG, 1, VCFHeaderLineType.Flag, "Whitelisted hotspot"));
        writer.writeHeader(writerheader);

        variants.sort(new VariantContextComparator(contigs));
        variants.forEach(writer::add);
        writer.close();

    }

    private boolean isPathogenicOrLikelyPathogenic(String clnsig) {
        return clnsig.startsWith("Pathogenic") || clnsig.equals("Likely_pathogenic");
    }

}
