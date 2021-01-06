package com.hartwig.hmftools.serve.extraction.exon.tools;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Random;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.serve.Knowledgebase;
import com.hartwig.hmftools.serve.extraction.exon.KnownExon;
import com.hartwig.hmftools.serve.extraction.exon.KnownExonFile;
import com.hartwig.hmftools.serve.extraction.hotspot.ProteinKeyFormatter;

import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.core.config.Configurator;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

public class ExonAnnotationToVCFConverter {
    private static final Logger LOGGER = LogManager.getLogger(ExonAnnotationToVCFConverter.class);
    private static final boolean LOG_DEBUG = true;

    public static void main(String[] args) throws IOException {
        LOGGER.info("Running SERVE exon checker");

        if (LOG_DEBUG) {
            Configurator.setRootLevel(Level.DEBUG);
        }

        String knownExonsTsv = System.getProperty("user.home") + "/hmf/tmp/serve/KnownExons.SERVE.37.tsv";
        String outputFile = System.getProperty("user.home") + "/hmf/tmp/exon.vcf.gz";

        List<KnownExon> exons = KnownExonFile.read(knownExonsTsv);
        LOGGER.info("The size of the file is {}", exons.size());

        VariantContextWriter writer = new VariantContextWriterBuilder().setOutputFile(outputFile)
                .setOutputFileType(VariantContextWriterBuilder.OutputType.BLOCK_COMPRESSED_VCF)
                .modifyOption(Options.INDEX_ON_THE_FLY, false)
                .build();

        VCFHeader header = new VCFHeader(Sets.newHashSet(), Lists.newArrayList());
        header.addMetaDataLine(new VCFInfoHeaderLine("input", VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "input"));
        header.addMetaDataLine(new VCFInfoHeaderLine("source",
                VCFHeaderLineCount.UNBOUNDED,
                VCFHeaderLineType.String,
                "sources"));

        writer.writeHeader(header);

        Random random = new Random();
        String randomAltBase = null;

        for (KnownExon exon : exons) {
            String chromosome = exon.annotation().chromosome();
            Long start = exon.annotation().start() + 5;
            Long end = exon.annotation().end() - 5;
            String gene = exon.annotation().gene();

            Long l = new Long(end - start + 1);
            Long bewtweenNumber = start + random.nextInt(l.intValue());
            List<Long> genomicPositions = Lists.newArrayList(start, bewtweenNumber, end);

            for (Long genomicPosition : genomicPositions) {
                String extractRefBaseOfPosition = extractRefBaseOfGenomicPosition(chromosome, genomicPosition);

                if (extractRefBaseOfPosition.equals("A")) {
                    randomAltBase = "T";
                } else if (extractRefBaseOfPosition.equals("C")) {
                    randomAltBase = "A";
                } else if (extractRefBaseOfPosition.equals("T")) {
                    randomAltBase = "G";
                } else if (extractRefBaseOfPosition.equals("G")) {
                    randomAltBase = "A";
                }

                extactAnnotationVariantExonIndex(extractRefBaseOfPosition,
                        randomAltBase,
                        chromosome,
                        genomicPosition,
                        exon.sources(),
                        gene,
                        exon.annotation().exonIndex(),
                        exon.annotation().exonEnsemblId(),
                        writer);

            }
        }
        writer.close();
        LOGGER.info("All exons are checked!");

        LOGGER.info("Done!");
    }

    private static String extractRefBaseOfGenomicPosition(@Nullable String chromosome, Long genomicPosition) throws IOException {
        IndexedFastaSequenceFile fastaSequenceFile =
                new IndexedFastaSequenceFile(new File("/Users/liekeschoenmaker/hmf/refgenome/Homo_sapiens.GRCh37.GATK.illumina.fasta"));
        return fastaSequenceFile.getSubsequenceAt(chromosome, genomicPosition, genomicPosition).getBaseString();
    }

    private static void extactAnnotationVariantExonIndex(@Nullable String extractRefBaseOfPosition, @Nullable String randomAltBase,
            @Nullable String chromosome, Long position,
            @NotNull Set<Knowledgebase> knowledgebases, @NotNull String gene, int exonIndex, @NotNull String transcript,
            @NotNull VariantContextWriter writer)  {

        generateVcfFileOfGenomicPosition(extractRefBaseOfPosition,
                randomAltBase,
                chromosome,
                position,
                knowledgebases,
                gene,
                exonIndex,
                transcript, writer);

    }

    private static void generateVcfFileOfGenomicPosition(@Nullable String extractRefBaseOfPosition, @Nullable String randomAltBase,
            @Nullable String chromosome, Long position, @NotNull Set<Knowledgebase> knowledgebases, @NotNull String gene,
            int exonIndex, @NotNull String transcript, @NotNull VariantContextWriter writer) {


        List<Allele> hotspotAlleles =
                Lists.newArrayList(Allele.create(extractRefBaseOfPosition, true), Allele.create(randomAltBase, false));

        VariantContext variantContext = new VariantContextBuilder().noGenotypes()
                .source("SERVE")
                .chr(chromosome)
                .start(position)
                .alleles(hotspotAlleles)
                .computeEndFromAlleles(hotspotAlleles, new Long(position).intValue())
                .attribute("source", Knowledgebase.commaSeparatedSourceString(knowledgebases))
                .attribute("input", ProteinKeyFormatter.toExonKey(gene, transcript, Integer.toString(exonIndex)))
                .make();

        LOGGER.debug(" Writing variant to VCF file'{}'", variantContext);
        writer.add(variantContext);
    }
}
