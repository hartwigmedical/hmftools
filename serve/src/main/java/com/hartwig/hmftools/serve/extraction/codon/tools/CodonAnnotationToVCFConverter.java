package com.hartwig.hmftools.serve.extraction.codon.tools;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.serve.Knowledgebase;
import com.hartwig.hmftools.serve.extraction.codon.KnownCodon;
import com.hartwig.hmftools.serve.extraction.codon.KnownCodonFile;
import com.hartwig.hmftools.serve.extraction.util.KeyFormatter;

import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.core.config.Configurator;
import org.apache.logging.log4j.util.Strings;
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

public class CodonAnnotationToVCFConverter {
    private static final Logger LOGGER = LogManager.getLogger(CodonAnnotationToVCFConverter.class);
    private static final boolean LOG_DEBUG = true;

    public static void main(String[] args) throws IOException {
        LOGGER.info("Running SERVE codon VCF converter");

        if (LOG_DEBUG) {
            Configurator.setRootLevel(Level.DEBUG);
        }

        String knownCodonsTsv = System.getProperty("user.home") + "/hmf/tmp/serve/KnownCodons.SERVE.37.tsv";
        String outputFile = System.getProperty("user.home") + "/hmf/tmp/codon.vcf.gz";

        List<KnownCodon> codons = KnownCodonFile.read(knownCodonsTsv);

        LOGGER.info("The size of the file is {}", codons.size());

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

        String randomAltBase = Strings.EMPTY;
        for (KnownCodon codon : codons) {
            Long start = codon.annotation().start();
            Long end = codon.annotation().end();
            Long bewtweenNumber = end > start ? start + 1 : start - 1;
            List<Long> genomicPositions = Lists.newArrayList(start, bewtweenNumber, end);
            String chromosome = codon.annotation().chromosome();

            for (Long genomicPosition : genomicPositions) {
                String extractRefBaseOfPosition = extractRefBaseOfGenomicPosition(chromosome, genomicPosition);

                if (extractRefBaseOfPosition.equals("A")) {
                    randomAltBase = "T";
                } else if (extractRefBaseOfPosition.equals("C")) {
                    randomAltBase = "A";
                } else if (extractRefBaseOfPosition.equals("T")) {
                    randomAltBase = "G";
                } else if (extractRefBaseOfPosition.equals("G")) {
                    randomAltBase = "C";
                }

                extactAnnotationVariantCodonIndex(extractRefBaseOfPosition,
                        randomAltBase,
                        chromosome,
                        genomicPosition,
                        codon.sources(),
                        codon.annotation().gene(),
                        codon.annotation().proteinAnnotation(),
                        codon.annotation().transcript(),
                        writer);
            }
        }

        writer.close();
        LOGGER.info("All codons are written to VCF file!");
        LOGGER.info("Done!");
    }

    private static String extractRefBaseOfGenomicPosition(@Nullable String chromosome, Long genomicPosition) throws IOException {
        IndexedFastaSequenceFile fastaSequenceFile =
                new IndexedFastaSequenceFile(new File("/Users/liekeschoenmaker/hmf/refgenome/Homo_sapiens.GRCh37.GATK.illumina.fasta"));
        return fastaSequenceFile.getSubsequenceAt(chromosome, genomicPosition, genomicPosition).getBaseString();

    }

    private static void extactAnnotationVariantCodonIndex(@NotNull String extractRefBaseOfPosition, @NotNull String randomAltBase,
            @Nullable String chromosome, Long position,
            @NotNull Set<Knowledgebase> knowledgebases, @NotNull String gene, @NotNull String proteinAnnotation, @NotNull String transcript,
             @NotNull VariantContextWriter writer)  {

        generateVcfFileOfGenomicPosition(extractRefBaseOfPosition,
                randomAltBase,
                chromosome,
                position,
                knowledgebases,
                gene,
                proteinAnnotation,
                transcript, writer);

    }

    private static void generateVcfFileOfGenomicPosition(@NotNull String extractRefBaseOfPosition, @NotNull String randomAltBase,
            @Nullable String chromosome, Long position, @NotNull Set<Knowledgebase> knowledgebases, @NotNull String gene,
            @NotNull String proteinAnnotation, @NotNull String transcript, @NotNull VariantContextWriter writer) {


        List<Allele> hotspotAlleles =
                Lists.newArrayList(Allele.create(extractRefBaseOfPosition, true), Allele.create(randomAltBase, false));

        VariantContext variantContext = new VariantContextBuilder().noGenotypes()
                .source("SERVE")
                .chr(chromosome)
                .start(position)
                .alleles(hotspotAlleles)
                .computeEndFromAlleles(hotspotAlleles, new Long(position).intValue())
                .attribute("source", Knowledgebase.commaSeparatedSourceString(knowledgebases))
                .attribute("input", KeyFormatter.toProteinKey(gene, transcript, proteinAnnotation))
                .make();

        LOGGER.debug(" Writing variant to VCF file'{}'", variantContext);
        writer.add(variantContext);
    }

}
