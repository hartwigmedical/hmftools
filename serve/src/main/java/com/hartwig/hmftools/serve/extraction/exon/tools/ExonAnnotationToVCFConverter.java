package com.hartwig.hmftools.serve.extraction.exon.tools;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Set;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.serve.Knowledgebase;
import com.hartwig.hmftools.serve.extraction.exon.KnownExon;
import com.hartwig.hmftools.serve.extraction.exon.KnownExonFile;
import com.hartwig.hmftools.serve.extraction.util.GenerateAltBase;
import com.hartwig.hmftools.serve.extraction.util.KeyFormatter;
import com.hartwig.hmftools.serve.extraction.util.VCFWriterFactory;

import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.core.config.Configurator;
import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;

public class ExonAnnotationToVCFConverter {

    private static final Logger LOGGER = LogManager.getLogger(ExonAnnotationToVCFConverter.class);
    private static final boolean LOG_DEBUG = true;

    public static void main(String[] args) throws IOException {
        LOGGER.info("Running SERVE exon VCF converter");

        if (LOG_DEBUG) {
            Configurator.setRootLevel(Level.DEBUG);
        }

        String knownExonsTsv = System.getProperty("user.home") + "/hmf/tmp/serve/KnownExons.SERVE.37.tsv";
        String outputFile = System.getProperty("user.home") + "/hmf/tmp/exons.vcf.gz";
        GenerateAltBase altBaseGenerator = new GenerateAltBase(RefGenomeVersion.V37,
                new IndexedFastaSequenceFile(new File(
                        System.getProperty("user.home") + "/hmf/refgenomes/grch37/Homo_sapiens.GRCh37.GATK.illumina.fasta")));

        List<KnownExon> exons = KnownExonFile.read(knownExonsTsv);
        LOGGER.info("The number of known exons in the known exon file is '{}'", exons.size());

        VariantContextWriter writer = VCFWriterFactory.generateVCFWriterWithInputAndSources(outputFile);

        for (KnownExon exon : exons) {
            String chromosome = exon.annotation().chromosome();
            long start = exon.annotation().start() + 5; // remove the first 5 slice position before exon
            long end = exon.annotation().end() - 5; // remove the last 5 splice position after exon
            long middle = start + Math.round((end - start) / 2D); // take middle position of exon
            List<Long> genomicPositions = Lists.newArrayList(start, middle, end);

            String gene = exon.annotation().gene();
            for (long position : genomicPositions) {
                String refBaseOfPosition = altBaseGenerator.extractRefBaseAtGenomicPosition(chromosome, position);
                String randomAltBase = altBaseGenerator.createAltForRefBase(chromosome, position);

                writeVariantToVCF(writer,
                        chromosome,
                        position,
                        refBaseOfPosition,
                        randomAltBase,
                        exon.sources(),
                        gene,
                        exon.annotation().exonIndex(),
                        exon.annotation().transcript());
            }
        }

        writer.close();

        LOGGER.info("All known exons are converted and written to VCF file!");
    }

    private static void writeVariantToVCF(@NotNull VariantContextWriter writer, @NotNull String chromosome, long position,
            @NotNull String ref, @NotNull String alt, @NotNull Set<Knowledgebase> knowledgebases, @NotNull String gene, int exonIndex,
            @NotNull String transcript) {
        List<Allele> alleles = Lists.newArrayList(Allele.create(ref, true), Allele.create(alt, false));

        VariantContext variantContext = new VariantContextBuilder().noGenotypes()
                .source("ExonChecker")
                .chr(chromosome)
                .start(position)
                .alleles(alleles)
                .computeEndFromAlleles(alleles, new Long(position).intValue())
                .attribute("source", Knowledgebase.toCommaSeparatedSourceString(knowledgebases))
                .attribute("input", KeyFormatter.toExonKey(gene, transcript, exonIndex))
                .make();

        LOGGER.debug(" Writing variant to VCF file'{}'", variantContext);
        writer.add(variantContext);
    }
}
