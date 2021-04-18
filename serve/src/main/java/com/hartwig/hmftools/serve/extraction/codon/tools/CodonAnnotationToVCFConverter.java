package com.hartwig.hmftools.serve.extraction.codon.tools;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Set;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.serve.Knowledgebase;
import com.hartwig.hmftools.serve.extraction.codon.KnownCodon;
import com.hartwig.hmftools.serve.extraction.codon.KnownCodonFile;
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

public class CodonAnnotationToVCFConverter {

    private static final Logger LOGGER = LogManager.getLogger(CodonAnnotationToVCFConverter.class);
    private static final boolean LOG_DEBUG = true;

    public static void main(String[] args) throws IOException {
        LOGGER.info("Running SERVE codon VCF converter");

        if (LOG_DEBUG) {
            Configurator.setRootLevel(Level.DEBUG);
        }

        String knownCodonsTsv = System.getProperty("user.home") + "/hmf/tmp/serve/KnownCodons.SERVE.37.tsv";
        String outputFile = System.getProperty("user.home") + "/hmf/tmp/codons.vcf.gz";
        GenerateAltBase altBaseGenerator = new GenerateAltBase(RefGenomeVersion.V37,
                new IndexedFastaSequenceFile(new File(
                        System.getProperty("user.home") + "/hmf/refgenomes/grch37/Homo_sapiens.GRCh37.GATK.illumina.fasta")));

        List<KnownCodon> codons = KnownCodonFile.read(knownCodonsTsv);

        LOGGER.info("The number of codons in known codon file is '{}'", codons.size());

        VariantContextWriter writer = VCFWriterFactory.generateVCFWriterWithInputAndSources(outputFile);

        for (KnownCodon codon : codons) {
            long start = codon.annotation().start();
            long end = codon.annotation().end();
            long middle = start + 1;
            List<Long> genomicPositions = Lists.newArrayList(start, middle, end);

            String chromosome = codon.annotation().chromosome();
            for (long position : genomicPositions) {
                String refBaseOfPosition = altBaseGenerator.extractRefBaseAtGenomicPosition(chromosome, position);
                String randomAltBase = altBaseGenerator.createAltForRefBase(chromosome, position);

                writeVariantToVCF(writer,
                        chromosome,
                        position,
                        refBaseOfPosition,
                        randomAltBase,
                        codon.sources(),
                        codon.annotation().gene(),
                        codon.annotation().transcript(),
                        codon.annotation().codonIndex());
            }
        }

        writer.close();

        LOGGER.info("All known codons are converted and written to VCF file!");
    }

    private static void writeVariantToVCF(@NotNull VariantContextWriter writer, @NotNull String chromosome, long position,
            @NotNull String ref, @NotNull String alt, @NotNull Set<Knowledgebase> knowledgebases, @NotNull String gene,
            @NotNull String transcript, int codonIndex) {
        List<Allele> alleles = Lists.newArrayList(Allele.create(ref, true), Allele.create(alt, false));

        VariantContext variantContext = new VariantContextBuilder().noGenotypes()
                .source("CodonChecker")
                .chr(chromosome)
                .start(position)
                .alleles(alleles)
                .computeEndFromAlleles(alleles, new Long(position).intValue())
                .attribute("source", Knowledgebase.toCommaSeparatedSourceString(knowledgebases))
                .attribute("input", KeyFormatter.toCodonKey(gene, transcript, codonIndex))
                .make();

        LOGGER.debug(" Writing variant to VCF file'{}'", variantContext);
        writer.add(variantContext);
    }
}
