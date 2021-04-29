package com.hartwig.hmftools.serve.extraction.codon.tools;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.serve.Knowledgebase;
import com.hartwig.hmftools.serve.ServeConfig;
import com.hartwig.hmftools.serve.ServeLocalConfigProvider;
import com.hartwig.hmftools.serve.extraction.KnownEvent;
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

        ServeConfig config = ServeLocalConfigProvider.create();

        String knownCodonsTsv = System.getProperty("user.home") + "/hmf/tmp/serve/KnownCodons.SERVE.37.tsv";
        String outputVcf = System.getProperty("user.home") + "/hmf/tmp/codons.vcf.gz";
        GenerateAltBase altBaseGenerator =
                new GenerateAltBase(RefGenomeVersion.V37, new IndexedFastaSequenceFile(new File(config.refGenome37FastaFile())));

        List<KnownCodon> codons = KnownCodonFile.read(knownCodonsTsv);

        LOGGER.info("The number of codons in known codon file is {}", codons.size());

        VariantContextWriter writer = VCFWriterFactory.generateVCFWriterWithInputAndSources(outputVcf, uniqueSourcesString(codons));

        for (KnownCodon codon : codons) {
            long start = codon.annotation().start();
            long end = codon.annotation().end();
            long middle = start + 1;
            List<Long> positions = Lists.newArrayList(start, middle, end);

            String chromosome = codon.annotation().chromosome();
            for (long position : positions) {
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

        LOGGER.info("All known codons are converted and written to '{}'", outputVcf);
    }

    private static void writeVariantToVCF(@NotNull VariantContextWriter writer, @NotNull String chromosome, long position,
            @NotNull String ref, @NotNull String alt, @NotNull Set<Knowledgebase> knowledgebases, @NotNull String gene,
            @NotNull String transcript, int codonIndex) {
        List<Allele> alleles = Lists.newArrayList(Allele.create(ref, true), Allele.create(alt, false));

        VariantContext variant = new VariantContextBuilder().noGenotypes()
                .source("CodonChecker")
                .chr(chromosome)
                .start(position)
                .alleles(alleles)
                .computeEndFromAlleles(alleles, new Long(position).intValue())
                .attribute(VCFWriterFactory.INPUT_FIELD, KeyFormatter.toCodonKey(gene, transcript, codonIndex))
                .attribute(VCFWriterFactory.SOURCES_FIELD, Knowledgebase.toCommaSeparatedSourceString(knowledgebases))
                .make();

        LOGGER.debug(" Writing '{}' to VCF file", variant);
        writer.add(variant);
    }

    @NotNull
    private static String uniqueSourcesString(@NotNull Iterable<? extends KnownEvent> events) {
        Set<Knowledgebase> sources = Sets.newHashSet();
        for (KnownEvent event : events) {
            sources.addAll(event.sources());
        }
        return Knowledgebase.toCommaSeparatedSourceString(sources);
    }
}
