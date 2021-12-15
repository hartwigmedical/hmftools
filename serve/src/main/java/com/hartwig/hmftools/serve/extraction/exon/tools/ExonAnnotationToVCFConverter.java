package com.hartwig.hmftools.serve.extraction.exon.tools;

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

        ServeConfig config = ServeLocalConfigProvider.create();
        IndexedFastaSequenceFile refSequence37 = new IndexedFastaSequenceFile(new File(config.refGenome37FastaFile()));

        String knownExonsTsv = System.getProperty("user.home") + "/hmf/tmp/KnownExons.SERVE.38.tsv";
        String outputVcf = System.getProperty("user.home") + "/hmf/tmp/exons.vcf.gz";
        GenerateAltBase altBaseGenerator = new GenerateAltBase(RefGenomeVersion.V37, refSequence37);

        List<KnownExon> exons = KnownExonFile.read(knownExonsTsv);
        LOGGER.info("The number of known exons in the known exon file is {}", exons.size());

        VariantContextWriter writer = VCFWriterFactory.openVCFWriter(outputVcf, uniqueSourcesString(exons));

        for (KnownExon exon : exons) {
            String chromosome = exon.annotation().chromosome();
            long start = exon.annotation().start() + 10; // remove the first 10 non-coding positions before exon
            long end = exon.annotation().end() - 10; // remove the last 10 non-coding positions after exon
            long middle = start + Math.round((end - start) / 2D); // take middle position of exon

            List<Long> positions = Lists.newArrayList(start, middle, end);
            for (long position : positions) {
                String refBaseOfPosition = altBaseGenerator.extractRefBaseAtGenomicPosition(chromosome, position);
                String randomAltBase = altBaseGenerator.createAltForRefBase(chromosome, position);

                writeVariantToVCF(writer,
                        chromosome,
                        position,
                        refBaseOfPosition,
                        randomAltBase,
                        exon.sources(),
                        exon.annotation().gene(),
                        exon.annotation().transcript(),
                        exon.annotation().rank());
            }
        }

        writer.close();

        LOGGER.info("All known exons are converted and written to '{}'", outputVcf);
    }

    private static void writeVariantToVCF(@NotNull VariantContextWriter writer, @NotNull String chromosome, long position,
            @NotNull String ref, @NotNull String alt, @NotNull Set<Knowledgebase> knowledgebases, @NotNull String gene, @NotNull String transcript,
            int exonRank) {
        List<Allele> alleles = Lists.newArrayList(Allele.create(ref, true), Allele.create(alt, false));

        VariantContext variant = new VariantContextBuilder().noGenotypes()
                .source("ExonChecker")
                .chr(chromosome)
                .start(position)
                .alleles(alleles)
                .computeEndFromAlleles(alleles, new Long(position).intValue())
                .attribute(VCFWriterFactory.INPUT_FIELD, KeyFormatter.toExonKey(gene, transcript, exonRank))
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
