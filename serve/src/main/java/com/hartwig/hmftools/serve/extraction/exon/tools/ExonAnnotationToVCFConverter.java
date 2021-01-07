package com.hartwig.hmftools.serve.extraction.exon.tools;

import java.io.IOException;
import java.util.List;
import java.util.Set;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.serve.Knowledgebase;
import com.hartwig.hmftools.serve.extraction.exon.KnownExon;
import com.hartwig.hmftools.serve.extraction.exon.KnownExonFile;
import com.hartwig.hmftools.serve.extraction.util.GenerateAltBase;
import com.hartwig.hmftools.serve.extraction.util.KeyFormatter;
import com.hartwig.hmftools.serve.extraction.util.VCFWriter;

import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.core.config.Configurator;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;

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

        VariantContextWriter writer = VCFWriter.generateVCFWriter(outputFile);

        for (KnownExon exon : exons) {
            String chromosome = exon.annotation().chromosome();
            Long start = exon.annotation().start() + 5; // remove the first 5 slice position before exon
            Long end = exon.annotation().end() - 5; // remove the last 5 splice postion after exon
            String gene = exon.annotation().gene();

            Long bewtweenNumber = start + Math.round((end - start)/2);
            List<Long> genomicPositions = Lists.newArrayList(start, bewtweenNumber, end);

            for (Long genomicPosition : genomicPositions) {
                String extractRefBaseOfPosition = GenerateAltBase.extractRefBaseOfGenomicPosition(chromosome, genomicPosition);
                String randomAltBase = GenerateAltBase.createAltOfRefBase(chromosome, genomicPosition);

                extactAnnotationVariantExonIndex(extractRefBaseOfPosition,
                        randomAltBase,
                        chromosome,
                        genomicPosition,
                        exon.sources(),
                        gene,
                        exon.annotation().exonIndex(),
                        exon.annotation().transcript(),
                        writer);

            }
        }
        writer.close();
        LOGGER.info("All exons are checked!");

        LOGGER.info("Done!");
    }



    private static void extactAnnotationVariantExonIndex(@NotNull String extractRefBaseOfPosition, @NotNull String randomAltBase,
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

    private static void generateVcfFileOfGenomicPosition(@NotNull String extractRefBaseOfPosition, @NotNull String randomAltBase,
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
                .attribute("input", KeyFormatter.toExonKey(gene, transcript, Integer.toString(exonIndex)))
                .make();

        LOGGER.debug(" Writing variant to VCF file'{}'", variantContext);
        writer.add(variantContext);
    }
}
