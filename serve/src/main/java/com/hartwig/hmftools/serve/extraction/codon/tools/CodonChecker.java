package com.hartwig.hmftools.serve.extraction.codon.tools;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Random;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.serve.Knowledgebase;
import com.hartwig.hmftools.common.variant.snpeff.SnpEffAnnotation;
import com.hartwig.hmftools.common.variant.snpeff.SnpEffAnnotationFactory;
import com.hartwig.hmftools.serve.extraction.codon.KnownCodon;
import com.hartwig.hmftools.serve.extraction.codon.KnownCodonFile;
import com.hartwig.hmftools.serve.extraction.hotspot.ProteinKeyFormatter;
import com.hartwig.hmftools.serve.extraction.util.MutationTypeFilter;

import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.core.config.Configurator;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFCodec;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

public class CodonChecker {
    private static final Logger LOGGER = LogManager.getLogger(CodonChecker.class);
    private static final boolean LOG_DEBUG = true;

    public static void main(String[] args) throws IOException {
        LOGGER.info("Running SERVE codon checker");

        if (LOG_DEBUG) {
            Configurator.setRootLevel(Level.DEBUG);
        }

        Random random = new Random();
        String knownCodonsTsv = System.getProperty("user.home") + "/hmf/tmp/serve/KnownCodons.SERVE.37.tsv";
        String outputFile = System.getProperty("user.home") + "/hmf/tmp/test.vcf.gz";

        List<KnownCodon> codons = KnownCodonFile.read(knownCodonsTsv);

        LOGGER.info("The size of the file is {}", codons.size());

        List<String> bases = Lists.newArrayList();
        String randomAltBase = null;
        for (KnownCodon codon : codons) {
            Long start = codon.annotation().start();
            Long end = codon.annotation().end();
            Long bewtweenNumber = end > start ? start + 1 : start - 1;
            List<Long> genomicPositions = Lists.newArrayList(start, bewtweenNumber, end);
            String chromosome = codon.annotation().chromosome();
            Integer codonIndex = codon.annotation().codonIndex();
            MutationTypeFilter mutationTypeFilter = codon.annotation().mutationType();

            for (Long genomicPosition : genomicPositions) {
                String extractRefBaseOfPosition = extractRefBaseOfGenomicPosition(chromosome, genomicPosition);

                if (extractRefBaseOfPosition.equals("A")) {
                    bases = Lists.newArrayList("C", "T", "G");
                    randomAltBase = bases.get(random.nextInt(3));
                } else if (extractRefBaseOfPosition.equals("C")) {
                    bases = Lists.newArrayList("A", "T", "G");
                    randomAltBase = bases.get(random.nextInt(3));
                } else if (extractRefBaseOfPosition.equals("T")) {
                    bases = Lists.newArrayList("C", "A", "G");
                    randomAltBase = bases.get(random.nextInt(3));
                } else if (extractRefBaseOfPosition.equals("G")) {
                    bases = Lists.newArrayList("C", "T", "A");
                    randomAltBase = bases.get(random.nextInt(3));
                }

                Integer annotatedVariantCodonIndex = extactAnnotationVariantCodonIndex(extractRefBaseOfPosition,
                        randomAltBase,
                        chromosome,
                        genomicPosition,
                        mutationTypeFilter,
                        codon.sources(),
                        codon.annotation().gene(),
                        codon.annotation().proteinAnnotation(),
                        codon.annotation().transcript(),
                        outputFile);

                if (!codonIndex.equals(annotatedVariantCodonIndex)) {
                    LOGGER.warn("Condon index of SERVE {} are not equals for annotated codon index {}",
                            codonIndex,
                            annotatedVariantCodonIndex);
                }
            }
        }

        LOGGER.info("All codons are checked!");
        LOGGER.info("Done!");
    }

    private static String extractRefBaseOfGenomicPosition(@Nullable String chromosome, Long genomicPosition) throws IOException {
        IndexedFastaSequenceFile fastaSequenceFile =
                new IndexedFastaSequenceFile(new File("/Users/liekeschoenmaker/hmf/refgenome/Homo_sapiens.GRCh37.GATK.illumina.fasta"));
        return fastaSequenceFile.getSubsequenceAt(chromosome, genomicPosition, genomicPosition).getBaseString();

    }

    private static Integer extactAnnotationVariantCodonIndex(@Nullable String extractRefBaseOfPosition, @Nullable String randomAltBase,
            @Nullable String chromosome, Long position, @NotNull MutationTypeFilter mutationTypeFilter,
            @NotNull Set<Knowledgebase> knowledgebases, @NotNull String gene, @NotNull String proteinAnnotation, @NotNull String transcript,
            @NotNull String outputFile) throws IOException {

        generateVcfFileOfGenomicPosition(extractRefBaseOfPosition,
                randomAltBase,
                chromosome,
                position,
                knowledgebases,
                gene,
                proteinAnnotation,
                transcript,
                outputFile);

        AbstractFeatureReader<VariantContext, LineIterator> reader =
                AbstractFeatureReader.getFeatureReader(outputFile, new VCFCodec(), false);
        for (VariantContext variant : reader.iterator()) {

            List<SnpEffAnnotation> annotations = SnpEffAnnotationFactory.fromContext(variant);
            LOGGER.info(annotations); //TODO fix list are empty
        }

        //TODO check mutationTyepfilter
        return 1;
    }

    private static void generateVcfFileOfGenomicPosition(@Nullable String extractRefBaseOfPosition, @Nullable String randomAltBase,
            @Nullable String chromosome, Long position, @NotNull Set<Knowledgebase> knowledgebases, @NotNull String gene,
            @NotNull String proteinAnnotation, @NotNull String transcript, @NotNull String outputFile) {
        VariantContextWriter writer = new VariantContextWriterBuilder().setOutputFile(outputFile)
                .setOutputFileType(VariantContextWriterBuilder.OutputType.BLOCK_COMPRESSED_VCF)
                .modifyOption(Options.INDEX_ON_THE_FLY, false)
                .build();

        VCFHeader header = new VCFHeader(Sets.newHashSet(), Lists.newArrayList());
        header.addMetaDataLine(new VCFInfoHeaderLine("input", VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "input"));
        header.addMetaDataLine(new VCFInfoHeaderLine("source",
                VCFHeaderLineCount.UNBOUNDED,
                VCFHeaderLineType.String,
                "sources [" + knowledgebases + "]"));

        writer.writeHeader(header);

        List<Allele> hotspotAlleles =
                Lists.newArrayList(Allele.create(extractRefBaseOfPosition, true), Allele.create(randomAltBase, false));

        VariantContext variantContext = new VariantContextBuilder().noGenotypes()
                .source("SERVE")
                .chr(chromosome)
                .start(position)
                .alleles(hotspotAlleles)
                .computeEndFromAlleles(hotspotAlleles, new Long(position).intValue())
                .attribute("source", Knowledgebase.commaSeparatedSourceString(knowledgebases))
                .attribute("input", ProteinKeyFormatter.toProteinKey(gene, transcript, proteinAnnotation))
                .make();

        LOGGER.debug(" Writing variant '{}'", variantContext);
        writer.add(variantContext);

        writer.close();
    }

}
