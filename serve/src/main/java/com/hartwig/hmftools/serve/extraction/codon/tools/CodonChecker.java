package com.hartwig.hmftools.serve.extraction.codon.tools;

import java.io.IOException;
import java.util.List;
import java.util.Random;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.serve.extraction.codon.KnownCodon;
import com.hartwig.hmftools.serve.extraction.codon.KnownCodonFile;
import com.hartwig.hmftools.serve.extraction.util.MutationTypeFilter;

import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.core.config.Configurator;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

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
        List<KnownCodon> codons = KnownCodonFile.read(knownCodonsTsv);

        LOGGER.info("The size of the file is {}", codons.size());

        List<String> bases = Lists.newArrayList();
        String randomAltBase = null;
        for (KnownCodon codon : codons) {
            Long start = codon.annotation().start();
            Long end = codon.annotation().end();
            Long bewtweenNumber = end > start ? start + 1 : start -1;
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

                Integer annotatedVariantCodonIndex =
                        extactAnnotationVariantCodonIndex(extractRefBaseOfPosition, randomAltBase, chromosome, genomicPosition, mutationTypeFilter);

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

    private static String extractRefBaseOfGenomicPosition(@Nullable String chromosome, Long genomicPosition) {
        String genomicPositionBase = chromosome + ":" + genomicPosition + "-" + genomicPosition;
        //TODO extract base with genomic posiiton
        return genomicPositionBase;
    }

    private static Integer extactAnnotationVariantCodonIndex(@Nullable String extractRefBaseOfPosition, @Nullable String randomAltBase,
            @Nullable String chromosome, Long position, @NotNull MutationTypeFilter mutationTypeFilter) {
        //TODO extract codon annotation

        //TODO check mutationTyepfilter
        return 1;
    }
}
