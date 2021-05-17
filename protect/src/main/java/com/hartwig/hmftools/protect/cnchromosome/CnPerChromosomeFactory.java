package com.hartwig.hmftools.protect.cnchromosome;

import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.Set;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeCoordinates;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumberFile;
import com.hartwig.hmftools.common.purple.segment.ChromosomeArm;
import com.hartwig.hmftools.common.serve.Knowledgebase;
import com.hartwig.hmftools.protect.purple.PurpleDataLoader;
import com.hartwig.hmftools.serve.extraction.copynumber.CopyNumberFunctions;
import com.hartwig.hmftools.serve.extraction.copynumber.CopyNumberType;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public final class CnPerChromosomeFactory {
    private static final Logger LOGGER = LogManager.getLogger(PurpleDataLoader.class);

    private CnPerChromosomeFactory() {
    }

    @NotNull
    public static Map<CopyNumberKey, Double> fromPurpleSomaticCopynumberTsv(@NotNull String purpleSomaticCopynumberTsv) throws IOException {
        List<PurpleCopyNumber> copyNumbers = PurpleCopyNumberFile.read(purpleSomaticCopynumberTsv);
        return extractCnPerChromosomeArm(copyNumbers);
    }

    @NotNull
    private static Map<CopyNumberKey, Double> extractCnPerChromosomeArm(@NotNull List<PurpleCopyNumber> copyNumbers) {
        Map<CopyNumberKey, Double> cnPerChromosomeArm = Maps.newHashMap();
        Map<CopyNumberKey, Double> cnPerChromosomePArm = Maps.newHashMap(determineCopyNumberPArm(copyNumbers));
        Map<CopyNumberKey, Double> cnPerChromosomeQArm = Maps.newHashMap(determineCopyNumberQArm(copyNumbers));

        cnPerChromosomeArm.putAll(cnPerChromosomePArm);
        cnPerChromosomeArm.putAll(cnPerChromosomeQArm);

        for (Map.Entry<CopyNumberKey, Double> entry :cnPerChromosomeArm.entrySet()) {
            LOGGER.info(entry.getKey().chromosome + " " + entry.getKey().chromosomeArm + " " + entry.getValue());
        }

        return cnPerChromosomeArm;
    }

    private static int getChromosomalArmLength(@NotNull String chromosome, @NotNull ChromosomeArm armType) {
        RefGenomeCoordinates refGenome = RefGenomeCoordinates.COORDS_37;
        HumanChromosome chr = HumanChromosome.fromString(chromosome);

        Long centromerePos = refGenome.centromeres().get(chr);

        if (centromerePos == null) {
            return 0;
        }

        if (armType == ChromosomeArm.P_ARM) {
            return centromerePos.intValue();
        }

        int chrLength = refGenome.lengths().get(chr).intValue();

        return chrLength - centromerePos.intValue();
    }

    private static Map<CopyNumberKey, Double> determineCopyNumberPArm(@NotNull List<PurpleCopyNumber> copyNumbers) {
        RefGenomeCoordinates refGenome = RefGenomeCoordinates.COORDS_37;
        Map<CopyNumberKey, Double> cnPerChromosomeArm = Maps.newHashMap();

        for (Map.Entry<Chromosome, Long> entry : refGenome.lengths().entrySet()) {
            final String chromosome = entry.getKey().toString();

            CopyNumberKey key = new CopyNumberKey(chromosome, ChromosomeArm.P_ARM);

            int chromosomeLength = getChromosomalArmLength(chromosome, ChromosomeArm.P_ARM);
            double copyNumberArm = 0;
            for (PurpleCopyNumber purpleCopyNumber : copyNumbers) {
                if (purpleCopyNumber.chromosome().equals(chromosome) && purpleCopyNumber.end() < chromosomeLength) {
                    double copyNumber = purpleCopyNumber.averageTumorCopyNumber();
                    long totalLengthSegment = (purpleCopyNumber.end() - purpleCopyNumber.start()) + 1;
                    copyNumberArm += copyNumber * totalLengthSegment / chromosomeLength;
                }
            }

            cnPerChromosomeArm.put(key, copyNumberArm);

        }

        return cnPerChromosomeArm;
    }

    private static Map<CopyNumberKey, Double> determineCopyNumberQArm(@NotNull List<PurpleCopyNumber> copyNumbers) {
        RefGenomeCoordinates refGenome = RefGenomeCoordinates.COORDS_37;
        Map<CopyNumberKey, Double> cnPerChromosomeArm = Maps.newHashMap();

        for (Map.Entry<Chromosome, Long> entry : refGenome.lengths().entrySet()) {
            final String chromosome = entry.getKey().toString();

            CopyNumberKey key = new CopyNumberKey(chromosome, ChromosomeArm.Q_ARM);

            int chromosomeLength = getChromosomalArmLength(chromosome, ChromosomeArm.Q_ARM);
            double copyNumberArm = 0;
            for (PurpleCopyNumber purpleCopyNumber : copyNumbers) {
                if (purpleCopyNumber.chromosome().equals(chromosome) && purpleCopyNumber.end() > chromosomeLength) {
                    double copyNumber = purpleCopyNumber.averageTumorCopyNumber();
                    long totalLengthSegment = (purpleCopyNumber.end() - purpleCopyNumber.start()) + 1;
                    copyNumberArm += copyNumber * totalLengthSegment / chromosomeLength;
                }
            }

            cnPerChromosomeArm.put(key, copyNumberArm);

        }

        return cnPerChromosomeArm;
    }

    public static class CopyNumberKey {

        @NotNull
        private final String chromosome;
        @NotNull
        private final ChromosomeArm chromosomeArm;

        public CopyNumberKey(@NotNull final String chromosome, @NotNull final ChromosomeArm chromosomeArm) {
            this.chromosome = chromosome;
            this.chromosomeArm = chromosomeArm;
        }

        @NotNull
        public String chromosome() {
            return chromosome;
        }

        @NotNull
        public ChromosomeArm chromosomeArm() {
            return chromosomeArm;
        }

        @Override
        public boolean equals(final Object o) {
            if (this == o) {
                return true;
            }
            if (o == null || getClass() != o.getClass()) {
                return false;
            }
            final CopyNumberKey that = (CopyNumberKey) o;
            return chromosome.equals(that.chromosome) && chromosomeArm == that.chromosomeArm;
        }

        @Override
        public int hashCode() {
            return Objects.hash(chromosome, chromosomeArm);
        }
    }
}