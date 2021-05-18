package com.hartwig.hmftools.protect.cnchromosome;

import java.io.IOException;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeCoordinates;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumberFile;
import com.hartwig.hmftools.common.purple.segment.ChromosomeArm;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public final class CnPerChromosomeFactory {

    private static final Logger LOGGER = LogManager.getLogger(CnPerChromosomeFactory.class);

    private static final RefGenomeCoordinates REF_GENOME_COORDINATES = RefGenomeCoordinates.COORDS_37;

    private CnPerChromosomeFactory() {
    }

    @NotNull
    public static Map<ChromosomeArmKey, Double> fromPurpleSomaticCopynumberTsv(@NotNull String purpleSomaticCopynumberTsv) throws IOException {
        List<PurpleCopyNumber> copyNumbers = PurpleCopyNumberFile.read(purpleSomaticCopynumberTsv);
        return extractCnPerChromosomeArm(copyNumbers);
    }

    @NotNull
    private static Map<ChromosomeArmKey, Double> extractCnPerChromosomeArm(@NotNull List<PurpleCopyNumber> copyNumbers) {
        Map<ChromosomeArmKey, Double> cnPerChromosomePArm = determineCopyNumberArm(copyNumbers, ChromosomeArm.P_ARM);
        Map<ChromosomeArmKey, Double> cnPerChromosomeQArm = determineCopyNumberArm(copyNumbers, ChromosomeArm.Q_ARM);

        Map<ChromosomeArmKey, Double> cnPerChromosomeArm = Maps.newHashMap();
        cnPerChromosomeArm.putAll(cnPerChromosomePArm);
        cnPerChromosomeArm.putAll(cnPerChromosomeQArm);

        for (Map.Entry<ChromosomeArmKey, Double> entry : cnPerChromosomeArm.entrySet()) {
            LOGGER.info("{}: {}", entry.getKey(), entry.getValue());
        }

        return cnPerChromosomeArm;
    }

    @NotNull
    private static Map<ChromosomeArmKey, Double> determineCopyNumberArm(@NotNull List<PurpleCopyNumber> copyNumbers,
            @NotNull ChromosomeArm chromosomeArm) {
        Map<ChromosomeArmKey, Double> cnPerChromosomeArm = Maps.newHashMap();

        for (Chromosome chr : REF_GENOME_COORDINATES.lengths().keySet()) {
            HumanChromosome chromosome = (HumanChromosome) chr;
            ChromosomeArmKey key = new ChromosomeArmKey(chromosome, chromosomeArm);

            int chromosomeLength = getChromosomalArmLength(chromosome, chromosomeArm);
            double copyNumberArm = 0;
            for (PurpleCopyNumber purpleCopyNumber : copyNumbers) {
                Chromosome copyNumberChromosome = HumanChromosome.fromString(purpleCopyNumber.chromosome());
                if (chromosomeArm == ChromosomeArm.P_ARM) {
                    if (copyNumberChromosome.equals(chromosome) && purpleCopyNumber.end() < chromosomeLength) {
                        double copyNumber = purpleCopyNumber.averageTumorCopyNumber();
                        long totalLengthSegment = (purpleCopyNumber.end() - purpleCopyNumber.start()) + 1;
                        copyNumberArm += (copyNumber * totalLengthSegment) / chromosomeLength;
                    }
                } else if (chromosomeArm == ChromosomeArm.Q_ARM) {
                    if (copyNumberChromosome.equals(chromosome) && purpleCopyNumber.end() > chromosomeLength) {
                        double copyNumber = purpleCopyNumber.averageTumorCopyNumber();
                        long totalLengthSegment = (purpleCopyNumber.end() - purpleCopyNumber.start()) + 1;
                        copyNumberArm += (copyNumber * totalLengthSegment) / chromosomeLength;
                    }
                }
            }

            cnPerChromosomeArm.put(key, copyNumberArm);
        }

        return cnPerChromosomeArm;
    }

    private static int getChromosomalArmLength(@NotNull Chromosome chromosome, @NotNull ChromosomeArm armType) {
        Long centromerePos = REF_GENOME_COORDINATES.centromeres().get(chromosome);

        if (centromerePos == null) {
            return 0;
        }

        if (armType == ChromosomeArm.P_ARM) {
            return centromerePos.intValue();
        }

        int chrLength = REF_GENOME_COORDINATES.lengths().get(chromosome).intValue();

        return chrLength - centromerePos.intValue();
    }
}