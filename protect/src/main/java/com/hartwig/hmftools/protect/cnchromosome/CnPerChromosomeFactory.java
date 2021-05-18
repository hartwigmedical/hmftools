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

    private CnPerChromosomeFactory() {
    }

    @NotNull
    public static Map<ChromosomeArmKey, Double> fromPurpleSomaticCopynumberTsv(@NotNull String purpleSomaticCopynumberTsv) throws IOException {
        List<PurpleCopyNumber> copyNumbers = PurpleCopyNumberFile.read(purpleSomaticCopynumberTsv);
        return extractCnPerChromosomeArm(copyNumbers);
    }

    @NotNull
    private static Map<ChromosomeArmKey, Double> extractCnPerChromosomeArm(@NotNull List<PurpleCopyNumber> copyNumbers) {
        Map<ChromosomeArmKey, Double> cnPerChromosomeArm = Maps.newHashMap();
        Map<ChromosomeArmKey, Double> cnPerChromosomePArm = Maps.newHashMap(determineCopyNumberArm(copyNumbers, ChromosomeArm.P_ARM));
        Map<ChromosomeArmKey, Double> cnPerChromosomeQArm = Maps.newHashMap(determineCopyNumberArm(copyNumbers, ChromosomeArm.Q_ARM));

        cnPerChromosomeArm.putAll(cnPerChromosomePArm);
        cnPerChromosomeArm.putAll(cnPerChromosomeQArm);

        for (Map.Entry<ChromosomeArmKey, Double> entry : cnPerChromosomeArm.entrySet()) {
            LOGGER.info("{}: {}", entry.getKey(), entry.getValue());
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

    @NotNull
    private static Map<ChromosomeArmKey, Double> determineCopyNumberArm(@NotNull List<PurpleCopyNumber> copyNumbers,
            @NotNull ChromosomeArm chromosomeArm) {
        RefGenomeCoordinates refGenome = RefGenomeCoordinates.COORDS_37;
        Map<ChromosomeArmKey, Double> cnPerChromosomeArm = Maps.newHashMap();

        for (Map.Entry<Chromosome, Long> entry : refGenome.lengths().entrySet()) {
            String chromosome = entry.getKey().toString();

            ChromosomeArmKey key = new ChromosomeArmKey(chromosome, chromosomeArm);

            int chromosomeLength = getChromosomalArmLength(chromosome, chromosomeArm);
            double copyNumberArm = 0;
            for (PurpleCopyNumber purpleCopyNumber : copyNumbers) {
                if (chromosomeArm == ChromosomeArm.P_ARM) {
                    if (purpleCopyNumber.chromosome().equals(chromosome) && purpleCopyNumber.end() < chromosomeLength) {
                        double copyNumber = purpleCopyNumber.averageTumorCopyNumber();
                        long totalLengthSegment = (purpleCopyNumber.end() - purpleCopyNumber.start()) + 1;
                        copyNumberArm += (copyNumber * totalLengthSegment) / chromosomeLength;
                    }
                } else if (chromosomeArm == ChromosomeArm.Q_ARM) {
                    if (purpleCopyNumber.chromosome().equals(chromosome) && purpleCopyNumber.end() > chromosomeLength) {
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

}