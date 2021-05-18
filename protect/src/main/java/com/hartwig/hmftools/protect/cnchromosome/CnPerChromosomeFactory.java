package com.hartwig.hmftools.protect.cnchromosome;

import java.io.IOException;
import java.util.List;
import java.util.Map;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeCoordinates;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.genome.region.GenomeRegions;
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
    public static Map<ChromosomeArmKey, Double> fromPurpleSomaticCopynumberTsv(@NotNull String purpleSomaticCopynumberTsv)
            throws IOException {
        List<PurpleCopyNumber> copyNumbers = PurpleCopyNumberFile.read(purpleSomaticCopynumberTsv);
        return extractCnPerChromosomeArm(copyNumbers);
    }

    @NotNull
    @VisibleForTesting
    static Map<ChromosomeArmKey, Double> extractCnPerChromosomeArm(@NotNull List<PurpleCopyNumber> copyNumbers) {
        Map<ChromosomeArmKey, Double> cnPerChromosomeArm = determineCopyNumberArm(copyNumbers);

        for (Map.Entry<ChromosomeArmKey, Double> entry : cnPerChromosomeArm.entrySet()) {
            LOGGER.info("{}: {}", entry.getKey(), entry.getValue());
        }

        return cnPerChromosomeArm;
    }

    @NotNull
    private static Map<ChromosomeArmKey, Double> determineCopyNumberArm(@NotNull List<PurpleCopyNumber> copyNumbers) {
        Map<ChromosomeArmKey, Double> cnPerChromosomeArm = Maps.newHashMap();

        for (Chromosome chr : REF_GENOME_COORDINATES.lengths().keySet()) {
            HumanChromosome chromosome = (HumanChromosome) chr;
            Map<ChromosomeArm, GenomeRegion> genomeRegion = determineArmRegion(chromosome);
            double copyNumberArm = 0;

            for (Map.Entry<ChromosomeArm, GenomeRegion> entry : genomeRegion.entrySet()) {
                for (PurpleCopyNumber purpleCopyNumber : copyNumbers) {
                    Chromosome copyNumberChromosome = HumanChromosome.fromString(purpleCopyNumber.chromosome());
                    if (entry.getKey() == ChromosomeArm.P_ARM) {
                        if (copyNumberChromosome.equals(chromosome) && purpleCopyNumber.end() <= entry.getValue().end()) {
                            double copyNumber = purpleCopyNumber.averageTumorCopyNumber();
                            long totalLengthSegment = (purpleCopyNumber.end() - purpleCopyNumber.start()) + 1;
                            copyNumberArm += (copyNumber * totalLengthSegment) / entry.getValue().bases();
                        }
                    } else if (entry.getKey() == ChromosomeArm.Q_ARM) {
                        if (copyNumberChromosome.equals(chromosome) && purpleCopyNumber.end() <= entry.getValue().end()) {
                            double copyNumber = purpleCopyNumber.averageTumorCopyNumber();
                            long totalLengthSegment = (purpleCopyNumber.end() - purpleCopyNumber.start()) + 1;
                            copyNumberArm += (copyNumber * totalLengthSegment) / entry.getValue().bases();
                        }
                    }
                }
                cnPerChromosomeArm.put(new ChromosomeArmKey(chromosome, entry.getKey()), copyNumberArm);
            }

        }

        return cnPerChromosomeArm;
    }

    @NotNull
    private static Map<ChromosomeArm, GenomeRegion> determineArmRegion(@NotNull Chromosome chromosome) {

        Long centromerePos = REF_GENOME_COORDINATES.centromeres().get(chromosome);
        int chrLength = REF_GENOME_COORDINATES.lengths().get(chromosome).intValue();
        long endChr = chrLength - (centromerePos +1);
        long startChr = centromerePos;
        Map<ChromosomeArm, GenomeRegion> chromosomeArmGenomeRegionMap = Maps.newHashMap();

        // The smallest part of a chromosome is the P arm.
        if (startChr < endChr) {
            chromosomeArmGenomeRegionMap.put(ChromosomeArm.P_ARM, GenomeRegions.create(chromosome.toString(), 1, centromerePos));
            chromosomeArmGenomeRegionMap.put(ChromosomeArm.Q_ARM, GenomeRegions.create(chromosome.toString(), centromerePos + 1, chrLength));
        } else  {
            chromosomeArmGenomeRegionMap.put(ChromosomeArm.Q_ARM, GenomeRegions.create(chromosome.toString(), 1, centromerePos));
            chromosomeArmGenomeRegionMap.put(ChromosomeArm.P_ARM, GenomeRegions.create(chromosome.toString(), centromerePos + 1, chrLength));
        }
        return chromosomeArmGenomeRegionMap;
    }

}