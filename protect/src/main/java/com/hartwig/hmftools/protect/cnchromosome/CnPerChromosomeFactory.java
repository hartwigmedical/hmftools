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

            GenomeRegion genomeRegion = determineArmRegion(chromosome);
            ChromosomeArm chromosomeArms = determineChromosomeArm(chromosome);

            double copyNumberArm = 0;
            for (PurpleCopyNumber purpleCopyNumber : copyNumbers) {
                Chromosome copyNumberChromosome = HumanChromosome.fromString(purpleCopyNumber.chromosome());
                if (chromosomeArms == ChromosomeArm.P_ARM) {
                    if (copyNumberChromosome.equals(chromosome) && purpleCopyNumber.end() < genomeRegion.end()) {
                        double copyNumber = purpleCopyNumber.averageTumorCopyNumber();
                        long totalLengthSegment = (purpleCopyNumber.end() - purpleCopyNumber.start()) + 1;
                        copyNumberArm += (copyNumber * totalLengthSegment) / genomeRegion.bases();
                    }
                } else if (chromosomeArms == ChromosomeArm.Q_ARM) {
                    if (copyNumberChromosome.equals(chromosome) && purpleCopyNumber.end() > genomeRegion.start()) {
                        double copyNumber = purpleCopyNumber.averageTumorCopyNumber();
                        long totalLengthSegment = (purpleCopyNumber.end() - purpleCopyNumber.start()) + 1;
                        copyNumberArm += (copyNumber * totalLengthSegment) / genomeRegion.bases();
                    }
                }
            }

            cnPerChromosomeArm.put(new ChromosomeArmKey(chromosome, chromosomeArm), copyNumberArm);
        }

        return cnPerChromosomeArm;
    }

    @NotNull
    private static GenomeRegion determineArmRegion(@NotNull Chromosome chromosome) {

        // The smallest part of a chromosome is the P arm.
        Long centromerePos = REF_GENOME_COORDINATES.centromeres().get(chromosome);
        int chrLength = REF_GENOME_COORDINATES.lengths().get(chromosome).intValue();
        long endChr = chrLength - centromerePos;
        long startChr = centromerePos;

        if (startChr < endChr) {  //P arm
            return GenomeRegions.create(chromosome.toString(), 1, centromerePos);
        } else  {  // q arm
            return GenomeRegions.create(chromosome.toString(), centromerePos + 1, chrLength);
        }

    }

    @NotNull
    private static ChromosomeArm determineChromosomeArm(@NotNull Chromosome chromosome) {

        // The smallest part of a chromosome is the P arm.
        Long centromerePos = REF_GENOME_COORDINATES.centromeres().get(chromosome);
        int chrLength = REF_GENOME_COORDINATES.lengths().get(chromosome).intValue();
        long endChr = chrLength - centromerePos;
        long startChr = centromerePos;

        if (startChr < endChr) {  //P arm
            return ChromosomeArm.P_ARM;
        } else  {  // q arm
            return ChromosomeArm.Q_ARM;
        }

    }
}