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

        return cnPerChromosomeArm;
    }

    @NotNull
    private static Map<ChromosomeArmKey, Double> determineCopyNumberArm(@NotNull List<PurpleCopyNumber> copyNumbers) {
        Map<ChromosomeArmKey, Double> cnPerChromosomeArm = Maps.newHashMap();

        for (Chromosome chr : REF_GENOME_COORDINATES.lengths().keySet()) {
            HumanChromosome chromosome = (HumanChromosome) chr;
            Map<ChromosomeArm, GenomeRegion> genomeRegion = determineArmRegion(chromosome);

            for (Map.Entry<ChromosomeArm, GenomeRegion> entry : genomeRegion.entrySet()) {
                GenomeRegion region = entry.getValue();

                double copyNumberArm = 0;
                for (PurpleCopyNumber purpleCopyNumber : copyNumbers) {
                    Chromosome copyNumberChromosome = HumanChromosome.fromString(purpleCopyNumber.chromosome());

                    if (copyNumberChromosome.equals(chromosome) && region.overlaps(purpleCopyNumber)) {
                        double copyNumber = purpleCopyNumber.averageTumorCopyNumber();
                        long totalLengthSegment = purpleCopyNumber.bases();
                        copyNumberArm += (copyNumber * totalLengthSegment) / region.bases();
                    }
                }
                cnPerChromosomeArm.put(new ChromosomeArmKey(chromosome, entry.getKey()), copyNumberArm);
            }
        }

        return cnPerChromosomeArm;
    }

    @NotNull
    private static Map<ChromosomeArm, GenomeRegion> determineArmRegion(@NotNull Chromosome chromosome) {
        long centromerePos = REF_GENOME_COORDINATES.centromeres().get(chromosome);
        long chrLength = REF_GENOME_COORDINATES.lengths().get(chromosome);

        long lengthPastCentromere = chrLength - (centromerePos + 1);
        long lengthBeforeCentromere = centromerePos;

        Map<ChromosomeArm, GenomeRegion> chromosomeArmGenomeRegionMap = Maps.newHashMap();

        // The smallest part of a chromosome is the P arm.
        if (lengthBeforeCentromere < lengthPastCentromere) {
            chromosomeArmGenomeRegionMap.put(ChromosomeArm.P_ARM, GenomeRegions.create(chromosome.toString(), 1, centromerePos));
            chromosomeArmGenomeRegionMap.put(ChromosomeArm.Q_ARM,
                    GenomeRegions.create(chromosome.toString(), centromerePos + 1, chrLength));
        } else {
            chromosomeArmGenomeRegionMap.put(ChromosomeArm.Q_ARM, GenomeRegions.create(chromosome.toString(), 1, centromerePos));
            chromosomeArmGenomeRegionMap.put(ChromosomeArm.P_ARM,
                    GenomeRegions.create(chromosome.toString(), centromerePos + 1, chrLength));
        }
        return chromosomeArmGenomeRegionMap;
    }

}