package com.hartwig.hmftools.common.purple.cnchromosome;

import java.io.IOException;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeCoordinates;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.genome.region.GenomeRegions;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumberFile;
import com.hartwig.hmftools.common.purple.segment.ChromosomeArm;

import org.jetbrains.annotations.NotNull;

public final class CnPerChromosomeFactory {

    private CnPerChromosomeFactory() {
    }

    @NotNull
    public static List<CnPerChromosomeArmData> fromPurpleSomaticCopynumberTsv(@NotNull String purpleSomaticCopynumberTsv,
            @NotNull RefGenomeCoordinates refGenomeCoordinates) throws IOException {
        List<PurpleCopyNumber> copyNumbers = PurpleCopyNumberFile.read(purpleSomaticCopynumberTsv);
        return extractCnPerChromosomeArm(copyNumbers, refGenomeCoordinates);
    }

    @NotNull
    public static List<CnPerChromosomeArmData> extractCnPerChromosomeArm(@NotNull List<PurpleCopyNumber> copyNumbers,
            @NotNull RefGenomeCoordinates refGenomeCoordinates) {
        List<CnPerChromosomeArmData> cnPerChromosomeArmData = Lists.newArrayList();
        for (Chromosome chr : refGenomeCoordinates.lengths().keySet()) {
            HumanChromosome chromosome = (HumanChromosome) chr;
            Map<ChromosomeArm, GenomeRegion> genomeRegion = determineArmRegions(chromosome, refGenomeCoordinates);

            for (Map.Entry<ChromosomeArm, GenomeRegion> entry : genomeRegion.entrySet()) {
                GenomeRegion arm = entry.getValue();

                double copyNumberArm = 0;
                for (PurpleCopyNumber purpleCopyNumber : copyNumbers) {
                    Chromosome copyNumberChromosome = HumanChromosome.fromString(purpleCopyNumber.chromosome());

                    if (copyNumberChromosome.equals(chromosome) && arm.overlaps(purpleCopyNumber)) {
                        double copyNumber = purpleCopyNumber.averageTumorCopyNumber();
                        long totalLengthSegment = purpleCopyNumber.bases();
                        copyNumberArm += (copyNumber * totalLengthSegment) / arm.bases();
                    }
                }

                if (copyNumberArm > 0) {
                    cnPerChromosomeArmData.add(ImmutableCnPerChromosomeArmData.builder()
                            .chromosome(chromosome)
                            .chromosomeArm(entry.getKey())
                            .copyNumber(copyNumberArm)
                            .build());
                }
            }
        }
        return cnPerChromosomeArmData;
    }

    @NotNull
    private static Map<ChromosomeArm, GenomeRegion> determineArmRegions(@NotNull Chromosome chromosome,
            @NotNull RefGenomeCoordinates refGenomeCoordinates) {
        long centromerePos = refGenomeCoordinates.centromeres().get(chromosome);
        long chrLength = refGenomeCoordinates.lengths().get(chromosome);

        Map<ChromosomeArm, GenomeRegion> chromosomeArmGenomeRegionMap = Maps.newHashMap();

        GenomeRegion partBeforeCentromere = GenomeRegions.create(chromosome.toString(), 1, centromerePos);
        GenomeRegion partAfterCentromere = GenomeRegions.create(chromosome.toString(), centromerePos + 1, chrLength);
        // The smallest part of a chromosome is the P arm.
        if (partBeforeCentromere.bases() < partAfterCentromere.bases()) {
            chromosomeArmGenomeRegionMap.put(ChromosomeArm.P_ARM, partBeforeCentromere);
            chromosomeArmGenomeRegionMap.put(ChromosomeArm.Q_ARM, partAfterCentromere);
        } else {
            chromosomeArmGenomeRegionMap.put(ChromosomeArm.P_ARM, partAfterCentromere);
            chromosomeArmGenomeRegionMap.put(ChromosomeArm.Q_ARM, partBeforeCentromere);
        }
        return chromosomeArmGenomeRegionMap;
    }
}