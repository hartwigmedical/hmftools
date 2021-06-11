package com.hartwig.hmftools.common.purple.cnchromosome;

import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
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
        Set<String> chromosomeSet = Sets.newHashSet();
        List<CnPerChromosomeArmData> CnPerChromosomeArmData = Lists.newArrayList();
        for (Chromosome chr : refGenomeCoordinates.lengths().keySet()) {
            HumanChromosome chromosome = (HumanChromosome) chr;
            Map<ChromosomeArm, GenomeRegion> genomeRegion = determineArmRegion(chromosome, refGenomeCoordinates);

            for (Map.Entry<ChromosomeArm, GenomeRegion> entry : genomeRegion.entrySet()) {
                GenomeRegion region = entry.getValue();

                double copyNumberArm = 0;
                for (PurpleCopyNumber purpleCopyNumber : copyNumbers) {
                    Chromosome copyNumberChromosome = HumanChromosome.fromString(purpleCopyNumber.chromosome());
                    chromosomeSet.add(purpleCopyNumber.chromosome());

                    if (copyNumberChromosome.equals(chromosome) && region.overlaps(purpleCopyNumber)) {
                        double copyNumber = purpleCopyNumber.averageTumorCopyNumber();
                        long totalLengthSegment = purpleCopyNumber.bases();
                        copyNumberArm += (copyNumber * totalLengthSegment) / region.bases();
                    }
                }

                // if chromosome not present in copyNumber file, the chromosome arm isn't calculated
                if (chromosomeSet.contains(chromosome.toString())) {
                    CnPerChromosomeArmData.add(ImmutableCnPerChromosomeArmData.builder()
                            .chromosome(chromosome)
                            .chromosomeArm(entry.getKey())
                            .copyNumber(copyNumberArm)
                            .build());
                }
            }
        }
        return CnPerChromosomeArmData;
    }

    @NotNull
    private static Map<ChromosomeArm, GenomeRegion> determineArmRegion(@NotNull Chromosome chromosome,
            @NotNull RefGenomeCoordinates refGenomeCoordinates) {
        long centromerePos = refGenomeCoordinates.centromeres().get(chromosome);
        long chrLength = refGenomeCoordinates.lengths().get(chromosome);

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