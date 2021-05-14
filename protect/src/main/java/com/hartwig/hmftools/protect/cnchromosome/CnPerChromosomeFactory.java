package com.hartwig.hmftools.protect.cnchromosome;

import java.io.IOException;
import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeCoordinates;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumberFile;
import com.hartwig.hmftools.common.purple.segment.ChromosomeArm;

import org.jetbrains.annotations.NotNull;

public final class CnPerChromosomeFactory {

    private CnPerChromosomeFactory() {
    }

    @NotNull
    public static CnPerChromosome fromPurpleSomaticCopynumberTsv(@NotNull String purpleSomaticCopynumberTsv) throws IOException {
        List<PurpleCopyNumber> copyNumbers = PurpleCopyNumberFile.read(purpleSomaticCopynumberTsv);

        return extractCnPerChromosomeArm(copyNumbers);
    }

    @NotNull
    static CnPerChromosome extractCnPerChromosomeArm(@NotNull List<PurpleCopyNumber> copyNumbers) {
        return ImmutableCnPerChromosome.builder()
                .chr1p(determineCopyNumberPArm(copyNumbers))
                .chr1q(determineCopyNumberQArm(copyNumbers))
                .chr2p(determineCopyNumberPArm(copyNumbers))
                .chr2q(determineCopyNumberQArm(copyNumbers))
                .chr3p(determineCopyNumberPArm(copyNumbers))
                .chr3q(determineCopyNumberQArm(copyNumbers))
                .chr4p(determineCopyNumberPArm(copyNumbers))
                .chr4q(determineCopyNumberQArm(copyNumbers))
                .chr5p(determineCopyNumberPArm(copyNumbers))
                .chr5q(determineCopyNumberQArm(copyNumbers))
                .chr6p(determineCopyNumberPArm(copyNumbers))
                .chr6q(determineCopyNumberQArm(copyNumbers))
                .chr7p(determineCopyNumberPArm(copyNumbers))
                .chr7q(determineCopyNumberQArm(copyNumbers))
                .chr8p(determineCopyNumberPArm(copyNumbers))
                .chr8q(determineCopyNumberQArm(copyNumbers))
                .chr9p(determineCopyNumberPArm(copyNumbers))
                .chr9q(determineCopyNumberQArm(copyNumbers))
                .chr10p(determineCopyNumberPArm(copyNumbers))
                .chr10q(determineCopyNumberQArm(copyNumbers))
                .chr11p(determineCopyNumberPArm(copyNumbers))
                .chr11q(determineCopyNumberQArm(copyNumbers))
                .chr12p(determineCopyNumberPArm(copyNumbers))
                .chr12q(determineCopyNumberQArm(copyNumbers))
                .chr13p(determineCopyNumberPArm(copyNumbers))
                .chr13q(determineCopyNumberQArm(copyNumbers))
                .chr14p(determineCopyNumberPArm(copyNumbers))
                .chr14q(determineCopyNumberQArm(copyNumbers))
                .chr15p(determineCopyNumberPArm(copyNumbers))
                .chr15q(determineCopyNumberQArm(copyNumbers))
                .chr16p(determineCopyNumberPArm(copyNumbers))
                .chr16q(determineCopyNumberQArm(copyNumbers))
                .chr17p(determineCopyNumberPArm(copyNumbers))
                .chr17q(determineCopyNumberQArm(copyNumbers))
                .chr18p(determineCopyNumberPArm(copyNumbers))
                .chr18q(determineCopyNumberQArm(copyNumbers))
                .chr19p(determineCopyNumberPArm(copyNumbers))
                .chr19q(determineCopyNumberQArm(copyNumbers))
                .chr20p(determineCopyNumberPArm(copyNumbers))
                .chr20q(determineCopyNumberQArm(copyNumbers))
                .chr21p(determineCopyNumberPArm(copyNumbers))
                .chr21q(determineCopyNumberQArm(copyNumbers))
                .chr22p(determineCopyNumberPArm(copyNumbers))
                .chr22q(determineCopyNumberQArm(copyNumbers))
                .chrXp(determineCopyNumberPArm(copyNumbers))
                .chrXq(determineCopyNumberQArm(copyNumbers))
                .chrYp(determineCopyNumberPArm(copyNumbers))
                .chrYq(determineCopyNumberQArm(copyNumbers))
                .build();
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

    private static double determineCopyNumberPArm(@NotNull List<PurpleCopyNumber> copyNumbers) {
        RefGenomeCoordinates refGenome = RefGenomeCoordinates.COORDS_37;

        for(Map.Entry<Chromosome,Long> entry : refGenome.lengths().entrySet()) {
            final String chromosome = entry.getKey().toString();
            int chromosomeLength = getChromosomalArmLength(chromosome, ChromosomeArm.P_ARM);
            double copyNumber = 0;
            double copyNumberArm = 0;
            for (PurpleCopyNumber purpleCopyNumber : copyNumbers) {
                if (purpleCopyNumber.chromosome().equals(chromosome) && purpleCopyNumber.end() < chromosomeLength) {
                    copyNumber += purpleCopyNumber.averageObservedBAF();
                    long totalLengthSegment = (purpleCopyNumber.end() - purpleCopyNumber.start()) - 1;
                    copyNumber += copyNumber * totalLengthSegment / chromosomeLength;
                }
            }

        }

        return 0.0;
    }

    private static double determineCopyNumberQArm(@NotNull List<PurpleCopyNumber> copyNumbers) {
        RefGenomeCoordinates refGenome = RefGenomeCoordinates.COORDS_37;

        for(Map.Entry<Chromosome,Long> entry : refGenome.lengths().entrySet()) {
            final String chromosome = entry.getKey().toString();
            int chromosomeLength = getChromosomalArmLength(chromosome, ChromosomeArm.Q_ARM);
            double copyNumber = 0;
            double copyNumberArm = 0;
            for (PurpleCopyNumber purpleCopyNumber : copyNumbers) {
                if (purpleCopyNumber.chromosome().equals(chromosome) && purpleCopyNumber.end() > chromosomeLength) {
                    copyNumber += purpleCopyNumber.averageObservedBAF();
                    long totalLengthSegment = (purpleCopyNumber.end() - purpleCopyNumber.start()) - 1;
                    copyNumber += copyNumber * totalLengthSegment / chromosomeLength;
                }
            }

        }

        return 0.0;
    }
}