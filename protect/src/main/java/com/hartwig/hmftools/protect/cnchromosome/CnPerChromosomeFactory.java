package com.hartwig.hmftools.protect.cnchromosome;

import java.io.IOException;
import java.util.List;

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
                .chr1p(determineCopyNumberPArm(copyNumbers, "1"))
                .chr1q(determineCopyNumberQArm(copyNumbers, "1"))
                .chr2p(determineCopyNumberPArm(copyNumbers, "2"))
                .chr2q(determineCopyNumberQArm(copyNumbers, "2"))
                .chr3p(determineCopyNumberPArm(copyNumbers, "3"))
                .chr3q(determineCopyNumberQArm(copyNumbers, "3"))
                .chr4p(determineCopyNumberPArm(copyNumbers, "4"))
                .chr4q(determineCopyNumberQArm(copyNumbers, "4"))
                .chr5p(determineCopyNumberPArm(copyNumbers, "5"))
                .chr5q(determineCopyNumberQArm(copyNumbers, "5"))
                .chr6p(determineCopyNumberPArm(copyNumbers, "6"))
                .chr6q(determineCopyNumberQArm(copyNumbers, "6"))
                .chr7p(determineCopyNumberPArm(copyNumbers, "7"))
                .chr7q(determineCopyNumberQArm(copyNumbers, "7"))
                .chr8p(determineCopyNumberPArm(copyNumbers, "8"))
                .chr8q(determineCopyNumberQArm(copyNumbers, "8"))
                .chr9p(determineCopyNumberPArm(copyNumbers, "9"))
                .chr9q(determineCopyNumberQArm(copyNumbers, "9"))
                .chr10p(determineCopyNumberPArm(copyNumbers, "10"))
                .chr10q(determineCopyNumberQArm(copyNumbers, "10"))
                .chr11p(determineCopyNumberPArm(copyNumbers, "11"))
                .chr11q(determineCopyNumberQArm(copyNumbers, "11"))
                .chr12p(determineCopyNumberPArm(copyNumbers, "12"))
                .chr12q(determineCopyNumberQArm(copyNumbers, "12"))
                .chr13p(determineCopyNumberPArm(copyNumbers, "13"))
                .chr13q(determineCopyNumberQArm(copyNumbers, "13"))
                .chr14p(determineCopyNumberPArm(copyNumbers, "14"))
                .chr14q(determineCopyNumberQArm(copyNumbers, "14"))
                .chr15p(determineCopyNumberPArm(copyNumbers, "15"))
                .chr15q(determineCopyNumberQArm(copyNumbers, "15"))
                .chr16p(determineCopyNumberPArm(copyNumbers, "16"))
                .chr16q(determineCopyNumberQArm(copyNumbers, "16"))
                .chr17p(determineCopyNumberPArm(copyNumbers, "17"))
                .chr17q(determineCopyNumberQArm(copyNumbers, "17"))
                .chr18p(determineCopyNumberPArm(copyNumbers, "18"))
                .chr18q(determineCopyNumberQArm(copyNumbers, "18"))
                .chr19p(determineCopyNumberPArm(copyNumbers, "19"))
                .chr19q(determineCopyNumberQArm(copyNumbers, "19"))
                .chr20p(determineCopyNumberPArm(copyNumbers, "20"))
                .chr20q(determineCopyNumberQArm(copyNumbers, "20"))
                .chr21p(determineCopyNumberPArm(copyNumbers, "21"))
                .chr21q(determineCopyNumberQArm(copyNumbers, "21"))
                .chr22p(determineCopyNumberPArm(copyNumbers, "22"))
                .chr22q(determineCopyNumberQArm(copyNumbers, "22"))
                .chrXp(determineCopyNumberPArm(copyNumbers, "X"))
                .chrXq(determineCopyNumberQArm(copyNumbers, "X"))
                .chrYp(determineCopyNumberPArm(copyNumbers, "Y"))
                .chrYq(determineCopyNumberQArm(copyNumbers, "Y"))
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

    private static double determineCopyNumberPArm(@NotNull List<PurpleCopyNumber> copyNumbers, @NotNull String chromosome) {
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
        return copyNumberArm;
    }

    private static double determineCopyNumberQArm(@NotNull List<PurpleCopyNumber> copyNumbers, @NotNull String chromosome) {
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
        return copyNumberArm;
    }
}