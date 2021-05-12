package com.hartwig.hmftools.patientreporter.cnchromosome;

import java.io.IOException;
import java.util.List;

import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeCoordinates;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumberFile;
import com.hartwig.hmftools.common.purple.segment.ChromosomeArm;

import org.jetbrains.annotations.NotNull;

public class CnPerChromosomeFactory {

    public CnPerChromosomeFactory() {

    }

    public static RefGenomeCoordinates refGenomeLengths() {
        return RefGenomeVersion.V38 == RefGenomeVersion.V38 ? RefGenomeCoordinates.COORDS_38 : RefGenomeCoordinates.COORDS_37;
    }

    public static int getChromosomalArmLength(final String chromosome, final ChromosomeArm armType) {
        final RefGenomeCoordinates refGenome = refGenomeLengths();
        final HumanChromosome chr = HumanChromosome.fromString(chromosome);

        final Long centromerePos = refGenome.centromeres().get(chr);

        if (centromerePos == null) {
            return 0;
        }

        if (armType == ChromosomeArm.P_ARM) {
            return centromerePos.intValue();
        }

        int chrLength = refGenome.lengths().get(chr).intValue();

        return chrLength - centromerePos.intValue();
    }

    public static double determineCopyNumberArm(@NotNull List<PurpleCopyNumber> copyNumbers, @NotNull String chromosome,
            @NotNull ChromosomeArm chromosomeArm) {
        int chromosomeLength = getChromosomalArmLength(chromosome, chromosomeArm);

        double copyNumber = 0;
        double copyNumberArm = 0;
        for (PurpleCopyNumber purpleCopyNumber : copyNumbers) {
            if (purpleCopyNumber.chromosome().equals(chromosome) && purpleCopyNumber.end() < chromosomeLength) {
                copyNumber += purpleCopyNumber.averageObservedBAF();
                long totalLengthSegemnt = (purpleCopyNumber.end() - purpleCopyNumber.start()) - 1;
                copyNumber += copyNumber * totalLengthSegemnt / chromosomeLength;
            }

        }
        return copyNumberArm;
    }

    public static CnPerChromosome extractCnPerChromsomsoemArm(@NotNull String purpleCnvSomaticTsv) throws IOException {

        List<PurpleCopyNumber> copyNumbers = PurpleCopyNumberFile.read(purpleCnvSomaticTsv);

        return ImmutableCnPerChromosome.builder()
                .chr1p(determineCopyNumberArm(copyNumbers, "1", ChromosomeArm.P_ARM))
                .chr1q(determineCopyNumberArm(copyNumbers, "1", ChromosomeArm.Q_ARM))
                .chr2p(determineCopyNumberArm(copyNumbers, "2", ChromosomeArm.P_ARM))
                .chr2q(determineCopyNumberArm(copyNumbers, "2", ChromosomeArm.Q_ARM))
                .chr3p(determineCopyNumberArm(copyNumbers, "3", ChromosomeArm.P_ARM))
                .chr3q(determineCopyNumberArm(copyNumbers, "3", ChromosomeArm.Q_ARM))
                .chr4p(determineCopyNumberArm(copyNumbers, "4", ChromosomeArm.P_ARM))
                .chr4q(determineCopyNumberArm(copyNumbers, "4", ChromosomeArm.Q_ARM))
                .chr5p(determineCopyNumberArm(copyNumbers, "5", ChromosomeArm.P_ARM))
                .chr5q(determineCopyNumberArm(copyNumbers, "5", ChromosomeArm.Q_ARM))
                .chr6p(determineCopyNumberArm(copyNumbers, "6", ChromosomeArm.P_ARM))
                .chr6q(determineCopyNumberArm(copyNumbers, "6", ChromosomeArm.Q_ARM))
                .chr7p(determineCopyNumberArm(copyNumbers, "7", ChromosomeArm.P_ARM))
                .chr7q(determineCopyNumberArm(copyNumbers, "7", ChromosomeArm.Q_ARM))
                .chr8p(determineCopyNumberArm(copyNumbers, "8", ChromosomeArm.P_ARM))
                .chr8q(determineCopyNumberArm(copyNumbers, "8", ChromosomeArm.Q_ARM))
                .chr9p(determineCopyNumberArm(copyNumbers, "9", ChromosomeArm.P_ARM))
                .chr9q(determineCopyNumberArm(copyNumbers, "9", ChromosomeArm.Q_ARM))
                .chr10p(determineCopyNumberArm(copyNumbers, "10", ChromosomeArm.P_ARM))
                .chr10q(determineCopyNumberArm(copyNumbers, "10", ChromosomeArm.Q_ARM))
                .chr11p(determineCopyNumberArm(copyNumbers, "11", ChromosomeArm.P_ARM))
                .chr11q(determineCopyNumberArm(copyNumbers, "11", ChromosomeArm.Q_ARM))
                .chr12p(determineCopyNumberArm(copyNumbers, "12", ChromosomeArm.P_ARM))
                .chr12q(determineCopyNumberArm(copyNumbers, "12", ChromosomeArm.Q_ARM))
                .chr13p(determineCopyNumberArm(copyNumbers, "13", ChromosomeArm.P_ARM))
                .chr13q(determineCopyNumberArm(copyNumbers, "13", ChromosomeArm.Q_ARM))
                .chr14p(determineCopyNumberArm(copyNumbers, "14", ChromosomeArm.P_ARM))
                .chr14q(determineCopyNumberArm(copyNumbers, "14", ChromosomeArm.Q_ARM))
                .chr15p(determineCopyNumberArm(copyNumbers, "15", ChromosomeArm.P_ARM))
                .chr15q(determineCopyNumberArm(copyNumbers, "15", ChromosomeArm.Q_ARM))
                .chr16p(determineCopyNumberArm(copyNumbers, "16", ChromosomeArm.P_ARM))
                .chr16q(determineCopyNumberArm(copyNumbers, "16", ChromosomeArm.Q_ARM))
                .chr17p(determineCopyNumberArm(copyNumbers, "17", ChromosomeArm.P_ARM))
                .chr17q(determineCopyNumberArm(copyNumbers, "17", ChromosomeArm.Q_ARM))
                .chr18p(determineCopyNumberArm(copyNumbers, "18", ChromosomeArm.P_ARM))
                .chr18q(determineCopyNumberArm(copyNumbers, "18", ChromosomeArm.Q_ARM))
                .chr19p(determineCopyNumberArm(copyNumbers, "19", ChromosomeArm.P_ARM))
                .chr19q(determineCopyNumberArm(copyNumbers, "19", ChromosomeArm.Q_ARM))
                .chr20p(determineCopyNumberArm(copyNumbers, "20", ChromosomeArm.P_ARM))
                .chr20q(determineCopyNumberArm(copyNumbers, "20", ChromosomeArm.Q_ARM))
                .chr21p(determineCopyNumberArm(copyNumbers, "21", ChromosomeArm.P_ARM))
                .chr21q(determineCopyNumberArm(copyNumbers, "21", ChromosomeArm.Q_ARM))
                .chr22p(determineCopyNumberArm(copyNumbers, "22", ChromosomeArm.P_ARM))
                .chr22q(determineCopyNumberArm(copyNumbers, "22", ChromosomeArm.Q_ARM))
                .chrXp(determineCopyNumberArm(copyNumbers, "X", ChromosomeArm.P_ARM))
                .chrXq(determineCopyNumberArm(copyNumbers, "X", ChromosomeArm.Q_ARM))
                .chrYp(determineCopyNumberArm(copyNumbers, "Y", ChromosomeArm.P_ARM))
                .chrYq(determineCopyNumberArm(copyNumbers, "Y", ChromosomeArm.Q_ARM))
                .build();
    }
}
