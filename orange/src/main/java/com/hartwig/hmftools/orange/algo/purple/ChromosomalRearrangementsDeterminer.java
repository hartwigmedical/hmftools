package com.hartwig.hmftools.orange.algo.purple;

import java.util.List;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeFunctions;
import com.hartwig.hmftools.common.purple.PurpleCopyNumber;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeCoordinates;
import com.hartwig.hmftools.datamodel.orange.OrangeRefGenomeVersion;

import org.jetbrains.annotations.NotNull;

public final class ChromosomalRearrangementsDeterminer
{
    @NotNull
    private final RefGenomeCoordinates refGenomeCoordinates;

    public ChromosomalRearrangementsDeterminer(@NotNull final OrangeRefGenomeVersion refGenomeVersion)
    {
        refGenomeCoordinates = getRefGenomeCoordinates(refGenomeVersion);
    }

    private static final double LOWER_THRESHOLD_TRISOMY = 2.8;
    private static final double UPPER_THRESHOLD_TRISOMY = 3.5;
    private static final double THRESHOLD_DELETION = 0.2;
    private static final double LOWER_LIMIT_PERCENTAGE_TRISOMY = 0.98;
    private static final double UPPER_LIMIT_PERCENTAGE_TRISOMY = 0.90;
    private static final double UPPER_LIMIT_PERCENTAGE_DELETION = 0.98;

    private static final String CHR_1 = "1";
    private static final String CHR_19 = "19";

    public boolean determine1qTrisomy(@NotNull List<PurpleCopyNumber> allSomaticCopyNumbers)
    {
        int centromere1 = refGenomeCoordinates.centromere(CHR_1);
        int length1QArm = refGenomeCoordinates.length(CHR_1) - centromere1;
        int lengthAboveLowerThreshold = 0;
        int lengthBelowUpperThreshold = 0;

        for(PurpleCopyNumber copyNumber : allSomaticCopyNumbers)
        {
            String chromosome = RefGenomeFunctions.stripChrPrefix(copyNumber.chromosome());
            if(chromosome.equals(CHR_1) && overlapsWithQArm(copyNumber, centromere1))
            {
                int length = copyNumber.start() > centromere1 ? copyNumber.length() : copyNumber.end() - centromere1 + 1;

                if(copyNumber.averageTumorCopyNumber() > LOWER_THRESHOLD_TRISOMY)
                {
                    lengthAboveLowerThreshold += length;
                }
                if(copyNumber.averageTumorCopyNumber() < UPPER_THRESHOLD_TRISOMY)
                {
                    lengthBelowUpperThreshold += length;
                }
            }
        }

        return (double) lengthAboveLowerThreshold / length1QArm > LOWER_LIMIT_PERCENTAGE_TRISOMY
                && (double) lengthBelowUpperThreshold / length1QArm > UPPER_LIMIT_PERCENTAGE_TRISOMY;
    }

    public boolean determine1p19qCodeletion(@NotNull List<PurpleCopyNumber> allSomaticCopyNumbers)
    {
        int centromere1 = refGenomeCoordinates.centromere(CHR_1);
        int centromere19 = refGenomeCoordinates.centromere(CHR_19);
        int length19QArm = refGenomeCoordinates.length(CHR_19) - centromere19;
        int lengthDeletion1p = 0;
        int lengthDeletionChr19 = 0;

        for(PurpleCopyNumber copyNumber : allSomaticCopyNumbers)
        {
            String chromosome = RefGenomeFunctions.stripChrPrefix(copyNumber.chromosome());
            if(chromosome.equals(CHR_1) && copyNumber.start() < centromere1 && copyNumber.minorAlleleCopyNumber() < THRESHOLD_DELETION)
            {
                int length = copyNumber.end() < centromere1 ? copyNumber.length() : centromere1 - copyNumber.start() + 1;
                lengthDeletion1p += length;
            }
            if(chromosome.equals(CHR_19) && overlapsWithQArm(copyNumber, centromere19)
                    && copyNumber.minorAlleleCopyNumber() < THRESHOLD_DELETION)
            {
                int length = copyNumber.start() > centromere19 ? copyNumber.length() : copyNumber.end() - centromere19 + 1;
                lengthDeletionChr19 += length;
            }
        }

        return (double) lengthDeletion1p / centromere1 > UPPER_LIMIT_PERCENTAGE_DELETION
                && (double) lengthDeletionChr19 / length19QArm > UPPER_LIMIT_PERCENTAGE_DELETION;
    }

    @NotNull
    private static RefGenomeCoordinates getRefGenomeCoordinates(@NotNull OrangeRefGenomeVersion refGenomeVersion)
    {
        return refGenomeVersion == OrangeRefGenomeVersion.V38 ? RefGenomeCoordinates.COORDS_38 : RefGenomeCoordinates.COORDS_37;
    }

    private static boolean overlapsWithQArm(@NotNull PurpleCopyNumber copyNumber, int centromere)
    {
        return copyNumber.start() > centromere || (copyNumber.start() < centromere && centromere < copyNumber.end());
    }
}
