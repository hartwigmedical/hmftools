package com.hartwig.hmftools.orange.algo.purple;

import java.util.List;

import com.hartwig.hmftools.common.purple.PurpleCopyNumber;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeCoordinates;
import com.hartwig.hmftools.datamodel.orange.OrangeRefGenomeVersion;

import org.jetbrains.annotations.NotNull;

final class ChromosomalRearrangementsDeterminer
{
    private final RefGenomeCoordinates refGenomeCoordinates;

    public ChromosomalRearrangementsDeterminer(@NotNull final OrangeRefGenomeVersion refGenomeVersion)
    {
        refGenomeCoordinates = getRefGenomeCoordinates(refGenomeVersion);
    }

    public static final double LOWER_THRESHOLD_TRISOMY = 2.8;
    public static final double UPPER_THRESHOLD_TRISOMY = 3.5;
    public static final double THRESHOLD_DELETION = 0.2;
    public static final double PCT_98 = 0.98;
    public static final double PCT_90 = 0.90;

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
            if(CHR_1.equals(copyNumber.chromosome()) && isQArm(copyNumber, centromere1))
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

        return (double) lengthAboveLowerThreshold / length1QArm > PCT_98 && (double) lengthBelowUpperThreshold / length1QArm > PCT_90;
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
            if(CHR_1.equals(copyNumber.chromosome()) && copyNumber.start() < centromere1
                    && copyNumber.minorAlleleCopyNumber() < THRESHOLD_DELETION)
            {
                int length = copyNumber.end() < centromere1 ? copyNumber.length() : centromere1 - copyNumber.start() + 1;
                lengthDeletion1p += length;
            }
            if(CHR_19.equals(copyNumber.chromosome()) && isQArm(copyNumber, centromere19)
                    && copyNumber.minorAlleleCopyNumber() < THRESHOLD_DELETION)
            {
                int length = copyNumber.start() > centromere19 ? copyNumber.length() : copyNumber.end() - centromere19 + 1;
                lengthDeletionChr19 += length;
            }
        }

        return (double) lengthDeletion1p / centromere1 > PCT_98 && (double) lengthDeletionChr19 / length19QArm > PCT_98;
    }

    @NotNull
    private static RefGenomeCoordinates getRefGenomeCoordinates(@NotNull OrangeRefGenomeVersion refGenomeVersion)
    {
        return refGenomeVersion == OrangeRefGenomeVersion.V38 ? RefGenomeCoordinates.COORDS_38 : RefGenomeCoordinates.COORDS_37;
    }

    private static boolean isQArm(@NotNull PurpleCopyNumber copyNumber, int centromere)
    {
        return copyNumber.start() > centromere || (copyNumber.start() < centromere && centromere < copyNumber.end());
    }
}
