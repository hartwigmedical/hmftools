package com.hartwig.hmftools.orange.algo.purple;

import java.util.List;
import java.util.OptionalDouble;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.chord.ChordStatus;
import com.hartwig.hmftools.common.purple.GeneCopyNumber;
import com.hartwig.hmftools.common.purple.GermlineDeletion;
import com.hartwig.hmftools.common.purple.GermlineStatus;
import com.hartwig.hmftools.common.purple.ImmutableGeneCopyNumber;
import com.hartwig.hmftools.common.variant.msi.MicrosatelliteStatus;
import com.hartwig.hmftools.datamodel.purple.Genes;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

final class LossOfHeterozygositySelector
{

    @NotNull
    public static List<GeneCopyNumber> selectHRDOrMSIGenesWithLOH(@NotNull List<GeneCopyNumber> allSomaticGeneCopyNumbers,
            @Nullable List<GermlineDeletion> allGermlineDeletions, @NotNull MicrosatelliteStatus microsatelliteStatus,
            @Nullable ChordStatus chordStatus)
    {
        List<GeneCopyNumber> suspectGeneCopyNumbersWithLOH = Lists.newArrayList();
        for(GeneCopyNumber geneCopyNumber : allSomaticGeneCopyNumbers)
        {
            boolean isRelevantHRD = Genes.HRD_GENES.contains(geneCopyNumber.geneName()) && chordStatus == ChordStatus.HR_DEFICIENT;
            boolean isRelevantMSI = Genes.MSI_GENES.contains(geneCopyNumber.geneName()) && microsatelliteStatus == MicrosatelliteStatus.MSI;
            boolean hasReportedGermlineDel = hasReportedGermlineDeletionWithTumorStatus(geneCopyNumber.geneName(),
                    GermlineStatus.HOM_DELETION, allGermlineDeletions);
            boolean fullyDeletedInTumor = geneCopyNumber.minCopyNumber() < 0.5;
            if((isRelevantHRD || isRelevantMSI) && !hasReportedGermlineDel && !fullyDeletedInTumor)
            {
                if(geneCopyNumber.minMinorAlleleCopyNumber() < 0.5)
                {
                    suspectGeneCopyNumbersWithLOH.add(geneCopyNumber);
                }
                else if(hasReportedGermlineDeletionWithTumorStatus(geneCopyNumber.geneName(), GermlineStatus.HET_DELETION,
                        allGermlineDeletions))
                {
                    suspectGeneCopyNumbersWithLOH.add(correctForGermlineImpact(geneCopyNumber, allGermlineDeletions));
                }
            }
        }
        return suspectGeneCopyNumbersWithLOH;
    }

    private static boolean hasReportedGermlineDeletionWithTumorStatus(@NotNull String geneName, @NotNull GermlineStatus tumorStatus,
            @Nullable List<GermlineDeletion> allGermlineDeletions)
    {
        if(allGermlineDeletions == null)
        {
            return false;
        }

        return allGermlineDeletions.stream().anyMatch(d -> isMatchingReportedGermlineDeletion(d, geneName, tumorStatus));
    }

    @NotNull
    private static GeneCopyNumber correctForGermlineImpact(@NotNull GeneCopyNumber geneCopyNumber,
            @NotNull List<GermlineDeletion> allGermlineDeletions)
    {
        return ImmutableGeneCopyNumber.builder()
                .from(geneCopyNumber)
                .minCopyNumber(adjustMinCopyNumberForGermlineImpact(geneCopyNumber, allGermlineDeletions))
                .minMinorAlleleCopyNumber(0D)
                .maxCopyNumber(adjustMaxCopyNumberForGermlineImpact(geneCopyNumber, allGermlineDeletions))
                .build();
    }

    private static double adjustMinCopyNumberForGermlineImpact(@NotNull GeneCopyNumber geneCopyNumber,
            @NotNull List<GermlineDeletion> allGermlineDeletions)
    {
        OptionalDouble minimumTumorCopyNumber = allGermlineDeletions.stream()
                .filter(d1 -> isMatchingReportedGermlineDeletion(d1, geneCopyNumber.geneName(), GermlineStatus.HET_DELETION))
                .mapToDouble(d -> d.TumorCopyNumber)
                .min();
        if(minimumTumorCopyNumber.isPresent())
        {
            return Math.min(minimumTumorCopyNumber.getAsDouble(), geneCopyNumber.minCopyNumber());
        }
        else
        {
            return geneCopyNumber.minCopyNumber();
        }
    }

    private static double adjustMaxCopyNumberForGermlineImpact(@NotNull GeneCopyNumber geneCopyNumber,
            @NotNull List<GermlineDeletion> allGermlineDeletions)
    {
        OptionalDouble maximumTumorCopyNumber = allGermlineDeletions.stream()
                .filter(d1 -> isMatchingReportedGermlineDeletion(d1, geneCopyNumber.geneName(), GermlineStatus.HET_DELETION))
                .mapToDouble(d -> d.TumorCopyNumber)
                .max();
        if(maximumTumorCopyNumber.isPresent())
        {
            return Math.max(maximumTumorCopyNumber.getAsDouble(), geneCopyNumber.maxCopyNumber());
        }
        else
        {
            return geneCopyNumber.maxCopyNumber();
        }
    }

    private static boolean isMatchingReportedGermlineDeletion(GermlineDeletion deletion, String geneName, GermlineStatus tumorStatus)
    {
        return deletion.GeneName.equals(geneName) && deletion.Reported && deletion.TumorStatus == tumorStatus;
    }
}
