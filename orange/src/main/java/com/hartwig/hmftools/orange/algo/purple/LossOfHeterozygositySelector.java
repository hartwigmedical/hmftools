package com.hartwig.hmftools.orange.algo.purple;

import java.util.List;
import java.util.OptionalDouble;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.chord.ChordStatus;
import com.hartwig.hmftools.common.purple.GeneCopyNumber;
import com.hartwig.hmftools.common.purple.GermlineDeletion;
import com.hartwig.hmftools.common.purple.GermlineStatus;
import com.hartwig.hmftools.common.purple.ImmutableGeneCopyNumber;
import com.hartwig.hmftools.common.variant.msi.MicrosatelliteStatus;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

final class LossOfHeterozygositySelector
{

    static final Set<String> HRD_GENES = Sets.newHashSet();
    static final Set<String> MSI_GENES = Sets.newHashSet();

    static
    {
        HRD_GENES.add("BRCA1");
        HRD_GENES.add("BRCA2");
        HRD_GENES.add("PALB2");
        HRD_GENES.add("RAD51C");
        HRD_GENES.add("RAD51B");

        MSI_GENES.add("MSH6");
        MSI_GENES.add("MSH2");
        MSI_GENES.add("MLH1");
        MSI_GENES.add("PMS2");
        MSI_GENES.add("EPCAM");
    }

    private LossOfHeterozygositySelector()
    {
    }

    @NotNull
    public static List<GeneCopyNumber> selectHRDOrMSIGenesWithLOH(@NotNull List<GeneCopyNumber> allSomaticGeneCopyNumbers,
            @Nullable List<GermlineDeletion> allGermlineDeletions, @NotNull MicrosatelliteStatus microsatelliteStatus,
            @Nullable ChordStatus chordStatus)
    {
        List<GeneCopyNumber> suspectGeneCopyNumbersWithLOH = Lists.newArrayList();
        for(GeneCopyNumber geneCopyNumber : allSomaticGeneCopyNumbers)
        {
            boolean isRelevantHRD = HRD_GENES.contains(geneCopyNumber.geneName()) && chordStatus == ChordStatus.HR_DEFICIENT;
            boolean isRelevantMSI = MSI_GENES.contains(geneCopyNumber.geneName()) && microsatelliteStatus == MicrosatelliteStatus.MSI;
            boolean hasReportedGermlineHomozygousDeletion = hasReportedGermlineHomDeletion(geneCopyNumber, allGermlineDeletions);

            if((isRelevantHRD || isRelevantMSI) && !hasReportedGermlineHomozygousDeletion)
            {
                if(hasLOH(geneCopyNumber))
                {
                    suspectGeneCopyNumbersWithLOH.add(geneCopyNumber);
                }
                else if(notFullyLostInTumor(geneCopyNumber) && hasReportedGermlineHetDeletion(geneCopyNumber, allGermlineDeletions))
                {
                    final GeneCopyNumber correctedGeneCopyNumber = correctForGermlineImpact(geneCopyNumber, allGermlineDeletions);
                    suspectGeneCopyNumbersWithLOH.add(correctedGeneCopyNumber);
                }
            }
        }
        return suspectGeneCopyNumbersWithLOH;
    }

    private static boolean hasLOH(@NotNull GeneCopyNumber geneCopyNumber)
    {
        return geneCopyNumber.minMinorAlleleCopyNumber() < 0.5 && notFullyLostInTumor(geneCopyNumber);
    }

    private static boolean notFullyLostInTumor(@NotNull GeneCopyNumber geneCopyNumber)
    {
        return geneCopyNumber.minCopyNumber() > 0.5;
    }

    private static boolean hasReportedGermlineHetDeletion(@NotNull GeneCopyNumber geneCopyNumber,
            @Nullable List<GermlineDeletion> allGermlineDeletions)
    {
        if(allGermlineDeletions == null)
        {
            return false;
        }

        return allGermlineDeletions.stream()
                .anyMatch(d -> isMatchingReportedGermlineDeletion(d, geneCopyNumber.geneName(), GermlineStatus.HET_DELETION));
    }

    private static boolean hasReportedGermlineHomDeletion(@NotNull GeneCopyNumber geneCopyNumber,
            @Nullable List<GermlineDeletion> allGermlineDeletions)
    {
        if(allGermlineDeletions == null)
        {
            return false;
        }

        return allGermlineDeletions.stream()
                .anyMatch(d -> isMatchingReportedGermlineDeletion(d, geneCopyNumber.geneName(), GermlineStatus.HOM_DELETION));
    }

    @NotNull
    private static GeneCopyNumber correctForGermlineImpact(@NotNull GeneCopyNumber geneCopyNumber,
            @NotNull List<GermlineDeletion> allGermlineDeletions)
    {
        double adjustedMinCopyNumber = adjustMinCopyNumberForGermlineImpact(geneCopyNumber, allGermlineDeletions);
        return ImmutableGeneCopyNumber.builder()
                .from(geneCopyNumber)
                .minCopyNumber(adjustedMinCopyNumber)
                .minMinorAlleleCopyNumber(0D)
                .build();
    }

    private static double adjustMinCopyNumberForGermlineImpact(GeneCopyNumber geneCopyNumber,
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

    private static boolean isMatchingReportedGermlineDeletion(GermlineDeletion deletion, String geneName, GermlineStatus tumorStatus)
    {
        return deletion.GeneName.equals(geneName) && deletion.Reported && deletion.TumorStatus == tumorStatus;
    }
}
