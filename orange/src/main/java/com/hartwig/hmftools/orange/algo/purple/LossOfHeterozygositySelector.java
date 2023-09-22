package com.hartwig.hmftools.orange.algo.purple;

import java.util.List;
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

            if(isRelevantHRD || isRelevantMSI)
            {
                if(hasLOH(geneCopyNumber))
                {
                    suspectGeneCopyNumbersWithLOH.add(geneCopyNumber);
                }
                else if(hasReportedGermlineHetDeletion(geneCopyNumber, allGermlineDeletions))
                {
                    suspectGeneCopyNumbersWithLOH.add(correctForGermlineImpact(geneCopyNumber));
                }
            }
        }
        return suspectGeneCopyNumbersWithLOH;
    }

    private static boolean hasLOH(@NotNull GeneCopyNumber geneCopyNumber)
    {
        return geneCopyNumber.minMinorAlleleCopyNumber() < 0.5 && geneCopyNumber.minCopyNumber() > 0.5;
    }

    private static boolean hasReportedGermlineHetDeletion(@NotNull GeneCopyNumber geneCopyNumber,
            @Nullable List<GermlineDeletion> allGermlineDeletions)
    {
        if(allGermlineDeletions == null)
        {
            return false;
        }

        GermlineDeletion deletion = findByGene(allGermlineDeletions, geneCopyNumber.geneName());
        if(deletion == null)
        {
            return false;
        }

        return deletion.Reported && deletion.TumorStatus == GermlineStatus.HET_DELETION;
    }

    @Nullable
    private static GermlineDeletion findByGene(@NotNull List<GermlineDeletion> deletions, @NotNull String geneNameToFind)
    {
        for(GermlineDeletion deletion : deletions)
        {
            if(deletion.GeneName.equals(geneNameToFind))
            {
                return deletion;
            }
        }
        return null;
    }

    @NotNull
    private static GeneCopyNumber correctForGermlineImpact(@NotNull GeneCopyNumber geneCopyNumber)
    {
        return ImmutableGeneCopyNumber.builder()
                .from(geneCopyNumber)
                .minCopyNumber(geneCopyNumber.minCopyNumber() - geneCopyNumber.minMinorAlleleCopyNumber())
                .minMinorAlleleCopyNumber(0D)
                .build();
    }
}
