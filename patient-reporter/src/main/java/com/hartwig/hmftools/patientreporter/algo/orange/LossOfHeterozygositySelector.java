package com.hartwig.hmftools.patientreporter.algo.orange;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.chord.ChordStatus;
import com.hartwig.hmftools.common.purple.GeneCopyNumber;
import com.hartwig.hmftools.common.variant.msi.MicrosatelliteStatus;
import com.hartwig.hmftools.patientreporter.util.Genes;

import org.jetbrains.annotations.NotNull;

public class LossOfHeterozygositySelector {

    private LossOfHeterozygositySelector() {
    }

    @NotNull
    public static List<GeneCopyNumber> selectMSIGenesWithLOH(@NotNull List<GeneCopyNumber> allSomaticGeneCopyNumbers,
            @NotNull MicrosatelliteStatus microsatelliteStatus) {
        List<GeneCopyNumber> suspectGeneCopyNumbersWithLOH = Lists.newArrayList();
        for (GeneCopyNumber geneCopyNumber : allSomaticGeneCopyNumbers) {
            if (hasLOH(geneCopyNumber)) {
                boolean isRelevantMSI = Genes.MSI_GENES.contains(geneCopyNumber.geneName()) && microsatelliteStatus == MicrosatelliteStatus.MSI;

                if (isRelevantMSI) {
                    suspectGeneCopyNumbersWithLOH.add(geneCopyNumber);
                }
            }
        }
        return suspectGeneCopyNumbersWithLOH;
    }

    @NotNull
    public static List<GeneCopyNumber> selectHRDGenesWithLOH(@NotNull List<GeneCopyNumber> allSomaticGeneCopyNumbers,
            @NotNull ChordStatus chordStatus) {
        List<GeneCopyNumber> suspectGeneCopyNumbersWithLOH = Lists.newArrayList();
        for (GeneCopyNumber geneCopyNumber : allSomaticGeneCopyNumbers) {
            if (hasLOH(geneCopyNumber)) {
                boolean isRelevantHRD = Genes.HRD_GENES.contains(geneCopyNumber.geneName()) && chordStatus == ChordStatus.HR_DEFICIENT;

                if (isRelevantHRD) {
                    suspectGeneCopyNumbersWithLOH.add(geneCopyNumber);
                }
            }
        }
        return suspectGeneCopyNumbersWithLOH;
    }

    private static boolean hasLOH(@NotNull GeneCopyNumber geneCopyNumber) {
        return geneCopyNumber.minMinorAlleleCopyNumber() < 0.5 && geneCopyNumber.minCopyNumber() > 0.5;
    }
}
