package com.hartwig.hmftools.orange.algo.selection;

import java.util.List;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.chord.ChordStatus;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumber;
import com.hartwig.hmftools.common.variant.msi.MicrosatelliteStatus;

import org.jetbrains.annotations.NotNull;

public final class LossOfHeterozygositySelector {

    private static final Set<String> HRD_GENES = Sets.newHashSet();
    private static final Set<String> MSI_GENES = Sets.newHashSet();

    static {
        HRD_GENES.add("BRCA1");
        HRD_GENES.add("BRCA2");
        HRD_GENES.add("PALB2");
        HRD_GENES.add("RAD51C");

        MSI_GENES.add("MSH6");
        MSI_GENES.add("MSH2");
        MSI_GENES.add("MLH1");
        MSI_GENES.add("PMS2");
        MSI_GENES.add("EPCAM");
    }

    private LossOfHeterozygositySelector() {
    }

    @NotNull
    public static List<GeneCopyNumber> selectHRDOrMSIGenesWithLOH(@NotNull List<GeneCopyNumber> allGeneCopyNumbers,
            @NotNull MicrosatelliteStatus microsatelliteStatus, @NotNull ChordStatus chordStatus) {
        List<GeneCopyNumber> reportable = Lists.newArrayList();
        for (GeneCopyNumber geneCopyNumber : allGeneCopyNumbers) {
            if (geneCopyNumber.minMinorAlleleCopyNumber() < 0.5 && geneCopyNumber.minCopyNumber() > 0.5) {
                boolean isRelevantHRD = HRD_GENES.contains(geneCopyNumber.geneName()) && chordStatus == ChordStatus.HR_DEFICIENT;
                boolean isRelevantMSI = MSI_GENES.contains(geneCopyNumber.geneName()) && microsatelliteStatus == MicrosatelliteStatus.MSI;

                if (isRelevantHRD || isRelevantMSI) {
                    reportable.add(geneCopyNumber);
                }
            }
        }
        return reportable;
    }
}
