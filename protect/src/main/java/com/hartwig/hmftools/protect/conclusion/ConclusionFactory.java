package com.hartwig.hmftools.protect.conclusion;

import java.util.List;

import com.hartwig.hmftools.common.purple.gene.GeneCopyNumber;
import com.hartwig.hmftools.common.variant.Variant;
import com.hartwig.hmftools.common.variant.structural.annotation.ReportableGeneFusion;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public class ConclusionFactory {

    private ConclusionFactory() {

    }

    public static String createConclusion(@NotNull String patientPrimaryTumorLocation, double tumorMTL, double tumorMTB, double tumorMSI,
            double chordScore, @NotNull List<ReportableGeneFusion> geneFusions, @NotNull List<GeneCopyNumber> geneCopyNumbers,
            @NotNull List<? extends Variant> passSomaticVariants) {
        String conclusion = Strings.EMPTY;
        String enter = "<enter>";

        String textTumorLocation = createTumorLocationSentense();
        conclusion = textTumorLocation + enter;
        return conclusion;
    }

    private static String createTumorLocationSentense() {
        String tumorLocation = Strings.EMPTY;
        ;
        return tumorLocation + "sample showing";
    }
}
