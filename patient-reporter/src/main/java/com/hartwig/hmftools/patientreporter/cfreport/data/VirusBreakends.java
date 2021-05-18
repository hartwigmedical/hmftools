package com.hartwig.hmftools.patientreporter.cfreport.data;

import java.util.List;
import java.util.Set;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.patientreporter.virusbreakend.ReportableVirusBreakend;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class VirusBreakends {

    private VirusBreakends() {

    }

    @NotNull
    public static String virusInterpretationSummary(@NotNull List<ReportableVirusBreakend> reportableVirusBreakends) {

        StringBuilder virusInterpretationSummary = new StringBuilder(",");

        for (ReportableVirusBreakend virusBreakend : reportableVirusBreakends) {
            if (determineInterpretation(virusBreakend.interpretations())) {
                virusInterpretationSummary.append(virusBreakend.interpretations() + " positive");
            } else {
                virusInterpretationSummary.append(virusBreakend.interpretations() + " negative");
            }
        }
        return virusInterpretationSummary.toString();

    }

    private static boolean determineInterpretation(@Nullable String interpretation) {
        Set<String> expectedInterpretations = Sets.newHashSet("EBV", "HPV", "MCV");
        if (interpretation != null && expectedInterpretations.contains(interpretation)) {
            return true;
        }
        return false;
    }
}