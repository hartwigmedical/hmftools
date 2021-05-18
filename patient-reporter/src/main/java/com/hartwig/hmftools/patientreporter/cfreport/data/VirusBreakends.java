package com.hartwig.hmftools.patientreporter.cfreport.data;

import java.util.List;
import java.util.Set;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.patientreporter.virusbreakend.ReportableVirusBreakend;

import org.jetbrains.annotations.NotNull;

public final class VirusBreakends {

    private VirusBreakends() {
    }

    @NotNull
    public static String virusInterpretationSummary(@NotNull List<ReportableVirusBreakend> reportableVirusBreakends) {
        Set<String> positiveInterpretations = Sets.newHashSet();
        for (ReportableVirusBreakend virusBreakend : reportableVirusBreakends) {
            if (virusBreakend.interpretation() != null) {
                positiveInterpretations.add(virusBreakend.interpretation());
            }
        }

        // TODO Also add "negative" for viruses that are not found.
        StringBuilder virusInterpretationSummary = new StringBuilder(",");
        for (String positiveVirus : positiveInterpretations) {
            virusInterpretationSummary.append(positiveVirus + " positive");
        }
        return virusInterpretationSummary.toString();
    }
}