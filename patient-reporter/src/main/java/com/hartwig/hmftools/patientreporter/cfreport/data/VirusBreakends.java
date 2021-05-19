package com.hartwig.hmftools.patientreporter.cfreport.data;

import java.util.List;
import java.util.Set;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.patientreporter.virusbreakend.ReportableVirusBreakend;

import org.jetbrains.annotations.NotNull;

public final class VirusBreakends {

    private static final Set<String> INTERPRETATIONS_TO_EVALUATE = Sets.newHashSet("HPV", "EBV", "MCV");

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

        Set<String> negativeInterpretations = Sets.newHashSet();
        for (String possibleViruses : INTERPRETATIONS_TO_EVALUATE) {
            if (!positiveInterpretations.contains(possibleViruses)) {
                negativeInterpretations.add(possibleViruses);
            }
        }

        StringBuilder virusInterpretationSummary = new StringBuilder();
        for (String positiveVirus : positiveInterpretations) {
            virusInterpretationSummary.append(positiveVirus + " positive");
            virusInterpretationSummary.append(", ");

        }

        for (String negativeVirus : negativeInterpretations) {
            virusInterpretationSummary.append(negativeVirus + " negative");
            virusInterpretationSummary.append(", ");
        }

        String virusInterpretationSummaryString = virusInterpretationSummary.toString();
        //TODO remove last comma
        return virusInterpretationSummaryString;
    }
}