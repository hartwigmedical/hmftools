package com.hartwig.hmftools.patientreporter.cfreport.data;

import java.util.Collection;
import java.util.List;
import java.util.Set;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.patientreporter.virusbreakend.ReportableVirusBreakend;

import org.jetbrains.annotations.NotNull;

public final class VirusBreakends {

    private VirusBreakends() {
    }

    @NotNull
    public static String virusInterpretationSummary(@NotNull List<ReportableVirusBreakend> reportableVirusBreakends,
            @NotNull Collection<String> interpretationVirus) {
        Set<String> positiveInterpretations = Sets.newHashSet();
        Set<String> negativeInterpretations = Sets.newHashSet();
        Set<String> virusInterpretationUnique = Sets.newHashSet();
        StringBuilder virusInterpretationSummary = new StringBuilder(",");

        for (ReportableVirusBreakend virusBreakend : reportableVirusBreakends) {
            if (virusBreakend.interpretation() != null) {
                positiveInterpretations.add(virusBreakend.interpretation());
            }
        }

        for (String possibleVirussen: interpretationVirus) {
            if (!positiveInterpretations.contains(possibleVirussen)) {
                negativeInterpretations.add(possibleVirussen);
            }
        }

        for (String positiveVirus : positiveInterpretations) {
            virusInterpretationUnique.add(positiveVirus + " positive");
        }

        for (String negativeVirus : negativeInterpretations) {
            virusInterpretationUnique.add(negativeVirus + " negative");
        }


        for (String virusInterpretation : virusInterpretationUnique) {
            virusInterpretationSummary.append(virusInterpretation);
        }

        return virusInterpretationSummary.toString();
    }
}