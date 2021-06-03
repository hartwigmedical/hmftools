package com.hartwig.hmftools.patientreporter.cfreport.data;

import java.util.List;
import java.util.Set;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.virus.VirusInterpretation;
import com.hartwig.hmftools.patientreporter.virusbreakend.ReportableVirusBreakend;

import org.jetbrains.annotations.NotNull;

public final class VirusBreakends {

    private VirusBreakends() {
    }

    @NotNull
    public static Set<String> virusInterpretationSummary(@NotNull List<ReportableVirusBreakend> reportableVirusBreakends) {
        Set<VirusInterpretation> positiveInterpretations = Sets.newHashSet();
        for (ReportableVirusBreakend virusBreakend : reportableVirusBreakends) {
            if (virusBreakend.interpretation() != null) {
                positiveInterpretations.add(virusBreakend.interpretation());
            }
        }

        Set<VirusInterpretation> negativeInterpretations = Sets.newHashSet();
        for (VirusInterpretation possibleViruses : VirusInterpretation.values()) {
            if (!positiveInterpretations.contains(possibleViruses)) {
                negativeInterpretations.add(possibleViruses);
            }
        }

        Set<String> virusInterpretationSummary = Sets.newHashSet();
        for (VirusInterpretation positiveVirus : positiveInterpretations) {
            virusInterpretationSummary.add(positiveVirus.toString() + " positive");
        }

        for (VirusInterpretation negativeVirus : negativeInterpretations) {
            virusInterpretationSummary.add(negativeVirus.toString() + " negative");
        }

        return virusInterpretationSummary;
    }
}