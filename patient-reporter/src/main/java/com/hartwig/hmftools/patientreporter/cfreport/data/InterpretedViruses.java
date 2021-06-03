package com.hartwig.hmftools.patientreporter.cfreport.data;

import java.util.List;
import java.util.Set;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.virus.InterpretedVirus;
import com.hartwig.hmftools.common.virus.VirusInterpretation;

import org.jetbrains.annotations.NotNull;

public final class InterpretedViruses {

    private InterpretedViruses() {
    }

    @NotNull
    public static Set<String> virusInterpretationSummary(@NotNull List<InterpretedVirus> reportableViruses) {
        Set<VirusInterpretation> positiveInterpretations = Sets.newHashSet();
        for (InterpretedVirus virus : reportableViruses) {
            if (virus.interpretation() != null) {
                positiveInterpretations.add(virus.interpretation());
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