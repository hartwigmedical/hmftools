package com.hartwig.hmftools.patientreporter.cfreport.data;

import java.util.List;
import java.util.Set;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.virus.AnnotatedVirus;

import org.jetbrains.annotations.NotNull;

public final class ViralPresence {

    private ViralPresence() {
    }

    @NotNull
    public static Set<String> virusInterpretationSummary(@NotNull List<AnnotatedVirus> reportableViruses) {
        Set<String> positiveInterpretations = Sets.newHashSet();
        Set<String> negativeInterpretations = Sets.newHashSet();
        Set<String> virusInterpretationSummary = Sets.newHashSet();

        for (AnnotatedVirus virus : reportableViruses) {
            if (virus.reportedSummary()) {
                positiveInterpretations.add(virus.interpretation());
            } else if (!virus.reportedSummary()) {
                negativeInterpretations.add(virus.interpretation());
            }
        }

        for (String positiveVirus : positiveInterpretations) {
            virusInterpretationSummary.add(positiveVirus + " positive");
        }

        for (String negativeVirus : negativeInterpretations) {
            virusInterpretationSummary.add(negativeVirus + " negative");
        }

        return virusInterpretationSummary;
    }
}