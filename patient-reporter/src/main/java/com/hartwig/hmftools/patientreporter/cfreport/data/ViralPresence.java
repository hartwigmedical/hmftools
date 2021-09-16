package com.hartwig.hmftools.patientreporter.cfreport.data;

import java.util.List;
import java.util.Set;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.virus.AnnotatedVirusV1;
import com.hartwig.hmftools.common.virus.VirusConstants;

import org.jetbrains.annotations.NotNull;

public final class ViralPresence {

    private ViralPresence() {
    }

    @NotNull
    public static Set<String> virusInterpretationSummary(@NotNull List<AnnotatedVirusV1> reportableViruses) {
        Set<String> positiveInterpretations = Sets.newHashSet();
        Set<String> negativeInterpretations = Sets.newHashSet();
        Set<String> virusInterpretationSummary = Sets.newHashSet();

        for (AnnotatedVirusV1 virus : reportableViruses) {
            VirusConstants virusConstants = VirusConstants.fromVirusName(virus.interpretation());
            if (virusConstants.reportVirusOnSummary()) {
                positiveInterpretations.add(virus.interpretation());
            } else if (!virusConstants.reportVirusOnSummary()) {
                negativeInterpretations.add(virus.interpretation());
            }
        }

        for (String virusSummary: VirusConstants.allViruses()) {
            if (!positiveInterpretations.contains(virusSummary)) {
                negativeInterpretations.add(virusSummary);
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