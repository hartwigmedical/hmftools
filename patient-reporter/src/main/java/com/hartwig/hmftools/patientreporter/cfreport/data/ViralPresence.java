package com.hartwig.hmftools.patientreporter.cfreport.data;

import java.util.List;
import java.util.Set;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.virus.AnnotatedVirus;
import com.hartwig.hmftools.common.virus.VirusConstants;

import org.jetbrains.annotations.NotNull;

public final class ViralPresence {

    private ViralPresence() {
    }

    @NotNull
    public static Set<String> virusInterpretationSummary(@NotNull List<AnnotatedVirus> reportableViruses) {
        Set<String> positiveInterpretations = Sets.newHashSet();
        Set<String> virusInterpretationSummary = Sets.newHashSet();

        for (AnnotatedVirus virus : reportableViruses) {
            String virusInterpretation = virus.interpretation();

            if (virus.interpretation() != null) {
                assert virusInterpretation != null;
                VirusConstants virusConstants = VirusConstants.fromVirusName(virusInterpretation);
                if (virusConstants.reportVirusOnSummary()) {
                    positiveInterpretations.add(virus.interpretation());
                }
            }
        }

        for (String positiveVirus : positiveInterpretations) {
            virusInterpretationSummary.add(positiveVirus + " positive");
        }

        return virusInterpretationSummary;
    }

    @NotNull
    public static String createViralCoverageString(double percentageCovered) {
        return percentageCovered + "%";
    }
}