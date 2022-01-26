package com.hartwig.hmftools.patientreporter.cfreport.data;

import java.util.List;
import java.util.Set;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.virus.AnnotatedVirus;
import com.hartwig.hmftools.common.virus.VirusLikelihoodType;

import org.jetbrains.annotations.NotNull;

public final class ViralPresence {

    private ViralPresence() {
    }

    @NotNull
    public static Set<String> virusInterpretationSummary(@NotNull List<AnnotatedVirus> reportableViruses) {
        Set<String> positiveHighDriverInterpretations = Sets.newHashSet();
        Set<String> virusInterpretationSummary = Sets.newHashSet();

        for (AnnotatedVirus virus : reportableViruses) {
            if (virus.virusDriverLikelihoodType() == VirusLikelihoodType.HIGH) {
                String virusInterpretation = virus.interpretation();
                if (virus.interpretation() != null) {
                    assert virusInterpretation != null;
                    positiveHighDriverInterpretations.add(virus.interpretation());
                }
            }
        }

        for (String positiveVirus : positiveHighDriverInterpretations) {
            virusInterpretationSummary.add(positiveVirus + " positive");
        }

        return virusInterpretationSummary;
    }

    @NotNull
    public static String createViralCoverageString(double percentageCovered) {
        return Math.round(percentageCovered) + "%";
    }

    @NotNull
    public static String createIntegrationSiteString(Integer integrations) {
        return integrations == 0 ? "Detected without integration sites" : Integer.toString(integrations);
    }
}