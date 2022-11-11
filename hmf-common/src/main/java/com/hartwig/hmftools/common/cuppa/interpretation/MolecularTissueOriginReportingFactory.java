package com.hartwig.hmftools.common.cuppa.interpretation;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public class MolecularTissueOriginReportingFactory {

    private MolecularTissueOriginReportingFactory() {
    }

    private static final String RESULTS_INCONCLUSIVE = "results inconclusive";

    @NotNull
    public static String curatedTumorLocation(@NotNull String cancerType) {
        if (cancerType.equals("Uterus: Endometrium")) {
            cancerType = "Endometrium";
        } else if (cancerType.equals("Colorectum/Appendix/SmallIntestine")) {
            cancerType = "Lower GI tract";
        }
        return cancerType;
    }

    @NotNull
    public static String interpretTumorLocation(double likelihood, @NotNull String cancerType) {
        // our cut-off is 80% likelihood. When this is below 80% then the results is inconclusive
        String interpretCancerType = Strings.EMPTY;
        if (likelihood <= 0.8) {
            interpretCancerType = RESULTS_INCONCLUSIVE;
        } else {
            interpretCancerType = cancerType;
        }
        return interpretCancerType;
    }

    @NotNull
    public static MolecularTissueOriginReporting createMolecularTissueOriginReportingData(@NotNull CuppaPrediction bestPrediction) {
        double likelihood = bestPrediction.likelihood();
        String cancerType = curatedTumorLocation(bestPrediction.cancerType());
        String interpretCancerType = interpretTumorLocation(likelihood, cancerType);
        Double interpretLikelihood = likelihood >= 0.8 ? likelihood : null;

        return ImmutableMolecularTissueOriginReporting.builder()
                .bestCancerType(cancerType)
                .bestLikelihood(likelihood)
                .interpretCancerType(interpretCancerType)
                .interpretLikelihood(interpretLikelihood)
                .build();
    }
}