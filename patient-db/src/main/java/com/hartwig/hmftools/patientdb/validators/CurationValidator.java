package com.hartwig.hmftools.patientdb.validators;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.ecrf.datamodel.ValidationFinding;
import com.hartwig.hmftools.common.ecrf.formstatus.FormStatusState;
import com.hartwig.hmftools.patientdb.curators.TreatmentCurator;
import com.hartwig.hmftools.patientdb.curators.TumorLocationCurator;

import org.jetbrains.annotations.NotNull;

public final class CurationValidator {

    private CurationValidator() {
    }

    @NotNull
    public static List<ValidationFinding> validateTumorLocationCurator(@NotNull TumorLocationCurator tumorLocationCurator) {
        final List<ValidationFinding> findings = Lists.newArrayList();

        for (String unusedTerm : tumorLocationCurator.unusedSearchTerms()) {
            findings.add(ValidationFinding.of("tumorLocationCuration",
                    "",
                    "", "tumor location search term not used",
                    FormStatusState.UNKNOWN,
                    false,
                    unusedTerm));
        }

        return findings;
    }

    @NotNull
    public static List<ValidationFinding> validateTreatmentCurator(@NotNull TreatmentCurator treatmentCurator) {
        final List<ValidationFinding> findings = Lists.newArrayList();

        for (String unusedTerm : treatmentCurator.unusedSearchTerms()) {
            findings.add(ValidationFinding.of("treatmentCuration",
                    "",
                    "",
                    "treatment search term not used",
                    FormStatusState.UNKNOWN,
                    false,
                    unusedTerm));
        }

        return findings;
    }
}
