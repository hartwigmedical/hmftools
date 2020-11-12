package com.hartwig.hmftools.patientdb.validators;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.ecrf.datamodel.ValidationFinding;
import com.hartwig.hmftools.common.ecrf.formstatus.FormStatus;
import com.hartwig.hmftools.patientdb.curators.CleanableCurator;
import com.hartwig.hmftools.patientdb.curators.PrimaryTumorCurator;
import com.hartwig.hmftools.patientdb.curators.TreatmentCurator;

import org.jetbrains.annotations.NotNull;

public final class CurationValidator {

    private CurationValidator() {
    }

    @NotNull
    public static List<ValidationFinding> validatePrimaryTumorCurator(@NotNull PrimaryTumorCurator primaryTumorCurator) {
        return validateCurator(primaryTumorCurator, "primaryTumorCuration", "Primary tumor search term not used");
    }

    @NotNull
    public static List<ValidationFinding> validateTreatmentCurator(@NotNull TreatmentCurator treatmentCurator) {
        return validateCurator(treatmentCurator, "treatmentCuration", "Treatment search term not used");
    }

    @NotNull
    private static List<ValidationFinding> validateCurator(@NotNull CleanableCurator curator, @NotNull String level,
            @NotNull String message) {
        List<ValidationFinding> findings = Lists.newArrayList();

        for (String unusedTerm : curator.unusedSearchTerms()) {
            findings.add(ValidationFinding.of(level, null, message, FormStatus.undefined(), unusedTerm));
        }

        return findings;
    }
}
