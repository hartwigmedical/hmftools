package com.hartwig.hmftools.patientdb.clinical.validators;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.patientdb.clinical.curators.CleanableCurator;
import com.hartwig.hmftools.patientdb.clinical.curators.PrimaryTumorCurator;
import com.hartwig.hmftools.patientdb.clinical.curators.TreatmentCurator;
import com.hartwig.hmftools.patientdb.clinical.ecrf.datamodel.ValidationFinding;
import com.hartwig.hmftools.patientdb.clinical.ecrf.formstatus.FormStatus;

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
