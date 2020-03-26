package com.hartwig.hmftools.patientdb.readers.wide;

import java.time.LocalDate;

import com.hartwig.hmftools.patientdb.curators.TumorLocationCurator;
import com.hartwig.hmftools.patientdb.data.BaselineData;
import com.hartwig.hmftools.patientdb.data.ImmutableBaselineData;

import org.jetbrains.annotations.NotNull;

public class BaselineReader {

    @NotNull
    private final TumorLocationCurator tumorLocationCurator;

    BaselineReader(@NotNull final TumorLocationCurator tumorLocationCurator) {
        this.tumorLocationCurator = tumorLocationCurator;
    }

    @NotNull
    BaselineData read() {
        return ImmutableBaselineData.builder()
                .registrationDate(null)
                .informedConsentStatus(null)
                .gender(null)
                .hospital(null)
                .birthYear(null)
                .curatedTumorLocation(null)
                .deathDate(null)
                .demographyStatus(null)
                .primaryTumorStatus(null)
                .informedConsentStatus(null)
                .eligibilityStatus(null)
                .selectionCriteriaStatus(null)
                .deathStatus(null)
                .build();
    }
}
