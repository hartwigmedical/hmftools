package com.hartwig.hmftools.patientdb.readers;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.ecrf.formstatus.FormStatus;
import com.hartwig.hmftools.patientdb.curators.PrimaryTumorCurator;
import com.hartwig.hmftools.patientdb.data.BaselineData;
import com.hartwig.hmftools.patientdb.data.CuratedPrimaryTumor;
import com.hartwig.hmftools.patientdb.data.ImmutableBaselineData;
import com.hartwig.hmftools.patientdb.data.ImmutablePreTreatmentData;
import com.hartwig.hmftools.patientdb.data.Patient;
import com.hartwig.hmftools.patientdb.data.PreTreatmentData;
import com.hartwig.hmftools.patientdb.data.SampleData;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class CorePatientReader {

    private static final Logger LOGGER = LogManager.getLogger(CorePatientReader.class);

    @NotNull
    private final PrimaryTumorCurator primaryTumorCurator;

    public CorePatientReader(@NotNull final PrimaryTumorCurator primaryTumorCurator) {
        this.primaryTumorCurator = primaryTumorCurator;
    }

    @NotNull
    public Patient read(@NotNull String patientIdentifier, @Nullable String limsPrimaryTumorLocation,
            @NotNull List<SampleData> sequencedSamples) {
        return new Patient(patientIdentifier,
                toBaselineData(limsPrimaryTumorLocation),
                noPreTreatmentData(),
                sequencedSamples,
                Lists.newArrayList(),
                Lists.newArrayList(),
                Lists.newArrayList(),
                Lists.newArrayList(),
                Lists.newArrayList(),
                Lists.newArrayList());
    }

    @NotNull
    private BaselineData toBaselineData(@Nullable String limsPrimaryTumorLocation) {
        CuratedPrimaryTumor curatedPrimaryTumor = primaryTumorCurator.search(limsPrimaryTumorLocation);
        if (curatedPrimaryTumor.location() == null && limsPrimaryTumorLocation != null
                && !limsPrimaryTumorLocation.isEmpty()) {
            LOGGER.warn("Could not curate CORE primary tumor '{}'", limsPrimaryTumorLocation);
        }

        return ImmutableBaselineData.builder()
                .registrationDate(null)
                .informedConsentDate(null)
                .gender(null)
                .hospital(null)
                .birthYear(null)
                .curatedPrimaryTumor(curatedPrimaryTumor)
                .deathDate(null)
                .demographyStatus(FormStatus.undefined())
                .primaryTumorStatus(FormStatus.undefined())
                .informedConsentStatus(FormStatus.undefined())
                .eligibilityStatus(FormStatus.undefined())
                .selectionCriteriaStatus(FormStatus.undefined())
                .deathStatus(FormStatus.undefined())
                .build();
    }

    @NotNull
    private static PreTreatmentData noPreTreatmentData() {
        return ImmutablePreTreatmentData.builder()
                .treatmentGiven(null)
                .radiotherapyGiven(null)
                .drugs(Lists.newArrayList())
                .formStatus(FormStatus.undefined())
                .build();
    }
}
