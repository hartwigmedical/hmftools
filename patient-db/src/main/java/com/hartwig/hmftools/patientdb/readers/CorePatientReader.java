package com.hartwig.hmftools.patientdb.readers;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.ecrf.formstatus.FormStatus;
import com.hartwig.hmftools.patientdb.curators.TumorLocationCurator;
import com.hartwig.hmftools.patientdb.data.BaselineData;
import com.hartwig.hmftools.patientdb.data.CuratedTumorLocation;
import com.hartwig.hmftools.patientdb.data.ImmutableBaselineData;
import com.hartwig.hmftools.patientdb.data.ImmutablePreTreatmentData;
import com.hartwig.hmftools.patientdb.data.Patient;
import com.hartwig.hmftools.patientdb.data.PreTreatmentData;
import com.hartwig.hmftools.patientdb.data.SampleData;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class CorePatientReader {

    @NotNull
    private final TumorLocationCurator tumorLocationCurator;

    public CorePatientReader(@NotNull final TumorLocationCurator tumorLocationCurator) {
        this.tumorLocationCurator = tumorLocationCurator;
    }

    @NotNull
    public Patient read(@NotNull String patientIdentifier, @Nullable String limsPrimaryTumorLocation,
            @NotNull List<SampleData> sequencedSamples) {
        return new Patient(patientIdentifier,
                toBaselineData(tumorLocationCurator.search(limsPrimaryTumorLocation)),
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
    private static BaselineData toBaselineData(@NotNull CuratedTumorLocation curatedTumorLocation) {
        return ImmutableBaselineData.of(null,
                null,
                null,
                null,
                null,
                curatedTumorLocation,
                null,
                FormStatus.undefined(),
                FormStatus.undefined(),
                FormStatus.undefined(),
                FormStatus.undefined(),
                FormStatus.undefined(),
                FormStatus.undefined());
    }

    @NotNull
    private static PreTreatmentData noPreTreatmentData() {
        return ImmutablePreTreatmentData.of(null, null, Lists.newArrayList(), FormStatus.undefined());
    }
}
