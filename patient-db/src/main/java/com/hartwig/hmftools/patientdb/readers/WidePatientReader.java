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
import com.hartwig.hmftools.patientdb.readers.wide.WideEcrfModel;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class WidePatientReader {

    @NotNull
    private final WideEcrfModel wideEcrfModel;
    @NotNull
    private final TumorLocationCurator tumorLocationCurator;

    public WidePatientReader(@NotNull final WideEcrfModel wideEcrfModel, @NotNull final TumorLocationCurator tumorLocationCurator) {
        this.wideEcrfModel = wideEcrfModel;
        this.tumorLocationCurator = tumorLocationCurator;
    }

    @NotNull
    public Patient read(@NotNull String patientIdentifier, @Nullable String primaryTumorLocation,
            @NotNull List<SampleData> sequencedSamples) {
        // TODO: Create the timeline based on the wideEcrfModel
        // TODO: Create patient object.
        return new Patient(patientIdentifier,
                toBaselineData(tumorLocationCurator.search(primaryTumorLocation)),
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
