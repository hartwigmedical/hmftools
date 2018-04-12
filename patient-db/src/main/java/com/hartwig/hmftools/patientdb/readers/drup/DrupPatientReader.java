package com.hartwig.hmftools.patientdb.readers.drup;

import java.util.List;

import com.hartwig.hmftools.common.ecrf.datamodel.EcrfPatient;
import com.hartwig.hmftools.patientdb.data.Patient;
import com.hartwig.hmftools.patientdb.data.SampleData;
import com.hartwig.hmftools.patientdb.readers.PatientReader;

import org.jetbrains.annotations.NotNull;

public class DrupPatientReader implements PatientReader {

    public DrupPatientReader() {
    }

    @NotNull
    @Override
    public Patient read(@NotNull final EcrfPatient ecrfPatient, @NotNull final List<SampleData> sequencedSamples) {
        return null;
    }
}
