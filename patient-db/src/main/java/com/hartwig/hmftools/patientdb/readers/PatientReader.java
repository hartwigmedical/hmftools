package com.hartwig.hmftools.patientdb.readers;

import java.io.IOException;
import java.util.List;

import com.hartwig.hmftools.common.ecrf.datamodel.EcrfPatient;
import com.hartwig.hmftools.patientdb.data.Patient;
import com.hartwig.hmftools.patientdb.data.SampleData;

import org.jetbrains.annotations.NotNull;

public interface PatientReader {

    @NotNull
    Patient read(@NotNull final EcrfPatient ecrfPatient, @NotNull final List<SampleData> sequencedSamples) throws IOException;
}
