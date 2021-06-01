package com.hartwig.hmftools.patientdb.clinical.readers;

import java.util.List;

import com.hartwig.hmftools.patientdb.clinical.datamodel.Patient;
import com.hartwig.hmftools.patientdb.clinical.datamodel.SampleData;
import com.hartwig.hmftools.patientdb.clinical.ecrf.datamodel.EcrfPatient;

import org.jetbrains.annotations.NotNull;

public interface EcrfPatientReader {

    @NotNull
    Patient read(@NotNull EcrfPatient ecrfPatient, @NotNull List<SampleData> sequencedSamples);
}
