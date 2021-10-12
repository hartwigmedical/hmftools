package com.hartwig.hmftools.patientdb.clinical.readers;

import java.io.IOException;
import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.patientdb.clinical.consents.ConsentConfig;
import com.hartwig.hmftools.patientdb.clinical.datamodel.Patient;
import com.hartwig.hmftools.patientdb.clinical.datamodel.SampleData;
import com.hartwig.hmftools.patientdb.clinical.ecrf.datamodel.EcrfPatient;

import org.jetbrains.annotations.NotNull;

public interface EcrfPatientReader {

    @NotNull
    Patient read(@NotNull EcrfPatient ecrfPatient, @NotNull List<SampleData> sequencedSamples,
            @NotNull Map<String, ConsentConfig> consentConfigMap, @NotNull String cohortId) throws IOException;
}
