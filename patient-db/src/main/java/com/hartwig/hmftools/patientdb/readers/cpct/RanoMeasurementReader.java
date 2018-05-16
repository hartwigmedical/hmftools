package com.hartwig.hmftools.patientdb.readers.cpct;

import java.time.LocalDate;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.ecrf.datamodel.EcrfForm;
import com.hartwig.hmftools.common.ecrf.datamodel.EcrfItemGroup;
import com.hartwig.hmftools.common.ecrf.datamodel.EcrfPatient;
import com.hartwig.hmftools.common.ecrf.datamodel.EcrfStudyEvent;
import com.hartwig.hmftools.patientdb.data.ImmutableRanoMeasurementData;
import com.hartwig.hmftools.patientdb.data.RanoMeasurementData;

import org.jetbrains.annotations.NotNull;

final class RanoMeasurementReader {
    private static final String STUDY_TREATMENT = "SE.TREATMENT";
    private static final String FORM_RESPONSERANO = "FRM.TMRANO";
    private static final String ITEMGROUP_RESPONSERANO = "GRP.TMRANO";

    private static final String FIELD_RANOTHERAPYGIVEN = "FLD.RNTMRNYN";
    private static final String FIELD_RANORESPONSEDATE = "FLD.RNTMDT";
    private static final String FIELD_RANOTARGETLESIONRESPONSE = "FLD.RNRSPTL";
    private static final String FIELD_RANONOTARGETLESIONRESPONSE = "FLD.RNRSPNTL";
    private static final String FIELD_RANOOVERALLRESPONSE = "FLD.RNRSPNTL";

    private RanoMeasurementReader() {
    }

    @NotNull
    static List<RanoMeasurementData> read(@NotNull final EcrfPatient patient) {
        final List<RanoMeasurementData> rano = Lists.newArrayList();
        for (final EcrfStudyEvent studyEvent : patient.studyEventsPerOID(STUDY_TREATMENT)) {
            for (final EcrfForm form : studyEvent.nonEmptyFormsPerOID(FORM_RESPONSERANO)) {
                for (final EcrfItemGroup responseGroup : form.nonEmptyItemGroupsPerOID(ITEMGROUP_RESPONSERANO)) {
                    String therapyGiven = responseGroup.readItemString(FIELD_RANOTHERAPYGIVEN);
                    LocalDate responseDate = responseGroup.readItemDate(FIELD_RANORESPONSEDATE);
                    String targetLesionResponse = responseGroup.readItemString(FIELD_RANOTARGETLESIONRESPONSE);
                    String noTargetLesionResponse = responseGroup.readItemString(FIELD_RANONOTARGETLESIONRESPONSE);
                    String overallResponse = responseGroup.readItemString(FIELD_RANOOVERALLRESPONSE);
                    rano.add(ImmutableRanoMeasurementData.of(patient.patientId(),
                            therapyGiven,
                            responseDate,
                            targetLesionResponse,
                            noTargetLesionResponse,
                            overallResponse,
                            form.status()));
                }
            }
        }
        return rano;
    }
}