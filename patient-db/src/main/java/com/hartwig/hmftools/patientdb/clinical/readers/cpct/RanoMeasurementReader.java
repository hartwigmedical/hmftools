package com.hartwig.hmftools.patientdb.clinical.readers.cpct;

import java.time.LocalDate;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.patientdb.clinical.datamodel.ImmutableRanoMeasurementData;
import com.hartwig.hmftools.patientdb.clinical.datamodel.RanoMeasurementData;
import com.hartwig.hmftools.patientdb.clinical.ecrf.datamodel.EcrfForm;
import com.hartwig.hmftools.patientdb.clinical.ecrf.datamodel.EcrfItemGroup;
import com.hartwig.hmftools.patientdb.clinical.ecrf.datamodel.EcrfPatient;
import com.hartwig.hmftools.patientdb.clinical.ecrf.datamodel.EcrfStudyEvent;

import org.jetbrains.annotations.NotNull;

final class RanoMeasurementReader {

    private static final String STUDY_TREATMENT = "SE.TREATMENT";
    private static final String FORM_RESPONSE_RANO = "FRM.TMRANO";
    private static final String ITEMGROUP_RESPONSE_RANO = "GRP.TMRANO.TMRANO";

    private static final String FIELD_THERAPY_GIVEN = "FLD.TMRANO.RNTMRNYN";
    private static final String FIELD_RESPONSE_DATE = "FLD.TMRANO.RNTMDT";
    private static final String FIELD_TARGET_LESION_RESPONSE = "FLD.TMRANO.RNRSPTL";
    private static final String FIELD_NO_TARGET_LESION_RESPONSE = "FLD.TMRANO.RNRSPNTL";
    private static final String FIELD_OVERALL_RESPONSE = "FLD.TMRANO.RNVRLL";

    private RanoMeasurementReader() {
    }

    @NotNull
    static List<RanoMeasurementData> read(@NotNull EcrfPatient patient) {
        List<RanoMeasurementData> rano = Lists.newArrayList();
        for (EcrfStudyEvent studyEvent : patient.studyEventsPerOID(STUDY_TREATMENT)) {
            for (EcrfForm form : studyEvent.nonEmptyFormsPerOID(FORM_RESPONSE_RANO)) {
                for (EcrfItemGroup responseGroup : form.nonEmptyItemGroupsPerOID(ITEMGROUP_RESPONSE_RANO)) {
                    String therapyGiven = responseGroup.readItemString(FIELD_THERAPY_GIVEN);
                    LocalDate responseDate = responseGroup.readItemDate(FIELD_RESPONSE_DATE);
                    String targetLesionResponse = responseGroup.readItemString(FIELD_TARGET_LESION_RESPONSE);
                    String noTargetLesionResponse = responseGroup.readItemString(FIELD_NO_TARGET_LESION_RESPONSE);
                    String overallResponse = responseGroup.readItemString(FIELD_OVERALL_RESPONSE);

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