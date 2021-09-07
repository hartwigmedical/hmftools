package com.hartwig.hmftools.patientdb.clinical.readers.cpct;

import java.time.LocalDate;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.patientdb.clinical.datamodel.BiopsyTreatmentResponseData;
import com.hartwig.hmftools.patientdb.clinical.ecrf.datamodel.EcrfForm;
import com.hartwig.hmftools.patientdb.clinical.ecrf.datamodel.EcrfItemGroup;
import com.hartwig.hmftools.patientdb.clinical.ecrf.datamodel.EcrfPatient;
import com.hartwig.hmftools.patientdb.clinical.ecrf.datamodel.EcrfStudyEvent;

import org.jetbrains.annotations.NotNull;

final class BiopsyTreatmentResponseReader {

    static final String STUDY_TREATMENT = "SE.TREATMENT";
    static final String FORM_TUMOR_MEASUREMENT = "FRM.TUMORMEASUREMENT";
    private static final String ITEMGROUP_MEASUREMENT = "GRP.TUMORMEASUREMENT.MEASUREMENT";
    private static final String FIELD_ASSESSMENT_DATE = "FLD.TUMORMEASUREMENT.ASSDTC";
    static final String ITEMGROUP_TUMOR_MEASUREMENT = "GRP.TUMORMEASUREMENT.TUMORMEASUREMENT";
    private static final String FIELD_RESPONSE_DATE = "FLD.TUMORMEASUREMENT.RESPONSEDTC";
    private static final String FIELD_BONE_ONLY_DISEASE = "FLD.TUMORMEASUREMENT.BONEYN";
    static final String FIELD_MEASUREMENT_DONE = "FLD.TUMORMEASUREMENT.TMYN";
    static final String FIELD_RESPONSE = "FLD.TUMORMEASUREMENT.BESTRESPON";

    private BiopsyTreatmentResponseReader() {
    }

    @NotNull
    static List<BiopsyTreatmentResponseData> read(@NotNull EcrfPatient patient) {
        List<BiopsyTreatmentResponseData> treatmentResponses = Lists.newArrayList();
        for (EcrfStudyEvent studyEvent : patient.studyEventsPerOID(STUDY_TREATMENT)) {
            for (EcrfForm form : studyEvent.nonEmptyFormsPerOID(FORM_TUMOR_MEASUREMENT)) {
                // There are generally multiple assessment dates per tumor measurement (one per target lesion)
                LocalDate assessmentDate = null;
                for (EcrfItemGroup itemGroup : form.nonEmptyItemGroupsPerOID(ITEMGROUP_MEASUREMENT)) {
                    LocalDate date = itemGroup.readItemDate(FIELD_ASSESSMENT_DATE);
                    if (date != null) {
                        assessmentDate = date;
                        break;
                    }
                }

                LocalDate responseDate = null;
                String measurementDone = null;
                String boneOnlyDisease = null;
                String response = null;
                for (EcrfItemGroup itemGroup : form.nonEmptyItemGroupsPerOID(ITEMGROUP_TUMOR_MEASUREMENT)) {
                    LocalDate responseDateValue = itemGroup.readItemDate(FIELD_RESPONSE_DATE);
                    if (responseDateValue != null) {
                        responseDate = responseDateValue;
                    }

                    String measurementDoneValue = itemGroup.readItemString(FIELD_MEASUREMENT_DONE);
                    if (measurementDoneValue != null) {
                        measurementDone = measurementDoneValue;
                    }

                    String boneOnlyDiseaseValue = itemGroup.readItemString(FIELD_BONE_ONLY_DISEASE);
                    if (boneOnlyDiseaseValue != null) {
                        boneOnlyDisease = boneOnlyDiseaseValue;
                    }

                    String responseValue = itemGroup.readItemString(FIELD_RESPONSE);
                    if (responseValue != null) {
                        response = responseValue;
                    }
                }

                BiopsyTreatmentResponseData responseData = BiopsyTreatmentResponseData.of(null,
                        assessmentDate,
                        responseDate,
                        response,
                        measurementDone,
                        boneOnlyDisease,
                        form.status());
                if (!isEmpty(responseData)) {
                    treatmentResponses.add(responseData);
                }
            }
        }
        return treatmentResponses;
    }

    private static boolean isEmpty(@NotNull BiopsyTreatmentResponseData response) {
        return (response.measurementDone() == null && response.boneOnlyDisease() == null && response.date() == null
                && response.response() == null);
    }
}
