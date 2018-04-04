package com.hartwig.hmftools.patientdb.readers.cpct;

import java.time.LocalDate;
import java.time.format.DateTimeFormatter;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.ecrf.datamodel.EcrfForm;
import com.hartwig.hmftools.common.ecrf.datamodel.EcrfItemGroup;
import com.hartwig.hmftools.common.ecrf.datamodel.EcrfPatient;
import com.hartwig.hmftools.common.ecrf.datamodel.EcrfStudyEvent;
import com.hartwig.hmftools.patientdb.data.BiopsyTreatmentResponseData;
import com.hartwig.hmftools.patientdb.data.ImmutableBiopsyTreatmentResponseData;

import org.jetbrains.annotations.NotNull;

public final class BiopsyTreatmentResponseReader {

    static final String STUDY_TREATMENT = "SE.TREATMENT";
    public static final String FORM_TUMOR_MEASUREMENT = "FRM.TUMORMEASUREMENT";
    private static final String ITEMGROUP_MEASUREMENT = "GRP.TUMORMEASUREMENT.MEASUREMENT";
    public static final String FIELD_ASSESSMENT_DATE = "FLD.TUMORMEASUREMENT.ASSDTC";
    static final String ITEMGROUP_TUMOR_MEASUREMENT = "GRP.TUMORMEASUREMENT.TUMORMEASUREMENT";
    public static final String FIELD_RESPONSE_DATE = "FLD.TUMORMEASUREMENT.RESPONSEDTC";
    private static final String FIELD_BONE_ONLY_DISEASE = "FLD.TUMORMEASUREMENT.BONEYN";
    public static final String FIELD_MEASUREMENT_DONE = "FLD.TUMORMEASUREMENT.TMYN";
    public static final String FIELD_RESPONSE = "FLD.TUMORMEASUREMENT.BESTRESPON";

    private static final DateTimeFormatter DATE_FORMATTER = DateTimeFormatter.ofPattern("yyyy-MM-dd");

    private BiopsyTreatmentResponseReader() {
    }

    @NotNull
    static List<BiopsyTreatmentResponseData> read(@NotNull final EcrfPatient patient) {
        final List<BiopsyTreatmentResponseData> treatmentResponses = Lists.newArrayList();
        for (final EcrfStudyEvent studyEvent : patient.studyEventsPerOID(STUDY_TREATMENT)) {
            for (final EcrfForm form : studyEvent.nonEmptyFormsPerOID(FORM_TUMOR_MEASUREMENT, false)) {
                // KODU: There are generally multiple assessment dates per tumor measurement (one per target lesion)
                LocalDate assessmentDate = null;
                for (final EcrfItemGroup itemGroup : form.nonEmptyItemGroupsPerOID(ITEMGROUP_MEASUREMENT, false)) {
                    final LocalDate date = itemGroup.readItemDate(FIELD_ASSESSMENT_DATE, 0, DATE_FORMATTER, false);
                    if (date != null) {
                        assessmentDate = date;
                        break;
                    }
                }
                LocalDate responseDate = null;
                String measurementDone = null;
                String boneOnlyDisease = null;
                String response = null;
                for (final EcrfItemGroup itemGroup : form.nonEmptyItemGroupsPerOID(ITEMGROUP_TUMOR_MEASUREMENT, false)) {
                    final LocalDate responseDateValue = itemGroup.readItemDate(FIELD_RESPONSE_DATE, 0, DATE_FORMATTER, false);
                    if (responseDateValue != null) {
                        responseDate = responseDateValue;
                    }

                    final String measurementDoneValue = itemGroup.readItemString(FIELD_MEASUREMENT_DONE, 0, false);
                    if (measurementDoneValue != null) {
                        measurementDone = measurementDoneValue;
                    }

                    final String boneOnlyDiseaseValue = itemGroup.readItemString(FIELD_BONE_ONLY_DISEASE, 0, false);
                    if (boneOnlyDiseaseValue != null) {
                        boneOnlyDisease = boneOnlyDiseaseValue;
                    }

                    final String responseValue = itemGroup.readItemString(FIELD_RESPONSE, 0, false);
                    if (responseValue != null) {
                        response = responseValue;
                    }
                }
                BiopsyTreatmentResponseData responseData = ImmutableBiopsyTreatmentResponseData.of(assessmentDate,
                        responseDate,
                        response,
                        measurementDone,
                        boneOnlyDisease,
                        form.status(),
                        form.locked());
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
