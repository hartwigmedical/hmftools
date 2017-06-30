package com.hartwig.hmftools.patientdb.readers;

import java.time.LocalDate;
import java.time.format.DateTimeFormatter;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.ecrf.datamodel.EcrfForm;
import com.hartwig.hmftools.common.ecrf.datamodel.EcrfItemGroup;
import com.hartwig.hmftools.common.ecrf.datamodel.EcrfPatient;
import com.hartwig.hmftools.common.ecrf.datamodel.EcrfStudyEvent;
import com.hartwig.hmftools.patientdb.data.BiopsyTreatmentResponseData;

import org.jetbrains.annotations.NotNull;

public final class BiopsyTreatmentResponseReader {

    private static final String STUDY_TREATMENT = "SE.TREATMENT";
    public static final String FORM_TUMOR_MEASUREMENT = "FRM.TUMORMEASUREMENT";
    private static final String ITEMGROUP_MEASUREMENT = "GRP.TUMORMEASUREMENT.MEASUREMENT";
    public static final String FIELD_ASSESSMENT_DATE = "FLD.TUMORMEASUREMENT.ASSDTC";
    private static final String ITEMGROUP_TUMOR_MEASUREMENT = "GRP.TUMORMEASUREMENT.TUMORMEASUREMENT";
    public static final String FIELD_RESPONSE_DATE = "FLD.TUMORMEASUREMENT.RESPONSEDTC";
    public static final String FIELD_MEASUREMENT_YN = "FLD.TUMORMEASUREMENT.TMYN";
    public static final String FIELD_RESPONSE = "FLD.TUMORMEASUREMENT.BESTRESPON";

    private static final DateTimeFormatter DATE_FORMATTER = DateTimeFormatter.ofPattern("yyyy-MM-dd");

    private BiopsyTreatmentResponseReader() {
    }

    @NotNull
    static List<BiopsyTreatmentResponseData> read(@NotNull final EcrfPatient patient) {
        final List<BiopsyTreatmentResponseData> treatmentResponses = Lists.newArrayList();
        for (final EcrfStudyEvent studyEvent : patient.studyEventsPerOID(STUDY_TREATMENT)) {
            for (final EcrfForm form : studyEvent.nonEmptyFormsPerOID(FORM_TUMOR_MEASUREMENT, true)) {
                LocalDate assessmentDate = null;
                for (final EcrfItemGroup itemGroup : form.nonEmptyItemGroupsPerOID(ITEMGROUP_MEASUREMENT, true)) {
                    final LocalDate date = itemGroup.readItemDate(FIELD_ASSESSMENT_DATE, 0, DATE_FORMATTER, true);
                    if (date != null) {
                        assessmentDate = date;
                        break;
                    }
                }
                LocalDate responseDate = null;
                String measurementDone = null;
                String response = null;
                for (final EcrfItemGroup itemGroup : form.nonEmptyItemGroupsPerOID(ITEMGROUP_TUMOR_MEASUREMENT, true)) {
                    final LocalDate date = itemGroup.readItemDate(FIELD_RESPONSE_DATE, 0, DATE_FORMATTER, true);
                    if (date != null) {
                        responseDate = date;
                    }
                    final String measurementDoneValue = itemGroup.readItemString(FIELD_MEASUREMENT_YN, 0, true);
                    if (measurementDoneValue != null) {
                        measurementDone = measurementDoneValue;
                    }
                    final String responseValue = itemGroup.readItemString(FIELD_RESPONSE, 0, true);
                    if (responseValue != null) {
                        response = responseValue;
                    }
                }
                treatmentResponses.add(new BiopsyTreatmentResponseData(assessmentDate, responseDate, response, measurementDone));
            }
        }
        return treatmentResponses;
    }
}
