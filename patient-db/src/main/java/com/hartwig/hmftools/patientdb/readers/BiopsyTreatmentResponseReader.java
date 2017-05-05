package com.hartwig.hmftools.patientdb.readers;

import java.time.LocalDate;
import java.time.format.DateTimeFormatter;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.ecrf.datamodel.EcrfForm;
import com.hartwig.hmftools.common.ecrf.datamodel.EcrfItemGroup;
import com.hartwig.hmftools.common.ecrf.datamodel.EcrfPatient;
import com.hartwig.hmftools.common.ecrf.datamodel.EcrfStudyEvent;
import com.hartwig.hmftools.patientdb.Utils;
import com.hartwig.hmftools.patientdb.data.BiopsyTreatmentResponseData;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class BiopsyTreatmentResponseReader {
    private static final Logger LOGGER = LogManager.getLogger(BiopsyTreatmentResponseReader.class);

    private static final String STUDY_TREATMENT = "SE.TREATMENT";
    private static final String FORM_TUMOR_MEASUREMENT = "FRM.TUMORMEASUREMENT";
    private static final String ITEMGROUP_MEASUREMENT = "GRP.TUMORMEASUREMENT.MEASUREMENT";
    private static final String FIELD_ASSESSMENT_DATE = "FLD.TUMORMEASUREMENT.ASSDTC";
    private static final String ITEMGROUP_TUMOR_MEASUREMENT = "GRP.TUMORMEASUREMENT.TUMORMEASUREMENT";
    private static final String FIELD_RESPONSE_DATE = "FLD.TUMORMEASUREMENT.RESPONSEDTC";
    private static final String FIELD_MEASUREMENT_YN = "FLD.TUMORMEASUREMENT.TMYN";
    private static final String FIELD_RESPONSE = "FLD.TUMORMEASUREMENT.BESTRESPON";

    private static final DateTimeFormatter dateFormatter = DateTimeFormatter.ofPattern("yyyy-MM-dd");

    @NotNull
    public static List<BiopsyTreatmentResponseData> read(@NotNull final EcrfPatient patient) {
        final List<BiopsyTreatmentResponseData> treatmentResponses = Lists.newArrayList();
        for (final EcrfStudyEvent studyEvent : patient.studyEventsPerOID(STUDY_TREATMENT)) {
            for (final EcrfForm form : studyEvent.formsPerOID(FORM_TUMOR_MEASUREMENT)) {
                if (form.isEmpty()) {
                    LOGGER.warn(
                            "Ignoring empty form: " + FORM_TUMOR_MEASUREMENT + " for patient " + patient.patientId());
                    continue;
                }
                LocalDate assessmentDate = null;
                for (final EcrfItemGroup itemGroup : form.itemGroupsPerOID(ITEMGROUP_MEASUREMENT)) {
                    for (final String assesmentDateValue : itemGroup.itemsPerOID(FIELD_ASSESSMENT_DATE)) {
                        final LocalDate date = Utils.getDate(assesmentDateValue, dateFormatter);
                        if (date != null) {
                            assessmentDate = date;
                            break;
                        }
                    }
                }
                LocalDate responseDate = null;
                String measurementDone = null;
                String response = null;
                for (final EcrfItemGroup itemGroup : form.itemGroupsPerOID(ITEMGROUP_TUMOR_MEASUREMENT)) {
                    if (itemGroup.isEmpty()) {
                        LOGGER.warn("Ignoring empty item group: " + ITEMGROUP_TUMOR_MEASUREMENT + " for patient "
                                + patient.patientId());
                        continue;
                    }
                    for (final String responseDateValue : itemGroup.itemsPerOID(FIELD_RESPONSE_DATE)) {
                        final LocalDate date = Utils.getDate(responseDateValue, dateFormatter);
                        if (date != null) {
                            responseDate = date;
                            break;
                        }
                    }
                    for (final String measurementDoneValue : itemGroup.itemsPerOID(FIELD_MEASUREMENT_YN)) {
                        if (measurementDoneValue != null) {
                            measurementDone = measurementDoneValue;
                            break;
                        }
                    }
                    for (final String responseValue : itemGroup.itemsPerOID(FIELD_RESPONSE)) {
                        if (responseValue != null) {
                            response = responseValue;
                            break;
                        }
                    }
                }
                if (assessmentDate == null && responseDate == null) {
                    LOGGER.warn("Both assessment date and treatment response date missing for patient: "
                            + patient.patientId());
                } else {
                    if (assessmentDate != null && responseDate != null && !assessmentDate.isEqual(responseDate)) {
                        LOGGER.warn("Assessment date(" + assessmentDate + ") differs from treatment response date("
                                + responseDate + ") for patient: " + patient.patientId());
                    }
                    treatmentResponses.add(
                            new BiopsyTreatmentResponseData(assessmentDate, responseDate, response, measurementDone));
                }
            }
        }
        return treatmentResponses;
    }
}
