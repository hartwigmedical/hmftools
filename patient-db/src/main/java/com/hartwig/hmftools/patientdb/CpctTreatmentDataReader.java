package com.hartwig.hmftools.patientdb;

import java.time.LocalDate;
import java.time.format.DateTimeFormatter;
import java.util.List;
import java.util.Optional;

import com.hartwig.hmftools.common.ecrf.datamodel.EcrfPatient;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

class CpctTreatmentDataReader {
    private static final String FIELD_TREATMENTNAME = "AFTERBIOPT.TRTAFTER.TRTAFTER.SYSREGPOST";
    private static final String FIELD_STARTDATE = "AFTERBIOPT.TRTAFTER.TRTAFTER.SYSSTDT";
    private static final String FIELD_ENDDATE = "AFTERBIOPT.TRTAFTER.TRTAFTER.SYSENDT";
    private static final String FIELD_RESPONSES = "TREATMENT.TUMORMEASUREMENT.TUMORMEASUREMENT.BESTRESPON";

    private static final Logger LOGGER = LogManager.getLogger(PatientDbRunner.class);
    private static final DateTimeFormatter dateFormatter = DateTimeFormatter.ofPattern("yyyy-MM-dd");

    @NotNull
    Optional<TreatmentData> read(@NotNull EcrfPatient patient) {
        final List<String> treatmentNames = GenericReader.getFieldValues(patient, FIELD_TREATMENTNAME);
        final List<String> startDates = GenericReader.getFieldValues(patient, FIELD_STARTDATE);
        final List<String> endDates = GenericReader.getFieldValues(patient, FIELD_ENDDATE);
        final List<String> responses = GenericReader.getFieldValues(patient, FIELD_RESPONSES);

        final String treatmentName = Utils.getElemAtIndex(treatmentNames, 0);
        final String initialResponse = Utils.getElemAtIndex(responses, 0);
        final String firstResponse = Utils.getElemAtIndex(responses, 1);

        if (initialResponse != null && !initialResponse.replaceAll("\\s", "").toLowerCase().equals("ne")
                && !initialResponse.replaceAll("\\s", "").toLowerCase().equals("nd")
                && initialResponse.replaceAll("\\s", "").length() != 0) {
            LOGGER.warn("first value for field " + FIELD_RESPONSES + " was " + initialResponse
                    + " instead of empty or NE for patient " + patient.patientId());
        }
        final LocalDate startDate = Utils.getDate(Utils.getElemAtIndex(startDates, 0), dateFormatter);
        if (startDate == null) {
            LOGGER.warn(FIELD_STARTDATE + " did not contain valid date at index 0 " + " for patient "
                    + patient.patientId());
        }
        final LocalDate endDate = Utils.getDate(Utils.getElemAtIndex(endDates, 0), dateFormatter);
        if (endDate == null) {
            LOGGER.warn(FIELD_ENDDATE + " did not contain valid date at index 0 " + " for patient  "
                    + patient.patientId());
        }
        if (startDate == null && endDate == null && (treatmentName == null
                || treatmentName.replaceAll("\\s", "").length() == 0) && firstResponse == null) {
            return Optional.empty();
        } else
            return Optional.of(new TreatmentData(startDate, endDate, treatmentName, firstResponse));
    }
}
