package com.hartwig.hmftools.patientdb;

import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.List;
import java.util.Optional;

import com.hartwig.hmftools.common.ecrf.datamodel.EcrfPatient;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class CpctTreatmentDataReader {
    private static final Logger LOGGER = LogManager.getLogger(PatientDbRunner.class);
    private static final DateFormat dateFormat = new SimpleDateFormat("yyyy-MM-dd");

    @NotNull
    Optional<TreatmentData> read(@NotNull EcrfPatient patient) {
        final List<String> treatmentNames = GenericReader.getFieldValues(patient,
                "AFTERBIOPT.TRTAFTER.TRTAFTER.SYSREGPOST");
        final List<String> startDates = GenericReader.getFieldValues(patient, "AFTERBIOPT.TRTAFTER.TRTAFTER.SYSSTDT");
        final List<String> endDates = GenericReader.getFieldValues(patient, "AFTERBIOPT.TRTAFTER.TRTAFTER.SYSENDT");
        final List<String> responses = GenericReader.getFieldValues(patient,
                "TREATMENT.TUMORMEASUREMENT.TUMORMEASUREMENT.BESTRESPON");

        final String treatmentName = Utils.getElemAtIndex(treatmentNames, 0);
        final String initialResponse = Utils.getElemAtIndex(responses, 0);
        final String firstResponse = Utils.getElemAtIndex(responses, 1);

        if (initialResponse != null && !initialResponse.replaceAll("\\s", "").toLowerCase().equals("ne")
                && !initialResponse.replaceAll("\\s", "").toLowerCase().equals("nd")
                && initialResponse.replaceAll("\\s", "").length() != 0) {
            LOGGER.warn("first value for field TREATMENT.TUMORMEASUREMENT.TUMORMEASUREMENT.BESTRESPON was "
                    + initialResponse + " instead of empty or NE for patient " + patient.patientId());
        }
        final Date startDate = Utils.getDate(Utils.getElemAtIndex(startDates, 0), dateFormat);
        if (startDate == null) {
            LOGGER.warn("BASELINE.PRETHERAPY.SYSTEMICTRT.SYSTEMICSTDTC did not contain valid date at index 0 "
                    + " for patient " + patient.patientId());
        }
        final Date endDate = Utils.getDate(Utils.getElemAtIndex(endDates, 0), dateFormat);
        if (endDate == null) {
            LOGGER.warn("BASELINE.PRETHERAPY.SYSTEMICTRT.SYSTEMICENDTC did not contain valid date at index 0 "
                    + " for patient  " + patient.patientId());
        }
        if (startDate == null && endDate == null && (treatmentName == null
                || treatmentName.replaceAll("\\s", "").length() == 0) && firstResponse == null) {
            return Optional.empty();
        } else
            return Optional.of(new TreatmentData(startDate, endDate, treatmentName, firstResponse));
    }
}
