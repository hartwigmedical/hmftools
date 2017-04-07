package com.hartwig.hmftools.patientdb.readers;

import java.time.LocalDate;
import java.time.format.DateTimeFormatter;
import java.util.List;
import java.util.Optional;

import com.hartwig.hmftools.common.ecrf.datamodel.EcrfPatient;
import com.hartwig.hmftools.patientdb.Utils;
import com.hartwig.hmftools.patientdb.data.TreatmentData;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

class CpctTreatmentDataReader {
    private static final Logger LOGGER = LogManager.getLogger(CpctTreatmentDataReader.class);

    private static final String FIELD_TREATMENTNAME = "AFTERBIOPT.TRTAFTER.TRTAFTER.SYSREGPOST";
    private static final String FIELD_STARTDATE = "AFTERBIOPT.TRTAFTER.TRTAFTER.SYSSTDT";
    private static final String FIELD_ENDDATE = "AFTERBIOPT.TRTAFTER.TRTAFTER.SYSENDT";
    private static final String FIELD_RESPONSES = "TREATMENT.TUMORMEASUREMENT.TUMORMEASUREMENT.BESTRESPON";
    private static final String FIELD_RADIO_STARTDATE = "AFTERBIOPT.TRTAFTER.TRTAFTER.RADIOSTDTC";
    private static final String FIELD_RADIO_ENDDATE = "AFTERBIOPT.TRTAFTER.TRTAFTER.RADIOENDTC";

    private static final DateTimeFormatter dateFormatter = DateTimeFormatter.ofPattern("yyyy-MM-dd");

    @NotNull
    Optional<TreatmentData> read(@NotNull final EcrfPatient patient) {
        final List<String> treatmentNames = GenericReader.getFieldValues(patient, FIELD_TREATMENTNAME);
        final List<LocalDate> startDates = GenericReader.getDateFieldValues(patient, FIELD_STARTDATE, dateFormatter);
        final List<LocalDate> endDates = GenericReader.getDateFieldValues(patient, FIELD_ENDDATE, dateFormatter);
        final List<String> responses = GenericReader.getFieldValues(patient, FIELD_RESPONSES);
        final List<LocalDate> radiotherapyStartDates = GenericReader.getDateFieldValues(patient, FIELD_RADIO_STARTDATE,
                dateFormatter);
        final List<LocalDate> radiotherapyEndDates = GenericReader.getDateFieldValues(patient, FIELD_RADIO_ENDDATE,
                dateFormatter);

        final String treatmentName = Utils.getElemAtIndex(treatmentNames, 0);
        final String initialResponse = Utils.getElemAtIndex(responses, 0);
        final String firstResponse = Utils.getElemAtIndex(responses, 1);

        if (initialResponse != null && !initialResponse.replaceAll("\\s", "").toLowerCase().equals("ne")
                && !initialResponse.replaceAll("\\s", "").toLowerCase().equals("nd")
                && initialResponse.replaceAll("\\s", "").length() != 0) {
            LOGGER.warn("first value for field " + FIELD_RESPONSES + " was " + initialResponse
                    + " instead of empty or NE or ND for patient " + patient.patientId());
        }
        final LocalDate startDate = Utils.getElemAtIndex(startDates, 0);
        final LocalDate endDate = Utils.getElemAtIndex(endDates, 0);
        final LocalDate radioTherapyStartDate = Utils.getElemAtIndex(radiotherapyStartDates, 0);
        final LocalDate radioTherapyEndDate = Utils.getElemAtIndex(radiotherapyEndDates, 0);
        if (Utils.anyNotNull(startDate, endDate, treatmentName, firstResponse, radioTherapyStartDate,
                radioTherapyEndDate)) {
            return Optional.of(
                    new TreatmentData(startDate, endDate, treatmentName, firstResponse, radioTherapyStartDate,
                            radioTherapyEndDate));
        } else {
            return Optional.empty();
        }
    }
}
