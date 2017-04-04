package com.hartwig.hmftools.patientdb;

import java.time.LocalDate;
import java.time.format.DateTimeFormatter;
import java.util.List;
import java.util.Optional;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.ecrf.datamodel.EcrfPatient;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

class CpctSystemicTherapyReader {
    private static final Logger LOGGER = LogManager.getLogger(CpctSystemicTherapyReader.class);

    private static final String FIELD_HADSYSTEMICTREATMENT = "BASELINE.PRETHERAPY.PRETHERAPY.SYSTEMIC";
    private static final String FIELD_SYSTEMICTYPE = "BASELINE.PRETHERAPY.SYSTEMICTRT.SYSTEMICTYPE";
    private static final String FIELD_SYSTEMICTYPEOTHER = "BASELINE.PRETHERAPY.SYSTEMICTRT.SYSTEMICTYPESP";
    private static final String FIELD_TREATMENT = "BASELINE.PRETHERAPY.SYSTEMICTRT.SYSTEMICREG";
    private static final String FIELD_RESPONSES = "BASELINE.PRETHERAPY.SYSTEMICTRT.SYSTEMICRESP";
    private static final String FIELD_STARTDATE = "BASELINE.PRETHERAPY.SYSTEMICTRT.SYSTEMICSTDTC";
    private static final String FIELD_ENDDATE = "BASELINE.PRETHERAPY.SYSTEMICTRT.SYSTEMICENDTC";

    private static final DateTimeFormatter dateFormatter = DateTimeFormatter.ofPattern("yyyy-MM-dd");

    @NotNull
    Optional<List<SystemicTherapyData>> read(@NotNull final EcrfPatient patient) {
        final String hadSystemicTreatment = GenericReader.getField(patient, FIELD_HADSYSTEMICTREATMENT);
        if (hadSystemicTreatment != null && hadSystemicTreatment.toLowerCase().equals("no")) {
            return Optional.empty();
        } else {
            final List<SystemicTherapyData> data = readData(patient);
            if (data.size() == 0) {
                LOGGER.warn(patient.patientId() + " had value " + hadSystemicTreatment + " instead of No for field "
                        + FIELD_HADSYSTEMICTREATMENT + " but systemic treatment fields contained no values.");
                return Optional.empty();
            } else {
                return Optional.of(data);
            }
        }
    }

    @NotNull
    private List<SystemicTherapyData> readData(@NotNull final EcrfPatient patient) {
        final List<String> types = GenericReader.getFieldValuesWithOthers(patient, FIELD_SYSTEMICTYPE,
                FIELD_SYSTEMICTYPEOTHER);
        final List<String> treatments = GenericReader.getFieldValues(patient, FIELD_TREATMENT);
        final List<String> responses = GenericReader.getFieldValues(patient, FIELD_RESPONSES);
        final List<String> startDates = GenericReader.getFieldValues(patient, FIELD_STARTDATE);
        final List<String> endDates = GenericReader.getFieldValues(patient, FIELD_ENDDATE);

        final int maxLength = Utils.getMaxLength(
                Lists.newArrayList(types, treatments, responses, startDates, endDates),
                "Not all systemic therapy fields contain the same number of values.");

        final List<SystemicTherapyData> therapies = Lists.newArrayList();
        for (int index = 0; index < maxLength; index++) {
            final String type = Utils.getElemAtIndex(types, index);
            final String treatment = Utils.getElemAtIndex(treatments, index);
            final String response = Utils.getElemAtIndex(responses, index);
            final LocalDate startDate = Utils.getDate(Utils.getElemAtIndex(startDates, index), dateFormatter);
            if (startDate == null) {
                LOGGER.warn(FIELD_STARTDATE + " did not contain valid date at index " + index + " for patient "
                        + patient.patientId());
            }
            final LocalDate endDate = Utils.getDate(Utils.getElemAtIndex(endDates, index), dateFormatter);
            if (endDate == null) {
                LOGGER.warn(FIELD_ENDDATE + " did not contain valid date at index " + index + " for patient  "
                        + patient.patientId());
            }
            therapies.add(new SystemicTherapyData(startDate, endDate, type, treatment, response));
        }
        return therapies;
    }

}
