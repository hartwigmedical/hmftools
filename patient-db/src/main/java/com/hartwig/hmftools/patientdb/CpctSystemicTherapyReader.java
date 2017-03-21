package com.hartwig.hmftools.patientdb;

import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.List;
import java.util.Optional;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.ecrf.datamodel.EcrfPatient;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

class CpctSystemicTherapyReader {
    private static final Logger LOGGER = LogManager.getLogger(PatientDbRunner.class);
    private static final DateFormat dateFormat = new SimpleDateFormat("yyyy-MM-dd");

    @NotNull
    Optional<List<SystemicTherapyData>> read(@NotNull EcrfPatient patient) {
        final String hadSystemicTreatment = GenericReader.getField(patient, "BASELINE.PRETHERAPY.PRETHERAPY.SYSTEMIC");
        if (hadSystemicTreatment != null && hadSystemicTreatment.toLowerCase().equals("no")) {
            return Optional.empty();
        } else {
            final List<SystemicTherapyData> data = readData(patient);
            if (data.size() == 0) {
                return Optional.empty();
            } else {
                return Optional.of(data);
            }
        }
    }

    @NotNull
    private List<SystemicTherapyData> readData(@NotNull EcrfPatient patient) {
        final List<String> types = GenericReader.getFieldValues(patient,
                "BASELINE.PRETHERAPY.SYSTEMICTRT.SYSTEMICTYPE");
        final List<String> treatments = GenericReader.getFieldValues(patient,
                "BASELINE.PRETHERAPY.SYSTEMICTRT.SYSTEMICREG");
        final List<String> responses = GenericReader.getFieldValues(patient,
                "BASELINE.PRETHERAPY.SYSTEMICTRT.SYSTEMICRESP");
        final List<String> startDates = GenericReader.getFieldValues(patient,
                "BASELINE.PRETHERAPY.SYSTEMICTRT.SYSTEMICSTDTC");
        final List<String> endDates = GenericReader.getFieldValues(patient,
                "BASELINE.PRETHERAPY.SYSTEMICTRT.SYSTEMICENDTC");
        final int maxLength = Utils.getMaxLength(
                Lists.newArrayList(types, treatments, responses, startDates, endDates),
                "Not all systemic therapy fields contain the same number of values.");

        final List<SystemicTherapyData> therapies = Lists.newArrayList();
        for (int index = 0; index < maxLength; index++) {
            final String type = Utils.getElemAtIndex(types, index);
            final String treatment = Utils.getElemAtIndex(treatments, index);
            final String response = Utils.getElemAtIndex(responses, index);
            final Date startDate = Utils.getDate(Utils.getElemAtIndex(startDates, index), dateFormat);
            if (startDate == null) {
                LOGGER.warn(
                        "BASELINE.PRETHERAPY.SYSTEMICTRT.SYSTEMICSTDTC did not contain valid date at index " + index
                                + " for patient " + patient.patientId());
            }
            final Date endDate = Utils.getDate(Utils.getElemAtIndex(endDates, index), dateFormat);
            if (endDate == null) {
                LOGGER.warn(
                        "BASELINE.PRETHERAPY.SYSTEMICTRT.SYSTEMICENDTC did not contain valid date at index " + index
                                + " for patient  " + patient.patientId());
            }
            therapies.add(new SystemicTherapyData(startDate, endDate, type, treatment, response));
        }
        return therapies;
    }
}
