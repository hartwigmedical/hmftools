package com.hartwig.hmftools.patientdb;

import java.text.DateFormat;
import java.text.ParseException;
import java.text.SimpleDateFormat;
import java.util.Collection;
import java.util.Date;
import java.util.Iterator;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.ecrf.datamodel.EcrfPatient;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

class CpctSystemicTherapyReader {
    private static final Logger LOGGER = LogManager.getLogger(PatientDbRunner.class);
    private static final DateFormat dateFormat = new SimpleDateFormat("yyyy-MM-dd");

    @NotNull
    List<SystemicTherapyData> read(@NotNull EcrfPatient patient) {
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
        final int maxLength = getMaxLength(Lists.newArrayList(types, treatments, responses, startDates, endDates));

        final List<SystemicTherapyData> therapies = Lists.newArrayList();
        for (int index = 0; index < maxLength; index++) {
            final String type = getElemAtIndex(types, index);
            final String treatment = getElemAtIndex(treatments, index);
            final String response = getElemAtIndex(responses, index);
            final Date startDate = getDate(getElemAtIndex(startDates, index));
            if (startDate == null) {
                LOGGER.warn(
                        "BASELINE.PRETHERAPY.SYSTEMICTRT.SYSTEMICSTDTC did not contain valid date at index " + index
                                + " for patient" + patient.patientId());
            }
            final Date endDate = getDate(getElemAtIndex(endDates, index));
            if (endDate == null) {
                LOGGER.warn(
                        "BASELINE.PRETHERAPY.SYSTEMICTRT.SYSTEMICENDTC did not contain valid date at index " + index
                                + " for patient" + patient.patientId());
            }
            therapies.add(new SystemicTherapyData(startDate, endDate, type, treatment, response));
        }
        return therapies;
    }

    private int getMaxLength(@NotNull Collection<List<String>> lists) {
        final Iterator<List<String>> listIterator = lists.iterator();
        int maxSize = -1;
        while (listIterator.hasNext()) {
            List<String> currentList = listIterator.next();
            if (currentList != null) {
                if (maxSize < currentList.size()) {
                    if (maxSize != -1) {
                        LOGGER.warn("Not all systemic therapy fields contain the same number of values.");
                    }
                    maxSize = currentList.size();
                }
            }
        }
        return maxSize;
    }

    @Nullable
    private Date getDate(@Nullable String dateFieldValue) {
        try {
            return dateFormat.parse(dateFieldValue);
        } catch (ParseException | NullPointerException e) {
            return null;
        }
    }

    @Nullable
    private String getElemAtIndex(@Nullable List<String> list, int index) {
        try {
            return list.get(index);
        } catch (NullPointerException | IndexOutOfBoundsException e) {
            return null;
        }
    }
}
