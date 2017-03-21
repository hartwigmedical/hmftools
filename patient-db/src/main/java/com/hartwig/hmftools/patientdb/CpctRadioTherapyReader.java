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

class CpctRadioTherapyReader {
    private static final Logger LOGGER = LogManager.getLogger(PatientDbRunner.class);
    private static final DateFormat dateFormat = new SimpleDateFormat("yyyy-MM-dd");

    @NotNull
    Optional<List<RadioTherapyData>> read(@NotNull EcrfPatient patient) {
        final String hadRadiotherapy = GenericReader.getField(patient, "BASELINE.PRETHERAPY.PRETHERAPY.RADIOTHER");
        if (hadRadiotherapy != null && hadRadiotherapy.toLowerCase().equals("no")) {
            return Optional.empty();
        } else {
            final List<RadioTherapyData> data = readData(patient);
            if (data.size() == 0) {
                return Optional.empty();
            } else {
                return Optional.of(data);
            }
        }
    }

    @NotNull
    private List<RadioTherapyData> readData(@NotNull EcrfPatient patient) {
        final List<String> endDates = GenericReader.getFieldValues(patient, "BASELINE.PRETHERAPY.RTP.RADIOTHERENDTC");
        final List<String> sites = GenericReader.getFieldValues(patient, "BASELINE.PRETHERAPY.RTP.RADIOTHERSITE");
        final int maxLength = Utils.getMaxLength(Lists.newArrayList(endDates, sites),
                "Not all radio therapy fields contain the same number of values.");

        final List<RadioTherapyData> therapies = Lists.newArrayList();
        for (int index = 0; index < maxLength; index++) {
            final String site = Utils.getElemAtIndex(sites, index);
            final Date endDate = Utils.getDate(Utils.getElemAtIndex(endDates, index), dateFormat);
            if (endDate == null) {
                LOGGER.warn(
                        "BASELINE.PRETHERAPY.SYSTEMICTRT.SYSTEMICENDTC did not contain valid date at index " + index
                                + " for patient  " + patient.patientId());
            }
            therapies.add(new RadioTherapyData(endDate, site));
        }
        return therapies;
    }

}
