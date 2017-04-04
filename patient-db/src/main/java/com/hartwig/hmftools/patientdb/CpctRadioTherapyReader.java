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

class CpctRadioTherapyReader {
    private static final Logger LOGGER = LogManager.getLogger(CpctRadioTherapyReader.class);

    private static final String FIELD_RADIOENDDATE = "BASELINE.PRETHERAPY.RTP.RADIOTHERENDTC";
    private static final String FIELD_RADIOSITE = "BASELINE.PRETHERAPY.RTP.RADIOTHERSITE";
    private static final String FIELD_HADRADIOTHERAPY = "BASELINE.PRETHERAPY.PRETHERAPY.RADIOTHER";

    private static final DateTimeFormatter dateFormatter = DateTimeFormatter.ofPattern("yyyy-MM-dd");

    @NotNull
    Optional<List<RadioTherapyData>> read(@NotNull final EcrfPatient patient) {
        final String hadRadiotherapy = GenericReader.getField(patient, FIELD_HADRADIOTHERAPY);
        if (hadRadiotherapy != null && hadRadiotherapy.toLowerCase().equals("no")) {
            return Optional.empty();
        } else {
            final List<RadioTherapyData> data = readData(patient);
            if (data.size() == 0) {
                LOGGER.warn(patient.patientId() + " had value " + hadRadiotherapy + " instead of No for field "
                        + FIELD_HADRADIOTHERAPY + " but radio therapy fields contained no values.");
                return Optional.empty();
            } else {
                return Optional.of(data);
            }
        }
    }

    @NotNull
    private List<RadioTherapyData> readData(@NotNull final EcrfPatient patient) {
        final List<String> endDates = GenericReader.getFieldValues(patient, FIELD_RADIOENDDATE);
        final List<String> sites = GenericReader.getFieldValues(patient, FIELD_RADIOSITE);
        final int maxLength = Utils.getMaxLength(Lists.newArrayList(endDates, sites),
                "Not all radio therapy fields contain the same number of values.");

        final List<RadioTherapyData> therapies = Lists.newArrayList();
        for (int index = 0; index < maxLength; index++) {
            final String site = Utils.getElemAtIndex(sites, index);
            final LocalDate endDate = Utils.getDate(Utils.getElemAtIndex(endDates, index), dateFormatter);
            if (endDate == null) {
                LOGGER.warn(FIELD_RADIOENDDATE + " did not contain valid date at index " + index + " for patient  "
                        + patient.patientId());
            }
            therapies.add(new RadioTherapyData(endDate, site));
        }
        return therapies;
    }

}
