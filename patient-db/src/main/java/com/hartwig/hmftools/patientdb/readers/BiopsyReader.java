package com.hartwig.hmftools.patientdb.readers;

import java.time.LocalDate;
import java.time.format.DateTimeFormatter;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.ecrf.datamodel.EcrfForm;
import com.hartwig.hmftools.common.ecrf.datamodel.EcrfItemGroup;
import com.hartwig.hmftools.common.ecrf.datamodel.EcrfPatient;
import com.hartwig.hmftools.common.ecrf.datamodel.EcrfStudyEvent;
import com.hartwig.hmftools.patientdb.data.BiopsyData;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

final class BiopsyReader {
    private static final Logger LOGGER = LogManager.getLogger(BiopsyReader.class);

    private static final String STUDY_BIOPSY = "SE.BIOPSY";
    private static final String FORM_BIOPS = "FRM.BIOPS";
    private static final String ITEMGROUP_BIOPSIES = "GRP.BIOPS.BIOPSIES";

    private static final String FIELD_DATE = "FLD.BIOPS.BIOPTDT";
    private static final String FIELD_LOCATION = "FLD.BIOPS.BILESSITE";
    private static final String FIELD_LOCATION_OTHER = "FLD.BIOPS.BIOTHLESSITE";

    private static final DateTimeFormatter DATE_FORMATTER = DateTimeFormatter.ofPattern("yyyy-MM-dd");

    private BiopsyReader() {
    }

    @NotNull
    static List<BiopsyData> read(@NotNull final EcrfPatient patient) {
        final List<BiopsyData> biopsies = Lists.newArrayList();
        for (final EcrfStudyEvent studyEvent : patient.studyEventsPerOID(STUDY_BIOPSY)) {
            for (final EcrfForm form : studyEvent.nonEmptyFormsPerOID(FORM_BIOPS, true)) {
                for (final EcrfItemGroup itemGroup : form.nonEmptyItemGroupsPerOID(ITEMGROUP_BIOPSIES, true)) {
                    final LocalDate date = itemGroup.readItemDate(FIELD_DATE, 0, DATE_FORMATTER, true);
                    final String location = itemGroup.readItemString(FIELD_LOCATION, 0, true);
                    if (location != null && location.trim().toLowerCase().startsWith("other")) {
                        final String location_other = itemGroup.readItemString(FIELD_LOCATION_OTHER, 0, true);
                        biopsies.add(new BiopsyData(date, location_other));
                    } else {
                        biopsies.add(new BiopsyData(date, location));
                    }
                }
            }
        }
        final Map<LocalDate, List<BiopsyData>> groups = biopsies.stream().filter(
                biopsy -> biopsy.date() != null).collect(Collectors.groupingBy(BiopsyData::date));
        for (final LocalDate biopsyDate : groups.keySet()) {
            final int biopsiesPerDate = groups.get(biopsyDate).size();
            if (biopsiesPerDate > 1) {
                LOGGER.warn(patient.patientId() + ": " + biopsiesPerDate + " biopsies with same date: " + biopsyDate);
            }
        }
        return biopsies;
    }
}
