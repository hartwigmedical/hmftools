package com.hartwig.hmftools.patientdb.readers;

import java.time.LocalDate;
import java.time.format.DateTimeFormatter;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.ecrf.datamodel.EcrfForm;
import com.hartwig.hmftools.common.ecrf.datamodel.EcrfItemGroup;
import com.hartwig.hmftools.common.ecrf.datamodel.EcrfPatient;
import com.hartwig.hmftools.common.ecrf.datamodel.EcrfStudyEvent;
import com.hartwig.hmftools.patientdb.data.BiopsyClinicalData;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class BiopsyClinicalDataReader {
    private static final Logger LOGGER = LogManager.getLogger(BiopsyClinicalDataReader.class);

    private static final String STUDY_BIOPSY = "SE.BIOPSY";
    private static final String FORM_BIOPS = "FRM.BIOPS";
    private static final String ITEMGROUP_BIOPSIES = "GRP.BIOPS.BIOPSIES";

    private static final String FIELD_DATE = "FLD.BIOPS.BIOPTDT";
    private static final String FIELD_LOCATION = "FLD.BIOPS.BILESSITE";
    //    private static final String FIELD_LOCATION_OTHER = "FLD.BIOPS.BIOTHLESSITE";

    private static final DateTimeFormatter dateFormatter = DateTimeFormatter.ofPattern("yyyy-MM-dd");

    @NotNull
    public static List<BiopsyClinicalData> read(@NotNull final EcrfPatient patient) {
        final List<BiopsyClinicalData> biopsies = Lists.newArrayList();
        for (final EcrfStudyEvent studyEvent : patient.studyEventsPerOID(STUDY_BIOPSY)) {
            for (final EcrfForm form : studyEvent.formsPerOID(FORM_BIOPS)) {
                if (form.isEmpty()) {
                    LOGGER.warn("Ignoring empty form: " + FORM_BIOPS + " for patient " + patient.patientId());
                    continue;
                }
                for (final EcrfItemGroup itemGroup : form.itemGroupsPerOID(ITEMGROUP_BIOPSIES)) {
                    if (itemGroup.isEmpty()) {
                        LOGGER.warn("Ignoring empty item group: " + ITEMGROUP_BIOPSIES + " for patient "
                                + patient.patientId());
                        continue;
                    }
                    final LocalDate date = itemGroup.readItemDate(FIELD_DATE, 0, dateFormatter);
                    final String location = itemGroup.readItemString(FIELD_LOCATION, 0);
                    if (date == null) {
                        LOGGER.warn("Found biopsy with empty date for patient: " + patient.patientId());
                    }
                    biopsies.add(new BiopsyClinicalData(date, location, null));
                }
            }
        }
        return biopsies;
    }
}
