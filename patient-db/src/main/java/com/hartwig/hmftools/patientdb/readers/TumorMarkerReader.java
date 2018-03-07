package com.hartwig.hmftools.patientdb.readers;

import java.time.LocalDate;
import java.time.format.DateTimeFormatter;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.ecrf.datamodel.EcrfForm;
import com.hartwig.hmftools.common.ecrf.datamodel.EcrfItemGroup;
import com.hartwig.hmftools.common.ecrf.datamodel.EcrfPatient;
import com.hartwig.hmftools.common.ecrf.datamodel.EcrfStudyEvent;
import com.hartwig.hmftools.patientdb.data.ImmutableTumorMarkerData;
import com.hartwig.hmftools.patientdb.data.TumorMarkerData;

import org.jetbrains.annotations.NotNull;

final class TumorMarkerReader {
    private static final String STUDY_TREATMENT = "SE.TREATMENT";
    private static final String FORM_RESPONSE = "FRM.RESPONSE";
    private static final String ITEMGROUP_RESPONSE = "GRP.RESPONSE.RESPONSE";

    private static final String FIELD_DATE = "FLD.RESPONSE.LBDTC_TUMORMARKERS";
    private static final String FIELD_MARKER = "FLD.RESPONSE.LBTERM_TUMORMARKERS";
    private static final String FIELD_MEASUREMENT = "FLD.RESPONSE.LBORRES_TUMORMARKERS";
    private static final String FIELD_UNIT = "FLD.RESPONSE.LBORRESU_TUMORMARKERS";

    private static final DateTimeFormatter DATE_FORMATTER = DateTimeFormatter.ofPattern("yyyy-MM-dd");

    private TumorMarkerReader() {
    }

    @NotNull
    static List<TumorMarkerData> read(@NotNull final EcrfPatient patient) {
        final List<TumorMarkerData> tumorMarkers = Lists.newArrayList();
        for (final EcrfStudyEvent studyEvent : patient.studyEventsPerOID(STUDY_TREATMENT)) {
            for (final EcrfForm form : studyEvent.nonEmptyFormsPerOID(FORM_RESPONSE, false)) {
                for (final EcrfItemGroup responseGroup : form.nonEmptyItemGroupsPerOID(ITEMGROUP_RESPONSE, false)) {
                    LocalDate date = responseGroup.readItemDate(FIELD_DATE, 0, DATE_FORMATTER, false);
                    String marker = responseGroup.readItemString(FIELD_MARKER, 0, false);
                    String measurement = responseGroup.readItemString(FIELD_MEASUREMENT, 0, false);
                    String unit = responseGroup.readItemString(FIELD_UNIT, 0, false);
                    tumorMarkers.add(ImmutableTumorMarkerData.of(patient.patientId(),
                            date,
                            marker,
                            measurement,
                            unit,
                            form.status(),
                            form.locked()));
                }
            }
        }
        return tumorMarkers;
    }
}
