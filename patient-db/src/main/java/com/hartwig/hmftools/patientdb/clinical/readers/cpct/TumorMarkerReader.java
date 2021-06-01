package com.hartwig.hmftools.patientdb.clinical.readers.cpct;

import java.time.LocalDate;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.patientdb.clinical.datamodel.ImmutableTumorMarkerData;
import com.hartwig.hmftools.patientdb.clinical.datamodel.TumorMarkerData;
import com.hartwig.hmftools.patientdb.clinical.ecrf.datamodel.EcrfForm;
import com.hartwig.hmftools.patientdb.clinical.ecrf.datamodel.EcrfItemGroup;
import com.hartwig.hmftools.patientdb.clinical.ecrf.datamodel.EcrfPatient;
import com.hartwig.hmftools.patientdb.clinical.ecrf.datamodel.EcrfStudyEvent;

import org.jetbrains.annotations.NotNull;

final class TumorMarkerReader {

    private static final String STUDY_TREATMENT = "SE.TREATMENT";
    private static final String FORM_RESPONSE = "FRM.RESPONSE";
    private static final String ITEMGROUP_RESPONSE = "GRP.RESPONSE.LB_TUMORMARKERS";

    private static final String FIELD_DATE = "FLD.RESPONSE.LBDTC_TUMORMARKERS";
    private static final String FIELD_MARKER = "FLD.RESPONSE.LBTERM_TUMORMARKERS";
    private static final String FIELD_MEASUREMENT = "FLD.RESPONSE.LBORRES_TUMORMARKERS";
    private static final String FIELD_UNIT = "FLD.RESPONSE.LBORRESU_TUMORMARKERS";

    private TumorMarkerReader() {
    }

    @NotNull
    static List<TumorMarkerData> read(@NotNull EcrfPatient patient) {
        List<TumorMarkerData> tumorMarkers = Lists.newArrayList();
        for (EcrfStudyEvent studyEvent : patient.studyEventsPerOID(STUDY_TREATMENT)) {
            for (EcrfForm form : studyEvent.nonEmptyFormsPerOID(FORM_RESPONSE)) {
                for (EcrfItemGroup responseGroup : form.nonEmptyItemGroupsPerOID(ITEMGROUP_RESPONSE)) {
                    LocalDate date = responseGroup.readItemDate(FIELD_DATE);
                    String marker = responseGroup.readItemString(FIELD_MARKER);
                    String measurement = responseGroup.readItemString(FIELD_MEASUREMENT);
                    String unit = responseGroup.readItemString(FIELD_UNIT);

                    tumorMarkers.add(ImmutableTumorMarkerData.of(patient.patientId(), date, marker, measurement, unit, form.status()));
                }
            }
        }

        return tumorMarkers;
    }
}
