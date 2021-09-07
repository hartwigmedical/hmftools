package com.hartwig.hmftools.patientdb.clinical.ecrf.datamodel;

import java.util.Collection;
import java.util.List;
import java.util.Map;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.Lists;
import com.google.common.collect.Multimap;

import org.jetbrains.annotations.NotNull;

public class EcrfPatient {

    @NotNull
    private final String patientId;
    @NotNull
    private final Map<String, List<EcrfStudyEvent>> studyEventsPerOID;
    @NotNull
    private final Multimap<String, EcrfDataField> fields;

    public EcrfPatient(@NotNull final String patientId, @NotNull final Map<String, List<EcrfStudyEvent>> studyEventsPerOID,
            @NotNull final List<EcrfDataField> fields) {
        this.patientId = patientId;
        this.studyEventsPerOID = studyEventsPerOID;
        this.fields = createDataFields(fields);
    }

    @NotNull
    public String patientId() {
        return patientId;
    }

    @NotNull
    public Collection<EcrfDataField> fields() {
        return fields.values();
    }

    @NotNull
    public Map<String, List<EcrfStudyEvent>> studyEventsPerOID() {
        return studyEventsPerOID;
    }

    @NotNull
    public List<EcrfStudyEvent> studyEventsPerOID(@NotNull String studyEventOID) {
        if (studyEventsPerOID().get(studyEventOID) == null) {
            return Lists.newArrayList();
        }
        return studyEventsPerOID.get(studyEventOID);
    }

    @NotNull
    private static Multimap<String, EcrfDataField> createDataFields(@NotNull List<EcrfDataField> ecrfDataFields) {
        Multimap<String, EcrfDataField> dataFieldMultimap = ArrayListMultimap.create();
        ecrfDataFields.forEach(field -> dataFieldMultimap.put(field.name(), field));
        return dataFieldMultimap;
    }
}
