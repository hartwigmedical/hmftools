package com.hartwig.hmftools.common.ecrf.datamodel;

import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.Lists;
import com.google.common.collect.Multimap;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

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

    @Nullable
    public List<String> fieldValuesByEcrfField(@NotNull final EcrfField field) {
        return fieldValuesByName(field.name());
    }

    @Nullable
    public List<String> fieldValuesByName(@NotNull final String fieldName) {
        final Collection<EcrfDataField> results = fields.get(fieldName);
        if (results != null && !results.isEmpty()) {
            return results.stream().map(EcrfDataField::itemValue).collect(Collectors.toList());
        }
        return null;
    }

    @NotNull
    public Map<String, List<EcrfStudyEvent>> studyEventsPerOID() {
        return studyEventsPerOID;
    }

    @NotNull
    public List<EcrfStudyEvent> studyEventsPerOID(@NotNull final String studyEventOID) {
        if (studyEventsPerOID().get(studyEventOID) == null) {
            return Lists.newArrayList();
        }
        return studyEventsPerOID.get(studyEventOID);
    }

    @NotNull
    private static Multimap<String, EcrfDataField> createDataFields(@NotNull final List<EcrfDataField> ecrfDataFields) {
        final Multimap<String, EcrfDataField> dataFieldMultimap = ArrayListMultimap.create();
        ecrfDataFields.forEach(field -> dataFieldMultimap.put(field.name(), field));
        return dataFieldMultimap;
    }
}
