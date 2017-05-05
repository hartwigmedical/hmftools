package com.hartwig.hmftools.common.ecrf.datamodel;

import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Lists;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class EcrfPatient {

    @NotNull
    private final String patientId;
    @NotNull
    private final Map<EcrfField, List<String>> fields;
    @NotNull
    private final Map<String, List<EcrfStudyEvent>> studyEventsPerOID;

    public EcrfPatient(@NotNull final String patientId, @NotNull final Map<EcrfField, List<String>> fields,
            @NotNull final Map<String, List<EcrfStudyEvent>> studyEventsPerOID) {
        this.patientId = patientId;
        this.fields = fields;
        this.studyEventsPerOID = studyEventsPerOID;
    }

    @NotNull
    public String patientId() {
        return patientId;
    }

    @NotNull
    public Set<EcrfField> fields() {
        return fields.keySet();
    }

    @Nullable
    public List<String> fieldValuesByEcrfField(@NotNull final EcrfField field) {
        return fields.get(field);
    }

    @Nullable
    public List<String> fieldValuesByName(@NotNull final String fieldName) {
        for (final EcrfField field : fields.keySet()) {
            if (field.name().equals(fieldName)) {
                return fieldValuesByEcrfField(field);
            }
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
}
