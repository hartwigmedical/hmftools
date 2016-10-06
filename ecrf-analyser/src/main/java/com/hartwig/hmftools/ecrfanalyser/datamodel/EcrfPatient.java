package com.hartwig.hmftools.ecrfanalyser.datamodel;

import java.util.Map;
import java.util.Set;

import org.jetbrains.annotations.NotNull;

public class EcrfPatient {

    @NotNull
    private final String patientId;
    @NotNull
    private final Map<EcrfField, String> fields;

    public EcrfPatient(@NotNull final String patientId, @NotNull final Map<EcrfField, String> fields) {
        this.patientId = patientId;
        this.fields = fields;
    }

    @NotNull
    public String patientId() {
        return patientId;
    }

    @NotNull
    public Set<EcrfField> fields() {
        return fields.keySet();
    }

    @NotNull
    public String fieldValue(@NotNull EcrfField field) {
        return fields.get(field);
    }
}
