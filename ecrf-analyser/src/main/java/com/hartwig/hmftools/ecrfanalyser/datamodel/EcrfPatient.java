package com.hartwig.hmftools.ecrfanalyser.datamodel;

import java.util.List;
import java.util.Map;
import java.util.Set;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class EcrfPatient {

    @NotNull
    private final String patientId;
    @NotNull
    private final Map<EcrfField, List<String>> fields;

    public EcrfPatient(@NotNull final String patientId, @NotNull final Map<EcrfField, List<String>> fields) {
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

    @Nullable
    public List<String> fieldValues(@NotNull EcrfField field) {
        return fields.get(field);
    }
}
