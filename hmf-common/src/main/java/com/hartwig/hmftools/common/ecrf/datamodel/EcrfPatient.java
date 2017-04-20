package com.hartwig.hmftools.common.ecrf.datamodel;

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
}
