package com.hartwig.hmftools.patientdb.clinical.matchers;

import java.util.List;

import com.hartwig.hmftools.patientdb.clinical.ecrf.datamodel.ValidationFinding;

import org.jetbrains.annotations.NotNull;

public class MatchResult<T> {

    @NotNull
    private final List<T> values;
    @NotNull
    private final List<ValidationFinding> findings;

    public MatchResult(@NotNull final List<T> values, @NotNull final List<ValidationFinding> findings) {
        this.values = values;
        this.findings = findings;
    }

    @NotNull
    public List<T> values() {
        return values;
    }

    @NotNull
    public List<ValidationFinding> findings() {
        return findings;
    }
}
