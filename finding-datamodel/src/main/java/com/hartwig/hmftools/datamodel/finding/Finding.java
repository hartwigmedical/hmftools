package com.hartwig.hmftools.datamodel.finding;

import java.util.Comparator;

import jakarta.validation.constraints.NotNull;

public interface Finding
{
    Comparator<Finding> COMPARATOR = Comparator.comparing(Finding::findingKey);

    @NotNull
    String findingKey();
}
