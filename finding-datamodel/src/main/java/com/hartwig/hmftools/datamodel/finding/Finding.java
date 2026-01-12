package com.hartwig.hmftools.datamodel.finding;

import java.util.Comparator;

import org.jetbrains.annotations.NotNull;

public interface Finding
{
    Comparator<Finding> COMPARATOR = Comparator.comparing(Finding::findingKey);

    @NotNull
    String findingKey();
}
