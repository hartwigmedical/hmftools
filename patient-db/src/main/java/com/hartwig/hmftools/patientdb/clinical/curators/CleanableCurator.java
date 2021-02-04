package com.hartwig.hmftools.patientdb.clinical.curators;

import java.util.Set;

import org.jetbrains.annotations.NotNull;

public interface CleanableCurator {

    @NotNull
    Set<String> unusedSearchTerms();
}
