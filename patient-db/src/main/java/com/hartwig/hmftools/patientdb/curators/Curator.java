package com.hartwig.hmftools.patientdb.curators;

import java.util.Set;

import org.jetbrains.annotations.NotNull;

public interface Curator {

    @NotNull
    Set<String> unusedSearchTerms();
}
