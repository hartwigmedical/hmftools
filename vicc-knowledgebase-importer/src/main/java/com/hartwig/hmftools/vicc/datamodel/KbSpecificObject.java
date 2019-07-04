package com.hartwig.hmftools.vicc.datamodel;

import org.jetbrains.annotations.Nullable;

public interface KbSpecificObject {

    @Nullable
    Cgi cgi();

    @Nullable
    BRCA brca();

    @Nullable
    Sage sage();

    @Nullable
    Pmkb pmkb();

    @Nullable
    Oncokb oncoKb();

    @Nullable
    Jax jax();

    @Nullable
    JaxTrials jaxTrials();

    @Nullable
    MolecularMatch molecularMatch();

    @Nullable
    MolecularMatchTrials molecularMatchTrials();

    @Nullable
    Civic civic();

}
