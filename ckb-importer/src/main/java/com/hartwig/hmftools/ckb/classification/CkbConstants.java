package com.hartwig.hmftools.ckb.classification;

import java.util.Set;

import com.google.common.collect.Sets;

public final class CkbConstants {

    public static final String NO_GENE = "-";

    public static final Set<String> NON_EXISTING_GENES = Sets.newHashSet();

    public static final Set<String> UNMAPPABLE_GENES = Sets.newHashSet();

    private CkbConstants() {
    }

    static {
        NON_EXISTING_GENES.add("MSI");
        NON_EXISTING_GENES.add("TMB");
        NON_EXISTING_GENES.add("Unknown");
        NON_EXISTING_GENES.add("TILs");

        UNMAPPABLE_GENES.add("CD24");
        UNMAPPABLE_GENES.add("EPOP");
    }
}
