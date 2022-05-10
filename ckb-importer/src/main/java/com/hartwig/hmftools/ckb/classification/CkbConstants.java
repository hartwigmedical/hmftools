package com.hartwig.hmftools.ckb.classification;

import java.util.Set;

import com.google.common.collect.Sets;

public final class CkbConstants {

    public static final String NO_GENE = "-";

    public static final String HRD_POSITIVE = "HRD pos";
    public static final String HRD_NEGATIVE = "HRD neg";
    public static final String MSI_HIGH = "MSI high";
    public static final String MSI_LOW = "MSI low";
    public static final String MSI_NEGATIVE = "MSI neg";
    public static final String TMB_HIGH = "TMB high";
    public static final String TMB_LOW = "TMB low";

    public static final Set<String> NON_EXISTING_GENES = Sets.newHashSet();

    public static final Set<String> UNMAPPABLE_GENES = Sets.newHashSet();

    private CkbConstants() {
    }

    static {
        NON_EXISTING_GENES.add("HRD");
        NON_EXISTING_GENES.add("MSI");
        NON_EXISTING_GENES.add("TMB");
        NON_EXISTING_GENES.add("Unknown");
        NON_EXISTING_GENES.add("TILs");

        UNMAPPABLE_GENES.add("CD24");
        UNMAPPABLE_GENES.add("EPOP");
    }
}
