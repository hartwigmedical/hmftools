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

        UNMAPPABLE_GENES.add("BCAR4");
        UNMAPPABLE_GENES.add("CD24");
        UNMAPPABLE_GENES.add("COX2");
        UNMAPPABLE_GENES.add("EPOP");
        UNMAPPABLE_GENES.add("FER1L5");
        UNMAPPABLE_GENES.add("IGHV3-21");
        UNMAPPABLE_GENES.add("IGHV4-34");
        UNMAPPABLE_GENES.add("MALAT1");
        UNMAPPABLE_GENES.add("MAP3K14");
        UNMAPPABLE_GENES.add("PVT1");
        UNMAPPABLE_GENES.add("RNF217-AS1");
        UNMAPPABLE_GENES.add("RYBP");
        UNMAPPABLE_GENES.add("SYN2");
        UNMAPPABLE_GENES.add("TERC");
        UNMAPPABLE_GENES.add("ZFTA");
    }
}
