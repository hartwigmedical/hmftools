package com.hartwig.hmftools.ckb.classification;

import java.util.Set;

import com.google.common.collect.Sets;

public final class CkbConstants {

    public static final String NO_GENE = "-";
    public static final Set<String> NON_EXISTING_GENES = Sets.newHashSet();

    // TODO Figure out how to get these genes into HMF v38 gene model
    public static final Set<String> UNMAPPABLE_GENES = Sets.newHashSet();

    // TODO Curate fusion legs in case the leg is not the gene on which the variant is defined.
    public static final Set<String> UNRESOLVABLE_FUSION_LEGS = Sets.newHashSet();

    // TODO Allow non-fusion events on genes that are not present in DNA directly.
    public static final Set<String> EXCLUSIVE_FUSION_GENES = Sets.newHashSet();

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

        UNRESOLVABLE_FUSION_LEGS.add("ACP3");
        UNRESOLVABLE_FUSION_LEGS.add("BCAR4");
        UNRESOLVABLE_FUSION_LEGS.add("CARS1");
        UNRESOLVABLE_FUSION_LEGS.add("CEP43");
        UNRESOLVABLE_FUSION_LEGS.add("COP1");
        UNRESOLVABLE_FUSION_LEGS.add("FYB1");
        UNRESOLVABLE_FUSION_LEGS.add("HMGN2P46");
        UNRESOLVABLE_FUSION_LEGS.add("LOC389473");
        UNRESOLVABLE_FUSION_LEGS.add("MALAT1");
        UNRESOLVABLE_FUSION_LEGS.add("MIR548F1");
        UNRESOLVABLE_FUSION_LEGS.add("MRTFB");
        UNRESOLVABLE_FUSION_LEGS.add("OGA");
        UNRESOLVABLE_FUSION_LEGS.add("PVT1");
        UNRESOLVABLE_FUSION_LEGS.add("RARS1");
        UNRESOLVABLE_FUSION_LEGS.add("RELCH");
        UNRESOLVABLE_FUSION_LEGS.add("RNF217-AS1");
        UNRESOLVABLE_FUSION_LEGS.add("SEPTIN7");
        UNRESOLVABLE_FUSION_LEGS.add("SEPTIN14");
        UNRESOLVABLE_FUSION_LEGS.add("SEPTIN6");
        UNRESOLVABLE_FUSION_LEGS.add("SEPTIN3");
        UNRESOLVABLE_FUSION_LEGS.add("SYN2");
        UNRESOLVABLE_FUSION_LEGS.add("ZFTA");
        UNRESOLVABLE_FUSION_LEGS.add("ZNF767P");

        EXCLUSIVE_FUSION_GENES.add("IGH");
        EXCLUSIVE_FUSION_GENES.add("IGK");
        EXCLUSIVE_FUSION_GENES.add("IGL");
    }
}
