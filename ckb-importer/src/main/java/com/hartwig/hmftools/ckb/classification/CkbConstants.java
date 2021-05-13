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
        UNMAPPABLE_GENES.add("HMGN2P46"); // Only appears as a leg in a fusion but not directly as a gene that holds a variant
        UNMAPPABLE_GENES.add("IGHV3-21");
        UNMAPPABLE_GENES.add("IGHV4-34");
        UNMAPPABLE_GENES.add("LOC389473"); // Only appears as a leg in a fusion but not directly as a gene that holds a variant
        UNMAPPABLE_GENES.add("MALAT1");
        UNMAPPABLE_GENES.add("MAP3K14");
        UNMAPPABLE_GENES.add("PVT1");
        UNMAPPABLE_GENES.add("RNF217-AS1");
        UNMAPPABLE_GENES.add("RYBP");
        UNMAPPABLE_GENES.add("SYN2");
        UNMAPPABLE_GENES.add("TERC");
        UNMAPPABLE_GENES.add("ZFTA");

        UNRESOLVABLE_FUSION_LEGS.add("ACP3"); // Map to ACPP
        UNRESOLVABLE_FUSION_LEGS.add("CARS1"); // Map to CARS
        UNRESOLVABLE_FUSION_LEGS.add("CEP43"); // Map to FGFR1OP
        UNRESOLVABLE_FUSION_LEGS.add("COP1"); // Map to RFWD2
        UNRESOLVABLE_FUSION_LEGS.add("FYB1"); // Map to FYB
        UNRESOLVABLE_FUSION_LEGS.add("MIR548F1"); // Gene exists in v38 but has no v37-v38 mapping
        UNRESOLVABLE_FUSION_LEGS.add("MRTFB"); // Map to MKL2
        UNRESOLVABLE_FUSION_LEGS.add("OGA"); // Map to MGEA5
        UNRESOLVABLE_FUSION_LEGS.add("RARS1"); // Map to RARS
        UNRESOLVABLE_FUSION_LEGS.add("RELCH"); // Map to KIAA1468
        UNRESOLVABLE_FUSION_LEGS.add("SEPTIN3"); // Map to SEPT3
        UNRESOLVABLE_FUSION_LEGS.add("SEPTIN6"); // Map to SEPT6
        UNRESOLVABLE_FUSION_LEGS.add("SEPTIN7"); // Map to SEPT7
        UNRESOLVABLE_FUSION_LEGS.add("SEPTIN14"); // Map to SEPT14
        UNRESOLVABLE_FUSION_LEGS.add("ZNF767P"); // Gene exists in v38 but has no v37-v38 mapping

        EXCLUSIVE_FUSION_GENES.add("IGH");
        EXCLUSIVE_FUSION_GENES.add("IGK");
        EXCLUSIVE_FUSION_GENES.add("IGL");
    }
}
