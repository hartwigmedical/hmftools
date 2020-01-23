package com.hartwig.hmftools.protect.conclusion;

import org.jetbrains.annotations.NotNull;

public enum AberrationGenSummary {
    NOTE_SIGNATURE_OR_SV_ANNOTATION("NOTE SIGNATURE/SV ANNOTATION"),
    LOW_PURITY("LOW PURITY"),
    NOT_DETECTABLE_TUMOR_PURITY("NOT DETECTABLE TUMOR PURITY"),
    NONE_FOUND("NONE FOUND"),
    HIGH_TMB_PASSENGER_VS_DRIVER("HIGH TMB - PASSENGER vs DRIVER"),
    PARTIAL_LOSS("PARTIAL LOSS"),
    HOMOZYGOUS_DISRUPTION("HOMOZYGOUS DISRUPTION"),
    DRIVER_VUS("DRIVER VUS"),
    NOT_BIALLELIC_TSG("NOT BIALLELIC (TSG)"),
    SPLICE_SITE_DAMAGING("SPLICE SITE DAMAGING"),
    GERMLINE_NOTIFICATION("GERMLINE NOTIFICATION"),
    SOMATIC_REMARK("SOMATIC REMARK"),
    AKT1_MUTATION("AKT1 mutation"),
    ALK_FUSION("ALK fusion"),
    ALK_FUSION_AND_ALK_RESISTANCE_MUTATION("ALK fusion + ALK resistance mutation"),
    APC_INACTIVATION("APC inactivation"),
    AR_AMPLIFICATION("AR amplification"),
    AR_MUTATION("AR mutation"),
    ARID1A_MUTATION("ARID1A mutation"),
    ACVR2A_INACTIVATION_OR_LOSS_CRC("ACVR2A inactivation/loss (CRC)"),
    ATM_INACTIVATION_OR_LOSS("ATM inactivation/loss"),
    AXL_AMPLIFICATION("AXL amplification"),
    B2M_INACTIVATION("B2M inactivation"),
    BAP1_INACTIVATION_OR_LOSS("BAP1 inactivation/loss"),
    BAP1_LOSS_MESOTHELIOMA("BAP1 loss (mesothelioma)"),
    BRAF_ACTIVATING_MUTATION("BRAF activating mutation"),
    BRAF_ACTIVATING_MUTATION_CRC("BRAF activating mutation (CRC)"),
    BRAF_ACTIVATING_MUTATION_EXON_11("BRAF activating mutation (exon 11)"),
    BRAF_CODON_594_MUTATION("BRAF codon 594 mutation"),
    BRAF_FUSION("BRAF fusion"),
    BRAF_INTERNAL_DELETION("BRAF internal deletion (exon X-X)"),
    BRCA1_INACTIVATION_OR_LOSS("BRCA1 inactivation/loss"),
    BRCA2_INACTIVATION_OR_LOSS("BRCA2 inactivation/loss"),
    CCND1_AMPLIFICATION("CCND1 amplification"),
    CCND2_AMPLIFICATION("CCND2 amplification"),
    CCNE1_AMPLIFICATION("CCNE1 amplification"),
    CDK12_INACTIVATION_OR_LOSS("CDK12 inactivation/loss"),
    CDK12_INACTIVATION_PROSTATE("CDK12 inactivation (prostate)"),
    CDK4_AMPLIFICATION("CDK4 amplification"),
    CDKN1A_INACTIVATION_OR_LOSS("CDKN1A inactivation/loss"),
    CDKN2A_INACTIVATION_OR_LOSS("CDKN2A inactivation/loss"),
    CHEK2_INACTIVATION_OR_LOSS("CHEK2 inactivation/loss"),
    CTNNB1_MUTATION_HOTSPOT("CTNNB1 mutation (hotspot)"),
    DDR2_AMPLIFICATION("DDR2 amplification"),
    EGFR_ACTIVATING_MUTATION("EGFR activating mutation"),
    EGFR_S768I("EGFR S768I"),
    EGFR_AMPLIFICATION("EGFR amplification"),
    EGFR_AMPLIFICATION_LUNG_AND_WILDTYPE("EGFR amplification (lung & wildtype EGFR)"),
    EGFR_EXON_20_MUTATION("EGFR exon 20 mutation (761-775)"),
    EGFR_RESISTANCE_MUTATION("EGFR resistance mutation"),
    EGFR_FUSION("EGFR fusion"),
    ERBB2_ACTIVATING_MUTATION_HOTSPOT("ERBB2 activating mutation (hotspot)"),
    ERBB2_ACTIVATING_MUTATION_EXON_20("ERBB2 activating mutation (exon 20)"),
    ERBB2_AMPLIFICATION("ERBB2 amplification"),
    ERBB3_AMPLIFICATION("ERBB3 amplification"),
    ERBB4_ACTIVATING_MUTATION("ERBB4 activating mutation"),
    ESR1_ACTIVATING_MUTATION("ESR1 activating mutation"),
    EZH2_ACTIVATING_MUTATION("EZH2 activating mutation"),
    EZH2_AMPLIFICATION("EZH2 amplification"),
    FAT1_INACTIVATION("FAT1 inactivation"),
    FBXW7_MUTATION("FBXW7 mutation"),
    FGF3_AMPLIFICATION("FGF3 amplification"),
    FGFR1_AMPLIFICATION("FGFR1 amplification"),
    FGFR1_FUSION("FGFR1 fusion"),
    FGFR2_FUSION("FGFR2 fusion"),
    FGFR2_ACTIVATING_MUTATION("FGFR2 activating mutation"),
    FGFR3_ACTIVATING_MUTATION("FGFR3 activating mutation"),
    FGFR3_FUSION("FGFR3 fusion"),
    FGFR4_FUSION("FGFR4 fusion"),
    FLT3_ACTIVATING_MUTATION("FLT3 activating mutation"),
    FLT3_AMPLIFICATION("FLT3 amplification"),
    GATA3_INACTIVATION("GATA3 inactivation"),
    HIGH_MTL("High TMB"),
    High_TMB_LUNG_MELANOMOA("High TMB (lung, melanoma)"),
    HR_DEFICIENT("HR-deficient"),
    HRAS_ACTIVATING_MUTATION_OTHER("HRAS activating mutation (other)"),
    IDH1_OR_2_MUTATION_HOTSPOT("IDH1/2 mutation (hotspot)"),
    JAK1_INACTIVATION_OR_LOSS("JAK1 inactivation/loss"),
    JAK2_INACTIVATION_OR_LOSS("JAK2 inactivation/loss"),
    JAK2_AMPLIFICATION("JAK2 amplification"),
    KDR_AMPLIFICATION("KDR amplification"),
    KEAP1_MUTATION("KEAP1 mutation"),
    KIT_ACTIVATING_MUTATION("KIT activating mutation"),
    KIT_AMPLIFICATION("KIT amplification"),
    KRAS_ACTIVATION_MUTATION_GI_TRACT("KRAS activating mutation (GI-tract)"),
    KRAS_ACTIVATING_MUTATION_LUNG("KRAS activating mutation (lung)"),
    KRAS_p_Gly12Cys_LUNG("KRAS p.Gly12Cys (lung)"),
    KRAS_ACTIVATING_MUTATION_OTHER("KRAS activating mutation (other)"),
    KRAS_AMPLIFICATION("KRAS amplification"),
    MACROD2_LOSS("MACROD2 loss"),
    MAP2K1_ACTIVATING_MUTATION("MAP2K1 activating mutation"),
    MAP2K4_INACTIVATING_OR_LOSS("MAP2K4 inactivation/loss"),
    MAP3K1_INACTIVATING_OR_LOSS("MAP3K1 inactivation/loss"),
    MCL1_AMPLIFICATION("MCL1 amplification"),
    MDM2_AMPLIFICATION("MDM2 amplification"),
    MET_AMPLIFICATION("MET amplification"),
    MET_EXON_14_SPLICING("MET exon 14 splicing"),
    MET_FUSION("MET fusion"),
    MEN1_MUTATION_NET("MEN1 mutation (NET)"),
    MSH2_INACTIVATION("MSH2 inactivation"),
    MSI("MSI"),
    MYC_AMPLIFICATION_BREAST_CANCER("MYC amplification (breast cancer)"),
    MYD88_ACTIVATING_MUTATION_L265P("MYD88 activating mutation (L265P)"),
    NF1_INACTIVATION_OR_MUTATION("NF1 inactivation/mutation"),
    NOTCH1("NOTCH1"),
    NRAS_ACTIVATING_MUTATION_CRC("NRAS activating mutation (CRC)"),
    NRAS_ACTIVATING_MUTATION_OTHER("NRAS activating mutation (other)"),
    NRAS_ACTIVATING_MUTATION_LUNG("NRAS activating mutation (lung)"),
    NRG1_ACTIVATING_MUTATION("NRG1 activating mutation"),
    NRG1_FUSION("NRG1 fusion"),
    NTRK1_AMPLIFICATION("NTRK1 amplfication"),
    NTRK2_FUSION("NTRK2 fusion"),
    NTRK3_FUSION("NTRK3 fusion"),
    NTRK3_ACTIVATION_MUTATION("NTRK3 activating mutation"),
    PDGFRA_AMPLIFICATION("PDGFRA amplification"),
    PIK3CA_ACTIVATING_MUTATION("PIK3CA activating mutation"),
    POLE_INACTIVATING_MUTATION("POLE inactivating mutation"),
    PTCH1_INACTIVATING_MUTATION("PTCH1 inactivating mutation"),
    PTEN_INACTIVATION_OR_LOSS("PTEN inactivation/loss"),
    RAC1_p_Pro29Ser_MELANOMA("RAC1 p.Pro29Ser (melanoma)"),
    RAD51B_LOSS("RAD51B loss"),
    RASA1_INACTIVATION("RASA1 inactivation"),
    RAF1_ACTIVATING_MUTATION("RAF1 activating mutation"),
    RAF1_AMPLIFICATION("RAF1 amplification"),
    RB1_AND_TP53_INACTIVATION_OR_LOSS_NSCLC("RB1 and TP53 inactivation/loss (NSCLC)"),
    RET_ACTIVATING_MUTATION("RET activating mutation"),
    RET_AMPLIFICATION("RET amplification"),
    RET_FUSION("RET fusion"),
    RNF43_INACTIVATION_MUTATION("RNF43 inactivating mutation"),
    ROS1_FUSION("ROS1 fusion"),
    RSF1_AMPLIFICATION("RSF1 amplification"),
    RSF1_AMPLIFICATION_BREAST_CANCER("RSF1 amplification (breast cancer)"),
    RSPO2_FUSION("RSPO2 fusion"),
    RSPO3_FUSION("RSPO3 fusion"),
    SMAD2_INACTIVATION_OR_LOSS("SMAD2 inactivation/loss"),
    SMAD3_INACTIVATION_OR_LOSS("SMAD3 inactivation/loss"),
    SMAD4_INACTIVATION_OR_LOSS("SMAD4 inactivation/loss"),
    SMAD4_INACTIVATION_OR_LOSS_CRC("SMAD4 inactivation/loss (CRC)"),
    SMARCA4_INACTIVATING_MUTATION("SMARCA4 inactivating mutation"),
    SMARCA4_LOSS("SMARCA4 loss"),
    STK11_INACTIVATION_OR_LOSS("STK11 inactivation/loss"),
    TERT_PROMOTOR_MUTATION("TERT promoter mutation"),
    TGFBR2_INACTIVATION_OR_LOSS("TGFBR2 inactivation/loss"),
    TMPRSS2_ERG_PROSTATE("TMPRSS2-ERG (prostate)"),
    TOP1_AMPLIFICATION("TOP1 amplification"),
    TOP2A_BREAST_CANCER("TOP2A (breast cancer)"),
    TP53_INACTIVATION("TP53 inactivation"),
    TP53_LOSS("TP53 loss"),
    TP53_SPLICE_SITE("TP53 splice site"),
    TYMS_AMPLIFICATION("TYMS amplification"),
    ZFP36L1_INACTIVATION("ZFP36L1 inactivation"),
    ZFP36L2_INACTIVATION("ZFP36L2 inactivation"),
    ZNF703_AMPLIFICATION("ZNF703 amplification");

    @NotNull
    private final String readableString;

    AberrationGenSummary(@NotNull final String readableString) {
        this.readableString = readableString;

    }

    @NotNull
    public String readableString() {
        return readableString;
    }

    @NotNull
    public static AberrationGenSummary fromString(@NotNull String event) {
        switch (event) {
            case "NOTE SIGNATURE/SV ANNOTATION":
                return NOTE_SIGNATURE_OR_SV_ANNOTATION;
            case "LOW PURITY":
                return LOW_PURITY;
            case "NOT DETECTABLE TUMOR PURITY":
                return NOT_DETECTABLE_TUMOR_PURITY;
            case "NONE FOUND":
                return NONE_FOUND;
            case "HIGH TMB - PASSENGER vs DRIVER":
                return HIGH_TMB_PASSENGER_VS_DRIVER;
            case "PARTIAL LOSS":
                return PARTIAL_LOSS;
            case "HOMOZYGOUS DISRUPTION":
                return HOMOZYGOUS_DISRUPTION;
            case "DRIVER VUS":
                return DRIVER_VUS;
            case "NOT BIALLELIC (TSG)":
                return NOT_BIALLELIC_TSG;
            case "SPLICE SITE DAMAGING":
                return SPLICE_SITE_DAMAGING;
            case "GERMLINE NOTIFICATION":
                return GERMLINE_NOTIFICATION;
            case "SOMATIC REMARK":
                return SOMATIC_REMARK;
            case "AKT1 mutation":
                return AKT1_MUTATION;
            case "ALK fusion":
                return ALK_FUSION;
            case "ALK fusion + ALK resistance mutation":
                return ALK_FUSION_AND_ALK_RESISTANCE_MUTATION;
            case "APC inactivation":
                return APC_INACTIVATION;
            case "AR amplification":
                return AR_AMPLIFICATION;
            case "AR mutation":
                return AR_MUTATION;
            case "ARID1A mutation":
                return ARID1A_MUTATION;
            case "ACVR2A inactivation/loss (CRC)":
                return ACVR2A_INACTIVATION_OR_LOSS_CRC;
            case "ATM inactivation/loss":
                return ATM_INACTIVATION_OR_LOSS;
            case "AXL amplification":
                return AXL_AMPLIFICATION;
            case "B2M inactivation":
                return B2M_INACTIVATION;
            case "BAP1 inactivation/loss":
                return BAP1_INACTIVATION_OR_LOSS;
            case "BAP1 loss (mesothelioma)":
                return BAP1_LOSS_MESOTHELIOMA;
            case "BRAF activating mutation":
                return BRAF_ACTIVATING_MUTATION;
            case "BRAF activating mutation (CRC)":
                return BRAF_ACTIVATING_MUTATION_CRC;
            case "BRAF activating mutation (exon 11)":
                return BRAF_ACTIVATING_MUTATION_EXON_11;
            case "BRAF codon 594 mutation":
                return BRAF_CODON_594_MUTATION;
            case "BRAF fusion":
                return BRAF_FUSION;
            case "BRAF internal deletion (exon X-X)":
                return BRAF_INTERNAL_DELETION;
            case "BRCA1 inactivation/loss":
                return BRCA1_INACTIVATION_OR_LOSS;
            case "BRCA2 inactivation/loss":
                return BRCA2_INACTIVATION_OR_LOSS;
            case "CCND1 amplification":
                return CCND1_AMPLIFICATION;
            case "CCND2 amplification":
                return CCND2_AMPLIFICATION;
            case "CCNE1 amplification":
                return CCNE1_AMPLIFICATION;
            case "CDK12 inactivation/loss":
                return CDK12_INACTIVATION_OR_LOSS;
            case "CDK12 inactivation (prostate)":
                return CDK12_INACTIVATION_PROSTATE;
            case "CDK4 amplification":
                return CDK4_AMPLIFICATION;
            case "CDKN1A inactivation/loss":
                return CDKN1A_INACTIVATION_OR_LOSS;
            case "CDKN2A inactivation/loss":
                return CDKN2A_INACTIVATION_OR_LOSS;
            case "CHEK2 inactivation/loss":
                return CHEK2_INACTIVATION_OR_LOSS;
            case "CTNNB1 mutation (hotspot)":
                return CTNNB1_MUTATION_HOTSPOT;
            case "DDR2 amplification":
                return DDR2_AMPLIFICATION;
            case "EGFR activating mutation":
                return EGFR_ACTIVATING_MUTATION;
            case "EGFR S768I":
                return EGFR_S768I;
            case "EGFR amplification":
                return EGFR_AMPLIFICATION;
            case "EGFR amplification (lung & wildtype EGFR)":
                return EGFR_AMPLIFICATION_LUNG_AND_WILDTYPE;
            case "EGFR exon 20 mutation (761-775)":
                return EGFR_EXON_20_MUTATION;
            case "EGFR resistance mutation":
                return EGFR_RESISTANCE_MUTATION;
            case "EGFR fusion":
                return EGFR_FUSION;
            case "ERBB2 activating mutation (hotspot)":
                return ERBB2_ACTIVATING_MUTATION_HOTSPOT;
            case "ERBB2 activating mutation (exon 20)":
                return ERBB2_ACTIVATING_MUTATION_EXON_20;
            case "ERBB2 amplification":
                return ERBB2_AMPLIFICATION;
            case "ERBB3 amplification":
                return ERBB3_AMPLIFICATION;
            case "ERBB4 activating mutation":
                return ERBB4_ACTIVATING_MUTATION;
            case "ESR1 activating mutation":
                return ESR1_ACTIVATING_MUTATION;
            case "EZH2 activating mutation":
                return EZH2_ACTIVATING_MUTATION;
            case "EZH2 amplification":
                return EZH2_AMPLIFICATION;
            case "FAT1 inactivation":
                return FAT1_INACTIVATION;
            case "FBXW7 mutation":
                return FBXW7_MUTATION;
            case "FGF3 amplification":
                return FGF3_AMPLIFICATION;
            case "FGFR1 amplification":
                return FGFR1_AMPLIFICATION;
            case "FGFR1 fusion":
                return FGFR1_FUSION;
            case "FGFR2 fusion":
                return FGFR2_FUSION;
            case "FGFR2 activating mutation":
                return FGFR2_ACTIVATING_MUTATION;
            case "FGFR3 activating mutation":
                return FGFR3_ACTIVATING_MUTATION;
            case "FGFR3 fusion":
                return FGFR3_FUSION;
            case "FGFR4 fusion":
                return FGFR4_FUSION;
            case "FLT3 activating mutation":
                return FLT3_ACTIVATING_MUTATION;
            case "FLT3 amplification":
                return FLT3_AMPLIFICATION;
            case "GATA3 inactivation":
                return GATA3_INACTIVATION;
            case "High TMB":
                return HIGH_MTL;
            case "High TMB (lung, melanoma)":
                return High_TMB_LUNG_MELANOMOA;
            case "HR-deficient":
                return HR_DEFICIENT;
            case "HRAS activating mutation (other)":
                return HRAS_ACTIVATING_MUTATION_OTHER;
            case "IDH1/2 mutation (hotspot)":
                return IDH1_OR_2_MUTATION_HOTSPOT;
            case "JAK1 inactivation/loss":
                return JAK1_INACTIVATION_OR_LOSS;
            case "JAK2 inactivation/loss":
                return JAK2_INACTIVATION_OR_LOSS;
            case "JAK2 amplification":
                return JAK2_AMPLIFICATION;
            case "KDR amplification":
                return KDR_AMPLIFICATION;
            case "KEAP1 mutation":
                return KEAP1_MUTATION;
            case "KIT activating mutation":
                return KIT_ACTIVATING_MUTATION;
            case "KIT amplification":
                return KIT_AMPLIFICATION;
            case "KRAS activating mutation (GI-tract)":
                return KRAS_ACTIVATION_MUTATION_GI_TRACT;
            case "KRAS activating mutation (lung)":
                return KRAS_ACTIVATING_MUTATION_LUNG;
            case "KRAS p.Gly12Cys (lung)":
                return KRAS_p_Gly12Cys_LUNG;
            case "KRAS activating mutation (other)":
                return KRAS_ACTIVATING_MUTATION_OTHER;
            case "KRAS amplification":
                return KRAS_AMPLIFICATION;
            case "MACROD2 loss":
                return MACROD2_LOSS;
            case "MAP2K1 activating mutation":
                return MAP2K1_ACTIVATING_MUTATION;
            case "MAP2K4 inactivation/loss":
                return MAP2K4_INACTIVATING_OR_LOSS;
            case "MAP3K1 inactivation/loss":
                return MAP3K1_INACTIVATING_OR_LOSS;
            case "MCL1 amplification":
                return MCL1_AMPLIFICATION;
            case "MDM2 amplification":
                return MDM2_AMPLIFICATION;
            case "MET amplification":
                return MET_AMPLIFICATION;
            case "MET exon 14 splicing":
                return MET_EXON_14_SPLICING;
            case "MET fusion":
                return MET_FUSION;
            case "MEN1 mutation (NET)":
                return MEN1_MUTATION_NET;
            case "MSH2 inactivation":
                return MSH2_INACTIVATION;
            case "MSI":
                return MSI;
            case "MYC amplification (breast cancer)":
                return MYC_AMPLIFICATION_BREAST_CANCER;
            case "MYD88 activating mutation (L265P)":
                return MYD88_ACTIVATING_MUTATION_L265P;
            case "NF1 inactivation/mutation":
                return NF1_INACTIVATION_OR_MUTATION;
            case "NOTCH1":
                return NOTCH1;
            case "NRAS activating mutation (CRC)":
                return NRAS_ACTIVATING_MUTATION_CRC;
            case "NRAS activating mutation (other)":
                return NRAS_ACTIVATING_MUTATION_OTHER;
            case "NRAS activating mutation (lung)":
                return NRAS_ACTIVATING_MUTATION_LUNG;
            case "NRG1 activating mutation":
                return NRG1_ACTIVATING_MUTATION;
            case "NRG1 fusion":
                return NRG1_FUSION;
            case "NTRK1 amplfication":
                return NTRK1_AMPLIFICATION;
            case "NTRK2 fusion":
                return NTRK2_FUSION;
            case "NTRK3 fusion":
                return NTRK3_FUSION;
            case "NTRK3 activating mutation":
                return NTRK3_ACTIVATION_MUTATION;
            case "PDGFRA amplification":
                return PDGFRA_AMPLIFICATION;
            case "PIK3CA activating mutation":
                return PIK3CA_ACTIVATING_MUTATION;
            case "POLE inactivating mutation":
                return POLE_INACTIVATING_MUTATION;
            case "PTCH1 inactivating mutation":
                return PTCH1_INACTIVATING_MUTATION;
            case "PTEN inactivation/loss":
                return PTEN_INACTIVATION_OR_LOSS;
            case "RAC1 p.Pro29Ser (melanoma)":
                return RAC1_p_Pro29Ser_MELANOMA;
            case "RAD51B loss":
                return RAD51B_LOSS;
            case "RASA1 inactivation":
                return RASA1_INACTIVATION;
            case "RAF1 activating mutation":
                return RAF1_ACTIVATING_MUTATION;
            case "RAF1 amplification":
                return RAF1_AMPLIFICATION;
            case "RB1 and TP53 inactivation/loss (NSCLC)":
                return RB1_AND_TP53_INACTIVATION_OR_LOSS_NSCLC;
            case "RET activating mutation":
                return RET_ACTIVATING_MUTATION;
            case "RET amplification":
                return RET_AMPLIFICATION;
            case "RET fusion":
                return RET_FUSION;
            case "RNF43 inactivating mutation":
                return RNF43_INACTIVATION_MUTATION;
            case "ROS1 fusion":
                return ROS1_FUSION;
            case "RSF1 amplification":
                return RSF1_AMPLIFICATION;
            case "RSF1 amplification (breast cancer)":
                return RSF1_AMPLIFICATION_BREAST_CANCER;
            case "RSPO2 fusion":
                return RSPO2_FUSION;
            case "RSPO3 fusion":
                return RSPO3_FUSION;
            case "SMAD2 inactivation/loss":
                return SMAD2_INACTIVATION_OR_LOSS;
            case "SMAD3 inactivation/loss":
                return SMAD3_INACTIVATION_OR_LOSS;
            case "SMAD4 inactivation/loss":
                return SMAD4_INACTIVATION_OR_LOSS;
            case "SMAD4 inactivation/loss (CRC)":
                return SMAD4_INACTIVATION_OR_LOSS_CRC;
            case "SMARCA4 inactivating mutation":
                return SMARCA4_INACTIVATING_MUTATION;
            case "SMARCA4 loss":
                return SMARCA4_LOSS;
            case "STK11 inactivation/loss":
                return STK11_INACTIVATION_OR_LOSS;
            case "TERT promoter mutation":
                return TERT_PROMOTOR_MUTATION;
            case "TGFBR2 inactivation/loss":
                return TGFBR2_INACTIVATION_OR_LOSS;
            case "TMPRSS2-ERG (prostate)":
                return TMPRSS2_ERG_PROSTATE;
            case "TOP1 amplification":
                return TOP1_AMPLIFICATION;
            case "TOP2A (breast cancer)":
                return TOP2A_BREAST_CANCER;
            case "TP53 inactivation":
                return TP53_INACTIVATION;
            case "TP53 loss":
                return TP53_LOSS;
            case "TP53 splice site":
                return TP53_SPLICE_SITE;
            case "TYMS amplification":
                return TYMS_AMPLIFICATION;
            case "ZFP36L1 inactivation":
                return ZFP36L1_INACTIVATION;
            case "ZFP36L2 inactivation":
                return ZFP36L2_INACTIVATION;
            case "ZNF703 amplification":
                return ZNF703_AMPLIFICATION;
            default:
                throw new IllegalArgumentException("Unrecognized gene event " + event);
        }
    }
}
