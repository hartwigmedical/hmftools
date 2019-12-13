package com.hartwig.hmftools.vicc.reader;

import static com.hartwig.hmftools.vicc.reader.JsonFunctions.string;

import java.util.List;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.gson.JsonObject;
import com.hartwig.hmftools.vicc.datamodel.brca.BRCA;
import com.hartwig.hmftools.vicc.datamodel.brca.BRCApart1;
import com.hartwig.hmftools.vicc.datamodel.brca.BRCApart2;
import com.hartwig.hmftools.vicc.datamodel.brca.ImmutableBRCA;
import com.hartwig.hmftools.vicc.datamodel.brca.ImmutableBRCApart1;
import com.hartwig.hmftools.vicc.datamodel.brca.ImmutableBRCApart2;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

final class BRCAObjectFactory {

    private static final Logger LOGGER = LogManager.getLogger(BRCAObjectFactory.class);

    private static final List<Integer> EXPECTED_BRCA_ELEMENT_SIZES = Lists.newArrayList(137);

    private BRCAObjectFactory() {
    }

    @NotNull
    static BRCA create(@NotNull JsonObject brcaObject) {
        Set<String> keysBRCA = brcaObject.keySet();

        if (!EXPECTED_BRCA_ELEMENT_SIZES.contains(keysBRCA.size())) {
            LOGGER.warn("Found {} in brca rather than the expected {}", keysBRCA.size(), EXPECTED_BRCA_ELEMENT_SIZES);
            LOGGER.warn(keysBRCA);
        }

        return ImmutableBRCA.builder().brcaPart1(createBRCAPart1(brcaObject)).brcaPart2(createBRCAPart2(brcaObject)).build();
    }

    @NotNull
    private static BRCApart1 createBRCAPart1(@NotNull JsonObject brcaObject) {
        return ImmutableBRCApart1.builder()
                .variantFrequencyLOVD(string(brcaObject, "Variant_frequency_LOVD"))
                .alleleFrequencyFINExAC(string(brcaObject, "Allele_frequency_FIN_ExAC"))
                .clinVarAccessionENIGMA(string(brcaObject, "ClinVarAccession_ENIGMA"))
                .homozygousCountAFRExAC(string(brcaObject, "Homozygous_count_AFR_ExAC"))
                .bxIdExAC(string(brcaObject, "BX_ID_ExAC"))
                .variantInLOVD(string(brcaObject, "Variant_in_LOVD"))
                .alleleFrequencyAFRExAC(string(brcaObject, "Allele_frequency_AFR_ExAC"))
                .chr(string(brcaObject, "Chr"))
                .bxIdENIGMA(string(brcaObject, "BX_ID_ENIGMA"))
                .cooccurrenceLRExLOVD(string(brcaObject, "Co_occurrence_LR_exLOVD"))
                .homozygousCountEASExAC(string(brcaObject, "Homozygous_count_EAS_ExAC"))
                .submitterClinVar(string(brcaObject, "Submitter_ClinVar"))
                .alleleFrequencyEASExAC(string(brcaObject, "Allele_frequency_EAS_ExAC"))
                .hg37End(string(brcaObject, "Hg37_End"))
                .submittersLOVD(string(brcaObject, "Submitters_LOVD"))
                .clinicalClassificationBIC(string(brcaObject, "Clinical_classification_BIC"))
                .homozygousCountNFEExAC(string(brcaObject, "Homozygous_count_NFE_ExAC"))
                .alleleCountSASExAC(string(brcaObject, "Allele_count_SAS_ExAC"))
                .methodClinVar(string(brcaObject, "Method_ClinVar"))
                .alleleCountNFEExAC(string(brcaObject, "Allele_count_NFE_ExAC"))
                .pathogenicityAll(string(brcaObject, "Pathogenicity_all"))
                .germlineOrSomaticBIC(string(brcaObject, "Germline_or_Somatic_BIC"))
                .homozygousCountSASExAC(string(brcaObject, "Homozygous_count_SAS_ExAC"))
                .bicNomenclature(string(brcaObject, "BIC_Nomenclature"))
                .assertionMethodENIGMA(string(brcaObject, "Assertion_method_ENIGMA"))
                .literatureSourceExLOVD(string(brcaObject, "Literature_source_exLOVD"))
                .changeTypeId(string(brcaObject, "Change_Type_id"))
                .collectionMethodENIGMA(string(brcaObject, "Collection_method_ENIGMA"))
                .sumFamilyLRExLOVD(string(brcaObject, "Sum_family_LR_exLOVD"))
                .hgvsCDNALOVD(string(brcaObject, "HGVS_cDNA_LOVD"))
                .homozygousCountFINExAC(string(brcaObject, "Homozygous_count_FIN_ExAC"))
                .easAlleleFrequency1000Genomes(string(brcaObject, "EAS_Allele_frequency_1000_Genomes"))
                .ethnicityBIC(string(brcaObject, "Ethnicity_BIC"))
                .individualsLOVD(string(brcaObject, "Individuals_LOVD"))
                .variantInExAC(string(brcaObject, "Variant_in_ExAC"))
                .urlENIGMA(string(brcaObject, "URL_ENIGMA"))
                .alleleOriginClinVar(string(brcaObject, "Allele_Origin_ClinVar"))
                .alleleFrequencyAMRExAC(string(brcaObject, "Allele_frequency_AMR_ExAC"))
                .variantIn1000Genomes(string(brcaObject, "Variant_in_1000_Genomes"))
                .afrAlleleFrequency1000Genomes(string(brcaObject, "AFR_Allele_frequency_1000_Genomes"))
                .bxIdExLOVD(string(brcaObject, "BX_ID_exLOVD"))
                .source(string(brcaObject, "Source"))
                .conditionIdValueENIGMA(string(brcaObject, "Condition_ID_value_ENIGMA"))
                .hgvsProtein(string(brcaObject, "HGVS_Protein"))
                .ref(string(brcaObject, "Ref"))
                .alleleNumberAFRExAC(string(brcaObject, "Allele_number_AFR_ExAC"))
                .alleleCountAFRExAC(string(brcaObject, "Allele_count_AFR_ExAC"))
                .bxIdLOVD(string(brcaObject, "BX_ID_LOVD"))
                .synonyms(string(brcaObject, "Synonyms"))
                .geneSymbol(string(brcaObject, "Gene_Symbol"))
                .commentOnClinicalSignificanceENIGMA(string(brcaObject, "Comment_on_clinical_significance_ENIGMA"))
                .missenseAnalysisPriorProbabilityExLOVD(string(brcaObject, "Missense_analysis_prior_probability_exLOVD"))
                .alleleNumberFINExAC(string(brcaObject, "Allele_number_FIN_ExAC"))
                .posteriorProbabilityExLOVD(string(brcaObject, "Posterior_probability_exLOVD"))
                .polyphenScore(string(brcaObject, "Polyphen_Score"))
                .referenceSequence(string(brcaObject, "Reference_Sequence"))
                .alleleCountEASExAC(string(brcaObject, "Allele_count_EAS_ExAC"))
                .hg38End(string(brcaObject, "Hg38_End"))
                .hgvsCDNA(string(brcaObject, "HGVS_cDNA"))
                .functionalAnalysisTechniqueLOVD(string(brcaObject, "Functional_analysis_technique_LOVD"))
                .sasAlleleFrequency1000Genomes(string(brcaObject, "SAS_Allele_frequency_1000_Genomes"))
                .rnaLOVD(string(brcaObject, "RNA_LOVD"))
                .combinedPriorProbabilityExLOVD(string(brcaObject, "Combined_prior_probablility_exLOVD"))
                .bxIdClinVar(string(brcaObject, "BX_ID_ClinVar"))
                .iarcClassExLOVD(string(brcaObject, "IARC_class_exLOVD"))
                .bxIdBIC(string(brcaObject, "BX_ID_BIC"))
                .siftPrediction(string(brcaObject, "Sift_Prediction"))
                .alleleNumberNFEExAC(string(brcaObject, "Allele_number_NFE_ExAC"))
                .alleleOriginENIGMA(string(brcaObject, "Allele_origin_ENIGMA"))
                .alleleNumberOTHExAC(string(brcaObject, "Allele_number_OTH_ExAC"))
                .hg36End(string(brcaObject, "Hg36_End"))
                .alleleFrequencySASExAC(string(brcaObject, "Allele_frequency_SAS_ExAC"))
                .dateLastUpdatedClinVar(string(brcaObject, "Date_Last_Updated_ClinVar"))
                .alleleNumberEASExAC(string(brcaObject, "Allele_number_EAS_ExAC"))
                .alleleFrequencyOTHExAC(string(brcaObject, "Allele_frequency_OTH_ExAC"))
                .sourceURL(string(brcaObject, "Source_URL"))
                .scvClinVar(string(brcaObject, "SCV_ClinVar"))
                .pathogenicityExpert(string(brcaObject, "Pathogenicity_expert"))
                .alleleFrequency1000Genomes(string(brcaObject, "Allele_frequency_1000_Genomes"))
                .functionalAnalysisResultLOVD(string(brcaObject, "Functional_analysis_result_LOVD"))
                .amrAlleleFrequency1000Genomes(string(brcaObject, "AMR_Allele_frequency_1000_Genomes"))
                .variantInESP(string(brcaObject, "Variant_in_ESP"))
                .variantInBIC(string(brcaObject, "Variant_in_BIC"))
                .clinicalSignificanceENIGMA(string(brcaObject, "Clinical_significance_ENIGMA"))
                .maxAlleleFrequency(string(brcaObject, "Max_Allele_Frequency"))
                .alleleCountAMRExAC(string(brcaObject, "Allele_count_AMR_ExAC"))
                .variantInENIGMA(string(brcaObject, "Variant_in_ENIGMA"))
                .bxIdESP(string(brcaObject, "BX_ID_ESP"))
                .patientNationalityBIC(string(brcaObject, "Patient_nationality_BIC"))
                .bxId1000Genomes(string(brcaObject, "BX_ID_1000_Genomes"))
                .genomicCoordinateHg37(string(brcaObject, "Genomic_Coordinate_hg37"))
                .genomicCoordinateHg36(string(brcaObject, "Genomic_Coordinate_hg36"))
                .eurAlleleFrequency1000Genomes(string(brcaObject, "EUR_Allele_frequency_1000_Genomes"))
                .numberOfFamilyMemberCarryingMutationBIC(string(brcaObject, "Number_of_family_member_carrying_mutation_BIC"))
                .segregationLRExLOVD(string(brcaObject, "Segregation_LR_exLOVD"))
                .alleleFrequency(string(brcaObject, "Allele_Frequency"))
                .minorAlleleFrequencyPercentESP(string(brcaObject, "Minor_allele_frequency_percent_ESP"))
                .alleleFrequencyExAC(string(brcaObject, "Allele_frequency_ExAC"))
                .mutationTypeBIC(string(brcaObject, "Mutation_type_BIC"))
                .assertionMethodCitationENIGMA(string(brcaObject, "Assertion_method_citation_ENIGMA"))
                .conditionIdTypeENIGMA(string(brcaObject, "Condition_ID_type_ENIGMA"))
                .alleleCountOTHExAC(string(brcaObject, "Allele_count_OTH_ExAC"))
                .hgvsProteinLOVD(string(brcaObject, "HGVS_protein_LOVD"))
                .variantInClinVar(string(brcaObject, "Variant_in_ClinVar"))
                .clinicalImportanceBIC(string(brcaObject, "Clinical_importance_BIC"))
                .discordant(string(brcaObject, "Discordant"))
                .build();
    }

    @NotNull
    private static BRCApart2 createBRCAPart2(@NotNull JsonObject brcaObject) {
        return ImmutableBRCApart2.builder()
                .alleleCountFINExAC(string(brcaObject, "Allele_count_FIN_ExAC"))
                .conditionCategoryENIGMA(string(brcaObject, "Condition_category_ENIGMA"))
                .alleleFrequencyESP(string(brcaObject, "Allele_Frequency_ESP"))
                .homozygousCountOTHExAC(string(brcaObject, "Homozygous_count_OTH_ExAC"))
                .geneticOriginLOVD(string(brcaObject, "Genetic_origin_LOVD"))
                .id(string(brcaObject, "id"))
                .homozygousCountAMRExAC(string(brcaObject, "Homozygous_count_AMR_ExAC"))
                .clinicalSignificanceClinVar(string(brcaObject, "Clinical_Significance_ClinVar"))
                .aaAlleleFrequencyESP(string(brcaObject, "AA_Allele_Frequency_ESP"))
                .proteinChange(string(brcaObject, "Protein_Change"))
                .variantInExLOVD(string(brcaObject, "Variant_in_exLOVD"))
                .eaAlleleFrequencyESP(string(brcaObject, "EA_Allele_Frequency_ESP"))
                .hgvsRNA(string(brcaObject, "HGVS_RNA"))
                .clinicalSignificanceCitationsENIGMA(string(brcaObject, "Clinical_significance_citations_ENIGMA"))
                .variantEffectLOVD(string(brcaObject, "Variant_effect_LOVD"))
                .polyphenPrediction(string(brcaObject, "Polyphen_Prediction"))
                .dataReleaseId(string(brcaObject, "Data_Release_id"))
                .hg37Start(string(brcaObject, "Hg37_Start"))
                .hg36Start(string(brcaObject, "Hg36_Start"))
                .siftScore(string(brcaObject, "Sift_Score"))
                .genomicCoordinateHg38(string(brcaObject, "Genomic_Coordinate_hg38"))
                .alt(string(brcaObject, "Alt"))
                .literatureCitationBIC(string(brcaObject, "Literature_citation_BIC"))
                .variantHaplotypeLOVD(string(brcaObject, "Variant_haplotype_LOVD"))
                .alleleFrequencyNFEExAC(string(brcaObject, "Allele_frequency_NFE_ExAC"))
                .hg38Start(string(brcaObject, "Hg38_Start"))
                .pos(string(brcaObject, "Pos"))
                .dateLastEvaluatedENIGMA(string(brcaObject, "Date_last_evaluated_ENIGMA"))
                .alleleNumberSASExAC(string(brcaObject, "Allele_number_SAS_ExAC"))
                .alleleNumberAMRExAC(string(brcaObject, "Allele_number_AMR_ExAC"))
                .dbIdLOVD(string(brcaObject, "DBID_LOVD"))
                .build();
    }
}
