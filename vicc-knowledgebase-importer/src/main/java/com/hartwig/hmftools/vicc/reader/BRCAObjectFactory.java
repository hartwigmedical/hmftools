package com.hartwig.hmftools.vicc.reader;

import static com.hartwig.hmftools.common.utils.json.JsonFunctions.string;

import com.google.gson.JsonObject;
import com.hartwig.hmftools.vicc.datamodel.brca.Brca;
import com.hartwig.hmftools.vicc.datamodel.brca.BrcaAnnotation1000Genomes;
import com.hartwig.hmftools.vicc.datamodel.brca.BrcaAnnotationBIC;
import com.hartwig.hmftools.vicc.datamodel.brca.BrcaAnnotationClinVar;
import com.hartwig.hmftools.vicc.datamodel.brca.BrcaAnnotationENIGMA;
import com.hartwig.hmftools.vicc.datamodel.brca.BrcaAnnotationESP;
import com.hartwig.hmftools.vicc.datamodel.brca.BrcaAnnotationExAC;
import com.hartwig.hmftools.vicc.datamodel.brca.BrcaAnnotationExLOVD;
import com.hartwig.hmftools.vicc.datamodel.brca.BrcaAnnotationLOVD;
import com.hartwig.hmftools.vicc.datamodel.brca.ImmutableBrca;
import com.hartwig.hmftools.vicc.datamodel.brca.ImmutableBrcaAnnotation1000Genomes;
import com.hartwig.hmftools.vicc.datamodel.brca.ImmutableBrcaAnnotationBIC;
import com.hartwig.hmftools.vicc.datamodel.brca.ImmutableBrcaAnnotationClinVar;
import com.hartwig.hmftools.vicc.datamodel.brca.ImmutableBrcaAnnotationENIGMA;
import com.hartwig.hmftools.vicc.datamodel.brca.ImmutableBrcaAnnotationESP;
import com.hartwig.hmftools.vicc.datamodel.brca.ImmutableBrcaAnnotationExAC;
import com.hartwig.hmftools.vicc.datamodel.brca.ImmutableBrcaAnnotationExLOVD;
import com.hartwig.hmftools.vicc.datamodel.brca.ImmutableBrcaAnnotationLOVD;

import org.jetbrains.annotations.NotNull;

final class BRCAObjectFactory {

    private BRCAObjectFactory() {
    }

    @NotNull
    static Brca create(@NotNull JsonObject brcaObject) {
        ViccDatamodelCheckerFactory.brcaEntryChecker().check(brcaObject);

        return ImmutableBrca.builder()
                .geneSymbol(string(brcaObject, "Gene_Symbol"))
                .chr(string(brcaObject, "Chr"))
                .pos(string(brcaObject, "Pos"))
                .ref(string(brcaObject, "Ref"))
                .alt(string(brcaObject, "Alt"))
                .genomicCoordinateHg36(string(brcaObject, "Genomic_Coordinate_hg36"))
                .hg36Start(string(brcaObject, "Hg36_Start"))
                .hg36End(string(brcaObject, "Hg36_End"))
                .genomicCoordinateHg37(string(brcaObject, "Genomic_Coordinate_hg37"))
                .hg37Start(string(brcaObject, "Hg37_Start"))
                .hg37End(string(brcaObject, "Hg37_End"))
                .genomicCoordinateHg38(string(brcaObject, "Genomic_Coordinate_hg38"))
                .hg38Start(string(brcaObject, "Hg38_Start"))
                .hg38End(string(brcaObject, "Hg38_End"))
                .proteinChange(string(brcaObject, "Protein_Change"))
                .referenceSequence(string(brcaObject, "Reference_Sequence"))
                .synonyms(string(brcaObject, "Synonyms"))
                .hgvsCDNA(string(brcaObject, "HGVS_cDNA"))
                .hgvsProtein(string(brcaObject, "HGVS_Protein"))
                .hgvsRNA(string(brcaObject, "HGVS_RNA"))
                .siftScore(string(brcaObject, "Sift_Score"))
                .siftPrediction(string(brcaObject, "Sift_Prediction"))
                .polyphenScore(string(brcaObject, "Polyphen_Score"))
                .polyphenPrediction(string(brcaObject, "Polyphen_Prediction"))
                .pathogenicityAll(string(brcaObject, "Pathogenicity_all"))
                .pathogenicityExpert(string(brcaObject, "Pathogenicity_expert"))
                .alleleFrequency(string(brcaObject, "Allele_Frequency"))
                .maxAlleleFrequency(string(brcaObject, "Max_Allele_Frequency"))
                .discordant(string(brcaObject, "Discordant"))
                .id(string(brcaObject, "id"))
                .changeTypeId(string(brcaObject, "Change_Type_id"))
                .dataReleaseId(string(brcaObject, "Data_Release_id"))
                .source(string(brcaObject, "Source"))
                .sourceURL(string(brcaObject, "Source_URL"))
                .annotation1000Genomes(createBRCAAnnotation1000Genomes(brcaObject))
                .annotationBIC(createBRCAAnnotationBIC(brcaObject))
                .annotationClinVar(createBRCAAnnotationClinVar(brcaObject))
                .annotationENIGMA(createBRCAAnnotationENIGMA(brcaObject))
                .annotationESP(createBRCAAnnotationESP(brcaObject))
                .annotationExAC(createBRCAAnnotationExAC(brcaObject))
                .annotationExLOVD(createBRCAAnnotationExLOVD(brcaObject))
                .annotationLOVD(createBRCAAnnotationLOVD(brcaObject))
                .build();
    }

    @NotNull
    private static BrcaAnnotation1000Genomes createBRCAAnnotation1000Genomes(@NotNull JsonObject brcaObject) {
        return ImmutableBrcaAnnotation1000Genomes.builder()
                .variantIn1000Genomes(string(brcaObject, "Variant_in_1000_Genomes"))
                .bxId(string(brcaObject, "BX_ID_1000_Genomes"))
                .alleleFrequency(string(brcaObject, "Allele_frequency_1000_Genomes"))
                .afrAlleleFrequency(string(brcaObject, "AFR_Allele_frequency_1000_Genomes"))
                .amrAlleleFrequency(string(brcaObject, "AMR_Allele_frequency_1000_Genomes"))
                .easAlleleFrequency(string(brcaObject, "EAS_Allele_frequency_1000_Genomes"))
                .eurAlleleFrequency(string(brcaObject, "EUR_Allele_frequency_1000_Genomes"))
                .sasAlleleFrequency(string(brcaObject, "SAS_Allele_frequency_1000_Genomes"))
                .build();
    }

    @NotNull
    private static BrcaAnnotationBIC createBRCAAnnotationBIC(@NotNull JsonObject brcaObject) {
        return ImmutableBrcaAnnotationBIC.builder()
                .variantInBIC(string(brcaObject, "Variant_in_BIC"))
                .bxId(string(brcaObject, "BX_ID_BIC"))
                .mutationType(string(brcaObject, "Mutation_type_BIC"))
                .clinicalClassification(string(brcaObject, "Clinical_classification_BIC"))
                .clinicalImportance(string(brcaObject, "Clinical_importance_BIC"))
                .nomenclature(string(brcaObject, "BIC_Nomenclature"))
                .ethnicity(string(brcaObject, "Ethnicity_BIC"))
                .patientNationality(string(brcaObject, "Patient_nationality_BIC"))
                .germlineOrSomatic(string(brcaObject, "Germline_or_Somatic_BIC"))
                .numberOfFamilyMemberCarryingMutation(string(brcaObject, "Number_of_family_member_carrying_mutation_BIC"))
                .literatureCitation(string(brcaObject, "Literature_citation_BIC"))
                .build();
    }

    @NotNull
    private static BrcaAnnotationClinVar createBRCAAnnotationClinVar(@NotNull JsonObject brcaObject) {
        return ImmutableBrcaAnnotationClinVar.builder()
                .variantInClinVar(string(brcaObject, "Variant_in_ClinVar"))
                .bxId(string(brcaObject, "BX_ID_ClinVar"))
                .clinicalSignificance(string(brcaObject, "Clinical_Significance_ClinVar"))
                .submitter(string(brcaObject, "Submitter_ClinVar"))
                .method(string(brcaObject, "Method_ClinVar"))
                .alleleOrigin(string(brcaObject, "Allele_Origin_ClinVar"))
                .scv(string(brcaObject, "SCV_ClinVar"))
                .dateLastUpdated(string(brcaObject, "Date_Last_Updated_ClinVar"))
                .build();
    }

    @NotNull
    private static BrcaAnnotationENIGMA createBRCAAnnotationENIGMA(@NotNull JsonObject brcaObject) {
        return ImmutableBrcaAnnotationENIGMA.builder()
                .variantInENIGMA(string(brcaObject, "Variant_in_ENIGMA"))
                .bxId(string(brcaObject, "BX_ID_ENIGMA"))
                .alleleOrigin(string(brcaObject, "Allele_origin_ENIGMA"))
                .clinVarAccession(string(brcaObject, "ClinVarAccession_ENIGMA"))
                .assertionMethod(string(brcaObject, "Assertion_method_ENIGMA"))
                .assertionMethodCitation(string(brcaObject, "Assertion_method_citation_ENIGMA"))
                .collectionMethod(string(brcaObject, "Collection_method_ENIGMA"))
                .conditionCategory(string(brcaObject, "Condition_category_ENIGMA"))
                .conditionIdValue(string(brcaObject, "Condition_ID_value_ENIGMA"))
                .conditionIdType(string(brcaObject, "Condition_ID_type_ENIGMA"))
                .clinicalSignificance(string(brcaObject, "Clinical_significance_ENIGMA"))
                .clinicalSignificanceCitations(string(brcaObject, "Clinical_significance_citations_ENIGMA"))
                .commentOnClinicalSignificance(string(brcaObject, "Comment_on_clinical_significance_ENIGMA"))
                .dateLastEvaluated(string(brcaObject, "Date_last_evaluated_ENIGMA"))
                .url(string(brcaObject, "URL_ENIGMA"))
                .build();
    }

    @NotNull
    private static BrcaAnnotationESP createBRCAAnnotationESP(@NotNull JsonObject brcaObject) {
        return ImmutableBrcaAnnotationESP.builder()
                .variantInESP(string(brcaObject, "Variant_in_ESP"))
                .bxId(string(brcaObject, "BX_ID_ESP"))
                .minorAlleleFrequencyPercent(string(brcaObject, "Minor_allele_frequency_percent_ESP"))
                .alleleFrequency(string(brcaObject, "Allele_Frequency_ESP"))
                .aaAlleleFrequency(string(brcaObject, "AA_Allele_Frequency_ESP"))
                .eaAlleleFrequency(string(brcaObject, "EA_Allele_Frequency_ESP"))
                .build();
    }

    @NotNull
    private static BrcaAnnotationExAC createBRCAAnnotationExAC(@NotNull JsonObject brcaObject) {
        return ImmutableBrcaAnnotationExAC.builder()
                .variantInExAC(string(brcaObject, "Variant_in_ExAC"))
                .bxId(string(brcaObject, "BX_ID_ExAC"))
                .alleleFrequency(string(brcaObject, "Allele_frequency_ExAC"))
                .alleleFrequencyAFR(string(brcaObject, "Allele_frequency_AFR_ExAC"))
                .alleleFrequencyAMR(string(brcaObject, "Allele_frequency_AMR_ExAC"))
                .alleleFrequencyEAS(string(brcaObject, "Allele_frequency_EAS_ExAC"))
                .alleleFrequencyFIN(string(brcaObject, "Allele_frequency_FIN_ExAC"))
                .alleleFrequencyNFE(string(brcaObject, "Allele_frequency_NFE_ExAC"))
                .alleleFrequencyOTH(string(brcaObject, "Allele_frequency_OTH_ExAC"))
                .alleleFrequencySAS(string(brcaObject, "Allele_frequency_SAS_ExAC"))
                .alleleNumberAFR(string(brcaObject, "Allele_number_AFR_ExAC"))
                .alleleNumberAMR(string(brcaObject, "Allele_number_AMR_ExAC"))
                .alleleNumberEAS(string(brcaObject, "Allele_number_EAS_ExAC"))
                .alleleNumberFIN(string(brcaObject, "Allele_number_FIN_ExAC"))
                .alleleNumberNFE(string(brcaObject, "Allele_number_NFE_ExAC"))
                .alleleNumberOTH(string(brcaObject, "Allele_number_OTH_ExAC"))
                .alleleNumberSAS(string(brcaObject, "Allele_number_SAS_ExAC"))
                .homozygousCountAFR(string(brcaObject, "Homozygous_count_AFR_ExAC"))
                .homozygousCountAMR(string(brcaObject, "Homozygous_count_AMR_ExAC"))
                .homozygousCountEAS(string(brcaObject, "Homozygous_count_EAS_ExAC"))
                .homozygousCountFIN(string(brcaObject, "Homozygous_count_FIN_ExAC"))
                .homozygousCountNFE(string(brcaObject, "Homozygous_count_NFE_ExAC"))
                .homozygousCountOTH(string(brcaObject, "Homozygous_count_OTH_ExAC"))
                .homozygousCountSAS(string(brcaObject, "Homozygous_count_SAS_ExAC"))
                .alleleCountAFR(string(brcaObject, "Allele_count_AFR_ExAC"))
                .alleleCountAMR(string(brcaObject, "Allele_count_AMR_ExAC"))
                .alleleCountEAS(string(brcaObject, "Allele_count_EAS_ExAC"))
                .alleleCountFIN(string(brcaObject, "Allele_count_FIN_ExAC"))
                .alleleCountNFE(string(brcaObject, "Allele_count_NFE_ExAC"))
                .alleleCountOTH(string(brcaObject, "Allele_count_OTH_ExAC"))
                .alleleCountSAS(string(brcaObject, "Allele_count_SAS_ExAC"))
                .build();
    }

    @NotNull
    private static BrcaAnnotationExLOVD createBRCAAnnotationExLOVD(@NotNull JsonObject brcaObject) {
        return ImmutableBrcaAnnotationExLOVD.builder()
                .variantInExLOVD(string(brcaObject, "Variant_in_exLOVD"))
                .bxId(string(brcaObject, "BX_ID_exLOVD"))
                .cooccurrenceLR(string(brcaObject, "Co_occurrence_LR_exLOVD"))
                .sumFamilyLR(string(brcaObject, "Sum_family_LR_exLOVD"))
                .segregationLR(string(brcaObject, "Segregation_LR_exLOVD"))
                .posteriorProbability(string(brcaObject, "Posterior_probability_exLOVD"))
                .missenseAnalysisPriorProbability(string(brcaObject, "Missense_analysis_prior_probability_exLOVD"))
                .combinedPriorProbability(string(brcaObject, "Combined_prior_probablility_exLOVD"))
                .iarcClass(string(brcaObject, "IARC_class_exLOVD"))
                .literatureSource(string(brcaObject, "Literature_source_exLOVD"))
                .build();
    }

    @NotNull
    private static BrcaAnnotationLOVD createBRCAAnnotationLOVD(@NotNull JsonObject brcaObject) {
        return ImmutableBrcaAnnotationLOVD.builder()
                .variantInLOVD(string(brcaObject, "Variant_in_LOVD"))
                .bxId(string(brcaObject, "BX_ID_LOVD"))
                .dbId(string(brcaObject, "DBID_LOVD"))
                .hgvsCDNA(string(brcaObject, "HGVS_cDNA_LOVD"))
                .hgvsProtein(string(brcaObject, "HGVS_protein_LOVD"))
                .rna(string(brcaObject, "RNA_LOVD"))
                .variantEffect(string(brcaObject, "Variant_effect_LOVD"))
                .variantFrequency(string(brcaObject, "Variant_frequency_LOVD"))
                .variantHaplotype(string(brcaObject, "Variant_haplotype_LOVD"))
                .geneticOrigin(string(brcaObject, "Genetic_origin_LOVD"))
                .functionalAnalysisTechnique(string(brcaObject, "Functional_analysis_technique_LOVD"))
                .functionalAnalysisResult(string(brcaObject, "Functional_analysis_result_LOVD"))
                .submitters(string(brcaObject, "Submitters_LOVD"))
                .individuals(string(brcaObject, "Individuals_LOVD"))
                .build();
    }
}
