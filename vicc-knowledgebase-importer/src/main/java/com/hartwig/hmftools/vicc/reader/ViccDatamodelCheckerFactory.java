package com.hartwig.hmftools.vicc.reader;

import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.utils.json.JsonDatamodelChecker;

import org.jetbrains.annotations.NotNull;

final class ViccDatamodelCheckerFactory {

    private ViccDatamodelCheckerFactory() {
    }

    @NotNull
    static JsonDatamodelChecker viccEntryChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("source", true);
        map.put("genes", true);
        map.put("gene_identifiers", true);
        map.put("feature_names", false);
        map.put("features", true);
        map.put("association", true);
        map.put("tags", true);
        map.put("dev_tags", true);
        map.put("cgi", false);
        map.put("brca", false);
        map.put("sage", false);
        map.put("pmkb", false);
        map.put("oncokb", false);
        map.put("jax", false);
        map.put("jax_trials", false);
        map.put("molecularmatch", false);
        map.put("molecularmatch_trials", false);
        map.put("civic", false);
        return new JsonDatamodelChecker("ViccEntry", map);
    }

    @NotNull
    static JsonDatamodelChecker geneIdentifierChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("symbol", true);
        map.put("entrez_id", true);
        map.put("ensembl_gene_id", true);
        return new JsonDatamodelChecker("GeneIdentifier", map);
    }

    @NotNull
    static JsonDatamodelChecker featureChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("name", true);
        map.put("biomarker_type", false);
        map.put("referenceName", false);
        map.put("chromosome", false);
        map.put("start", false);
        map.put("end", false);
        map.put("ref", false);
        map.put("alt", false);
        map.put("provenance", false);
        map.put("provenance_rule", false);
        map.put("geneSymbol", false);
        map.put("synonyms", false);
        map.put("entrez_id", false);
        map.put("sequence_ontology", false);
        map.put("links", false);
        map.put("description", false);
        map.put("attributes", false);
        map.put("info", false);
        return new JsonDatamodelChecker("Feature", map);
    }

    @NotNull
    static JsonDatamodelChecker featureInfoChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("germline_or_somatic", true);
        return new JsonDatamodelChecker("FeatureInfo", map);
    }

    @NotNull
    static JsonDatamodelChecker featureAttributeChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("amino_acid_change", true);
        map.put("germline", true);
        map.put("partner_gene", true);
        map.put("description", true);
        map.put("exons", true);
        map.put("notes", true);
        map.put("cosmic", true);
        map.put("effect", true);
        map.put("cnv_type", true);
        map.put("id", true);
        map.put("cytoband", true);
        map.put("variant_type", true);
        map.put("dna_change", true);
        map.put("codons", true);
        map.put("chromosome_based_cnv", true);
        map.put("transcript", true);
        map.put("description_type", true);
        map.put("chromosome", true);
        return new JsonDatamodelChecker("FeatureAttribute", map);
    }

    @NotNull
    static JsonDatamodelChecker featureAttributeObjectChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("string_value", true);
        return new JsonDatamodelChecker("FeatureAttributeObject", map);
    }

    @NotNull
    static JsonDatamodelChecker sequenceOntologyChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("hierarchy", false);
        map.put("soid", true);
        map.put("parent_soid", true);
        map.put("name", true);
        map.put("parent_name", true);
        return new JsonDatamodelChecker("SequenceOntology", map);
    }

    @NotNull
    static JsonDatamodelChecker associationChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("variant_name", false);
        map.put("evidence", true);
        map.put("evidence_level", false);
        map.put("evidence_label", false);
        map.put("response_type", false);
        map.put("drug_labels", false);
        map.put("source_link", false);
        map.put("publication_url", false);
        map.put("phenotype", false);
        map.put("description", true);
        map.put("environmentalContexts", false);
        map.put("oncogenic", false);
        return new JsonDatamodelChecker("Association", map);
    }

    @NotNull
    static JsonDatamodelChecker evidenceChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("info", true);
        map.put("evidenceType", true);
        map.put("description", true);
        return new JsonDatamodelChecker("Evidence", map);
    }

    @NotNull
    static JsonDatamodelChecker evidenceInfoChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("publications", true);
        return new JsonDatamodelChecker("EvidenceInfo", map);
    }

    @NotNull
    static JsonDatamodelChecker evidenceTypeChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("sourceName", true);
        map.put("id", false);
        return new JsonDatamodelChecker("EvidenceType", map);
    }

    @NotNull
    static JsonDatamodelChecker environmentalContextChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("term", false);
        map.put("description", true);
        map.put("taxonomy", false);
        map.put("source", false);
        map.put("usan_stem", false);
        map.put("approved_countries", false);
        map.put("toxicity", false);
        map.put("id", false);
        return new JsonDatamodelChecker("EnvironmentalContext", map);
    }

    @NotNull
    static JsonDatamodelChecker taxonomyChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("kingdom", true);
        map.put("direct-parent", true);
        map.put("class", true);
        map.put("subclass", false);
        map.put("superclass", true);
        return new JsonDatamodelChecker("Taxonomy", map);
    }

    @NotNull
    static JsonDatamodelChecker phenotypeChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("type", false);
        map.put("description", true);
        map.put("family", true);
        map.put("id", false);
        return new JsonDatamodelChecker("Phenotype", map);
    }

    @NotNull
    static JsonDatamodelChecker phenotypeTypeChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("source", true);
        map.put("term", true);
        map.put("id", true);
        return new JsonDatamodelChecker("PhenotypeType", map);
    }

    @NotNull
    static JsonDatamodelChecker brcaEntryChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("Gene_Symbol", true);
        map.put("Chr", true);
        map.put("Pos", true);
        map.put("Ref", true);
        map.put("Alt", true);
        map.put("Genomic_Coordinate_hg36", true);
        map.put("Hg36_Start", true);
        map.put("Hg36_End", true);
        map.put("Genomic_Coordinate_hg37", true);
        map.put("Hg37_Start", true);
        map.put("Hg37_End", true);
        map.put("Genomic_Coordinate_hg38", true);
        map.put("Hg38_Start", true);
        map.put("Hg38_End", true);
        map.put("Protein_Change", true);
        map.put("Reference_Sequence", true);
        map.put("Synonyms", true);
        map.put("HGVS_cDNA", true);
        map.put("HGVS_Protein", true);
        map.put("HGVS_RNA", true);
        map.put("Sift_Score", true);
        map.put("Sift_Prediction", true);
        map.put("Polyphen_Score", true);
        map.put("Polyphen_Prediction", true);
        map.put("Pathogenicity_all", true);
        map.put("Pathogenicity_expert", true);
        map.put("Allele_Frequency", true);
        map.put("Max_Allele_Frequency", true);
        map.put("Discordant", true);
        map.put("id", true);
        map.put("Change_Type_id", true);
        map.put("Data_Release_id", true);
        map.put("Source", true);
        map.put("Source_URL", true);

        map.put("Variant_in_1000_Genomes", true);
        map.put("BX_ID_1000_Genomes", true);
        map.put("Allele_frequency_1000_Genomes", true);
        map.put("AFR_Allele_frequency_1000_Genomes", true);
        map.put("AMR_Allele_frequency_1000_Genomes", true);
        map.put("EAS_Allele_frequency_1000_Genomes", true);
        map.put("EUR_Allele_frequency_1000_Genomes", true);
        map.put("SAS_Allele_frequency_1000_Genomes", true);

        map.put("Variant_in_BIC", true);
        map.put("BX_ID_BIC", true);
        map.put("Mutation_type_BIC", true);
        map.put("Clinical_classification_BIC", true);
        map.put("Clinical_importance_BIC", true);
        map.put("BIC_Nomenclature", true);
        map.put("Ethnicity_BIC", true);
        map.put("Patient_nationality_BIC", true);
        map.put("Germline_or_Somatic_BIC", true);
        map.put("Number_of_family_member_carrying_mutation_BIC", true);
        map.put("Literature_citation_BIC", true);

        map.put("Variant_in_ClinVar", true);
        map.put("BX_ID_ClinVar", true);
        map.put("Clinical_Significance_ClinVar", true);
        map.put("Submitter_ClinVar", true);
        map.put("Method_ClinVar", true);
        map.put("Allele_Origin_ClinVar", true);
        map.put("SCV_ClinVar", true);
        map.put("Date_Last_Updated_ClinVar", true);

        map.put("Variant_in_ENIGMA", true);
        map.put("BX_ID_ENIGMA", true);
        map.put("Allele_origin_ENIGMA", true);
        map.put("ClinVarAccession_ENIGMA", true);
        map.put("Assertion_method_ENIGMA", true);
        map.put("Assertion_method_citation_ENIGMA", true);
        map.put("Collection_method_ENIGMA", true);
        map.put("Condition_category_ENIGMA", true);
        map.put("Condition_ID_value_ENIGMA", true);
        map.put("Condition_ID_type_ENIGMA", true);
        map.put("Clinical_significance_ENIGMA", true);
        map.put("Clinical_significance_citations_ENIGMA", true);
        map.put("Comment_on_clinical_significance_ENIGMA", true);
        map.put("Date_last_evaluated_ENIGMA", true);
        map.put("URL_ENIGMA", true);

        map.put("Variant_in_ESP", true);
        map.put("BX_ID_ESP", true);
        map.put("Minor_allele_frequency_percent_ESP", true);
        map.put("Allele_Frequency_ESP", true);
        map.put("AA_Allele_Frequency_ESP", true);
        map.put("EA_Allele_Frequency_ESP", true);

        map.put("Variant_in_ExAC", true);
        map.put("BX_ID_ExAC", true);
        map.put("Allele_frequency_ExAC", true);
        map.put("Allele_frequency_AFR_ExAC", true);
        map.put("Allele_frequency_AMR_ExAC", true);
        map.put("Allele_frequency_EAS_ExAC", true);
        map.put("Allele_frequency_FIN_ExAC", true);
        map.put("Allele_frequency_NFE_ExAC", true);
        map.put("Allele_frequency_OTH_ExAC", true);
        map.put("Allele_frequency_SAS_ExAC", true);
        map.put("Allele_number_AFR_ExAC", true);
        map.put("Allele_number_AMR_ExAC", true);
        map.put("Allele_number_EAS_ExAC", true);
        map.put("Allele_number_FIN_ExAC", true);
        map.put("Allele_number_NFE_ExAC", true);
        map.put("Allele_number_OTH_ExAC", true);
        map.put("Allele_number_SAS_ExAC", true);
        map.put("Homozygous_count_AFR_ExAC", true);
        map.put("Homozygous_count_AMR_ExAC", true);
        map.put("Homozygous_count_EAS_ExAC", true);
        map.put("Homozygous_count_FIN_ExAC", true);
        map.put("Homozygous_count_NFE_ExAC", true);
        map.put("Homozygous_count_OTH_ExAC", true);
        map.put("Homozygous_count_SAS_ExAC", true);
        map.put("Allele_count_AFR_ExAC", true);
        map.put("Allele_count_AMR_ExAC", true);
        map.put("Allele_count_EAS_ExAC", true);
        map.put("Allele_count_FIN_ExAC", true);
        map.put("Allele_count_NFE_ExAC", true);
        map.put("Allele_count_OTH_ExAC", true);
        map.put("Allele_count_SAS_ExAC", true);

        map.put("Variant_in_exLOVD", true);
        map.put("BX_ID_exLOVD", true);
        map.put("Co_occurrence_LR_exLOVD", true);
        map.put("Sum_family_LR_exLOVD", true);
        map.put("Segregation_LR_exLOVD", true);
        map.put("Posterior_probability_exLOVD", true);
        map.put("Missense_analysis_prior_probability_exLOVD", true);
        map.put("Combined_prior_probablility_exLOVD", true);
        map.put("IARC_class_exLOVD", true);
        map.put("Literature_source_exLOVD", true);

        map.put("Variant_in_LOVD", true);
        map.put("BX_ID_LOVD", true);
        map.put("DBID_LOVD", true);
        map.put("HGVS_cDNA_LOVD", true);
        map.put("HGVS_protein_LOVD", true);
        map.put("RNA_LOVD", true);
        map.put("Variant_effect_LOVD", true);
        map.put("Variant_frequency_LOVD", true);
        map.put("Variant_haplotype_LOVD", true);
        map.put("Genetic_origin_LOVD", true);
        map.put("Functional_analysis_technique_LOVD", true);
        map.put("Functional_analysis_result_LOVD", true);
        map.put("Submitters_LOVD", true);
        map.put("Individuals_LOVD", true);

        return new JsonDatamodelChecker("BrcaEntry", map);
    }

    @NotNull
    static JsonDatamodelChecker civicEntryChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("entrez_id", true);
        map.put("entrez_name", true);
        map.put("name", true);
        map.put("type", true);
        map.put("coordinates", true);
        map.put("sources", true);
        map.put("variant_groups", false);
        map.put("variant_types", true);
        map.put("civic_actionability_score", true);
        map.put("clinvar_entries", true);
        map.put("lifecycle_actions", true);
        map.put("variant_aliases", true);
        map.put("allele_registry_id", true);
        map.put("provisional_values", true);
        map.put("gene_id", true);
        map.put("evidence_items", true);
        map.put("assertions", true);
        map.put("hgvs_expressions", true);
        map.put("errors", true);
        map.put("id", true);
        map.put("description", true);
        return new JsonDatamodelChecker("CivicEntry", map);
    }

    @NotNull
    static JsonDatamodelChecker civicCoordinatesChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("chromosome", true);
        map.put("start", true);
        map.put("stop", true);
        map.put("reference_bases", true);
        map.put("variant_bases", true);
        map.put("representative_transcript", true);
        map.put("ensembl_version", true);
        map.put("reference_build", false);
        map.put("chromosome2", true);
        map.put("start2", true);
        map.put("stop2", true);
        map.put("representative_transcript2", true);
        return new JsonDatamodelChecker("CivicCoordinates", map);
    }

    @NotNull
    static JsonDatamodelChecker civicSourceChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("status", true);
        map.put("open_access", true);
        map.put("name", true);
        map.put("journal", true);
        map.put("citation", true);
        map.put("pmc_id", true);
        map.put("full_journal_title", true);
        map.put("source_url", true);
        map.put("clinical_trials", true);
        map.put("pubmed_id", true);
        map.put("is_review", true);
        map.put("publication_date", true);
        map.put("id", true);
        return new JsonDatamodelChecker("CivicSource", map);
    }

    @NotNull
    static JsonDatamodelChecker civicProvisionalValueChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("description", false);
        return new JsonDatamodelChecker("CivicProvisionalValue", map);
    }

    @NotNull
    static JsonDatamodelChecker civicProvisionalValueDescriptionChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("revision_id", false);
        map.put("value", false);
        return new JsonDatamodelChecker("CivicProvisionalValueDescription", map);
    }

    @NotNull
    static JsonDatamodelChecker civicVariantGroupChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("name", true);
        map.put("type", true);
        map.put("description", true);
        map.put("variants", true);
        map.put("id", true);
        return new JsonDatamodelChecker("CivicVariantGroup", map);
    }

    @NotNull
    static JsonDatamodelChecker civicVariantChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("entrez_id", true);
        map.put("entrez_name", true);
        map.put("name", true);
        map.put("type", true);
        map.put("variant_types", true);
        map.put("civic_actionability_score", false);
        map.put("coordinates", false);
        map.put("id", true);
        map.put("gene_id", true);
        map.put("description", true);
        return new JsonDatamodelChecker("CivicVariant", map);
    }

    @NotNull
    static JsonDatamodelChecker civicVariantTypeChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("name", true);
        map.put("display_name", true);
        map.put("description", true);
        map.put("url", true);
        map.put("so_id", true);
        map.put("id", true);
        return new JsonDatamodelChecker("CivicVariantType", map);
    }

    @NotNull
    static JsonDatamodelChecker civicEvidenceItemChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("name", true);
        map.put("type", true);
        map.put("status", true);
        map.put("rating", true);
        map.put("evidence_type", true);
        map.put("evidence_level", true);
        map.put("evidence_direction", true);
        map.put("drug_interaction_type", true);
        map.put("drugs", true);
        map.put("disease", true);
        map.put("variant_origin", true);
        map.put("source", true);
        map.put("clinical_significance", true);
        map.put("open_change_count", true);
        map.put("description", true);
        map.put("variant_id", false);
        map.put("id", true);
        return new JsonDatamodelChecker("CivicEvidenceItem", map);
    }

    @NotNull
    static JsonDatamodelChecker civicClinicalTrialChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("name", true);
        map.put("nct_id", true);
        map.put("clinical_trial_url", true);
        map.put("description", true);
        return new JsonDatamodelChecker("CivicClinicalTrial", map);
    }

    @NotNull
    static JsonDatamodelChecker civicPublicationDateChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("year", false);
        map.put("month", false);
        map.put("day", false);
        return new JsonDatamodelChecker("CivicPublicationDate", map);
    }

    @NotNull
    static JsonDatamodelChecker civicDiseaseChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("name", true);
        map.put("display_name", true);
        map.put("doid", true);
        map.put("url", true);
        map.put("id", true);
        return new JsonDatamodelChecker("CivicDisease", map);
    }

    @NotNull
    static JsonDatamodelChecker civicDrugChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("name", true);
        map.put("pubchem_id", true);
        map.put("id", true);
        return new JsonDatamodelChecker("CivicDrug", map);
    }

    @NotNull
    static JsonDatamodelChecker civicLifecycleActionsChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("last_commented_on", false);
        map.put("last_modified", false);
        map.put("last_reviewed", false);
        return new JsonDatamodelChecker("CivicLifecycleActions", map);
    }

    @NotNull
    static JsonDatamodelChecker civicLastCommentedOnChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("timestamp", true);
        map.put("user", true);
        return new JsonDatamodelChecker("CivicLastCommentedOn", map);
    }

    @NotNull
    static JsonDatamodelChecker civicLastModifiedChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("timestamp", true);
        map.put("user", true);
        return new JsonDatamodelChecker("CivicLastModified", map);
    }

    @NotNull
    static JsonDatamodelChecker civicLastReviewedChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("timestamp", true);
        map.put("user", true);
        return new JsonDatamodelChecker("CivicLastReviewed", map);
    }

    @NotNull
    static JsonDatamodelChecker civicUserChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("username", true);
        map.put("name", true);
        map.put("display_name", true);
        map.put("role", true);
        map.put("organization", true);
        map.put("affiliation", true);
        map.put("featured_expert", true);
        map.put("area_of_expertise", true);
        map.put("bio", true);
        map.put("url", true);
        map.put("created_at", true);
        map.put("last_seen_at", true);
        map.put("avatars", true);
        map.put("avatar_url", true);
        map.put("twitter_handle", true);
        map.put("facebook_profile", true);
        map.put("linkedin_profile", true);
        map.put("orcid", true);
        map.put("signup_complete", true);
        map.put("accepted_license", true);
        map.put("id", true);
        return new JsonDatamodelChecker("CivicUser", map);
    }

    @NotNull
    static JsonDatamodelChecker civicOrganizationChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("name", false);
        map.put("url", false);
        map.put("profile_image", false);
        map.put("id", false);
        map.put("description", false);
        return new JsonDatamodelChecker("CivicOrganization", map);
    }

    @NotNull
    static JsonDatamodelChecker civicProfileImageChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("x14", true);
        map.put("x32", true);
        map.put("x64", true);
        map.put("x128", true);
        map.put("x256", true);
        return new JsonDatamodelChecker("CivicProfileImage", map);
    }

    @NotNull
    static JsonDatamodelChecker civicAvatarsChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("x14", true);
        map.put("x32", true);
        map.put("x64", true);
        map.put("x128", true);
        return new JsonDatamodelChecker("CivicAvatars", map);
    }

    @NotNull
    static JsonDatamodelChecker jaxEntryChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("responseType", true);
        map.put("approvalStatus", true);
        map.put("molecularProfile", true);
        map.put("therapy", true);
        map.put("evidenceType", true);
        map.put("indication", true);
        map.put("efficacyEvidence", true);
        map.put("references", true);
        map.put("id", true);
        return new JsonDatamodelChecker("JaxEntry", map);
    }

    @NotNull
    static JsonDatamodelChecker jaxMolecularProfileChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("profileName", true);
        map.put("id", true);
        return new JsonDatamodelChecker("JaxMolecularProfile", map);
    }

    @NotNull
    static JsonDatamodelChecker jaxTherapyChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("therapyName", true);
        map.put("id", true);
        return new JsonDatamodelChecker("JaxTherapy", map);
    }

    @NotNull
    static JsonDatamodelChecker jaxIndicationChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("source", true);
        map.put("id", true);
        map.put("name", true);
        return new JsonDatamodelChecker("JaxIndication", map);
    }

    @NotNull
    static JsonDatamodelChecker jaxReferenceChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("url", true);
        map.put("id", true);
        map.put("pubMedId", true);
        map.put("title", true);
        return new JsonDatamodelChecker("JaxReference", map);
    }

    @NotNull
    static JsonDatamodelChecker oncoKbEntryChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("biological", false);
        map.put("clinical", false);
        return new JsonDatamodelChecker("OncoKbEntry", map);
    }

    @NotNull
    static JsonDatamodelChecker oncoKbBiologicalChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("mutationEffectPmids", true);
        map.put("Isoform", true);
        map.put("variant", true);
        map.put("Entrez Gene ID", true);
        map.put("oncogenic", true);
        map.put("mutationEffect", true);
        map.put("RefSeq", true);
        map.put("gene", true);
        map.put("mutationEffectAbstracts", true);
        return new JsonDatamodelChecker("OncoKbBiological", map);
    }

    @NotNull
    static JsonDatamodelChecker oncoKbClinicalChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("RefSeq", true);
        map.put("level", true);
        map.put("Isoform", true);
        map.put("variant", true);
        map.put("Entrez Gene ID", true);
        map.put("drugPmids", true);
        map.put("cancerType", true);
        map.put("drug", true);
        map.put("gene", true);
        map.put("level_label", true);
        map.put("drugAbstracts", true);
        return new JsonDatamodelChecker("OncoKbClinical", map);
    }

    @NotNull
    static JsonDatamodelChecker oncoKbDrugsAbstractChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("text", true);
        map.put("link", true);
        return new JsonDatamodelChecker("OncoKbDrugsAbstract", map);
    }

    @NotNull
    static JsonDatamodelChecker oncoKbVariantChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("variantResidues", true);
        map.put("proteinStart", true);
        map.put("name", true);
        map.put("proteinEnd", true);
        map.put("refResidues", true);
        map.put("alteration", true);
        map.put("consequence", true);
        map.put("gene", true);
        return new JsonDatamodelChecker("OncoKbVariant", map);
    }

    @NotNull
    static JsonDatamodelChecker oncoKbConsequenceChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("term", true);
        map.put("description", true);
        map.put("isGenerallyTruncating", true);
        return new JsonDatamodelChecker("OncoKbConsequence", map);
    }

    @NotNull
    static JsonDatamodelChecker oncoKbGeneChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("oncogene", true);
        map.put("name", true);
        map.put("hugoSymbol", true);
        map.put("curatedRefSeq", true);
        map.put("entrezGeneId", true);
        map.put("geneAliases", true);
        map.put("tsg", true);
        map.put("curatedIsoform", true);
        return new JsonDatamodelChecker("OncoKbGene", map);
    }

    @NotNull
    static JsonDatamodelChecker pmkbEntryChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("tumor", true);
        map.put("tissues", true);
        map.put("variant", true);
        return new JsonDatamodelChecker("PmkbEntry", map);
    }

    @NotNull
    static JsonDatamodelChecker pmkbTumorChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("name", true);
        map.put("id", true);
        return new JsonDatamodelChecker("PmkbTumor", map);
    }

    @NotNull
    static JsonDatamodelChecker pmkbTissueChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("name", true);
        map.put("id", true);
        return new JsonDatamodelChecker("PmkbTissue", map);
    }

    @NotNull
    static JsonDatamodelChecker pmkbVariantChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("name", true);
        map.put("coordinates", true);
        map.put("chromosome", true);
        map.put("cytoband", true);
        map.put("gene", true);
        map.put("transcript", true);
        map.put("effect", true);
        map.put("codons", true);
        map.put("exons", true);
        map.put("dna_change", true);
        map.put("amino_acid_change", true);
        map.put("germline", true);
        map.put("partner_gene", true);
        map.put("cnv_type", true);
        map.put("chromosome_based_cnv", true);
        map.put("variant_type", true);
        map.put("cosmic", true);
        map.put("description", true);
        map.put("description_type", true);
        map.put("notes", true);
        map.put("id", true);
        return new JsonDatamodelChecker("PmkbVariant", map);
    }

    @NotNull
    static JsonDatamodelChecker pmkbGeneChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("name", true);
        map.put("created_at", true);
        map.put("updated_at", true);
        map.put("active_ind", true);
        map.put("description", true);
        map.put("external_id", true);
        map.put("id", true);
        return new JsonDatamodelChecker("PmkbGene", map);
    }

    @NotNull
    static JsonDatamodelChecker cgiEntryChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("Gene", true);
        map.put("Biomarker", true);
        map.put("Alteration", true);
        map.put("Alteration type", true);
        map.put("transcript", true);
        map.put("individual_mutation", true);
        map.put("gDNA", true);
        map.put("cDNA", true);
        map.put("info", true);
        map.put("region", true);
        map.put("strand", true);
        map.put("Association", true);
        map.put("Drug", true);
        map.put("Drug family", true);
        map.put("Drug full name", true);
        map.put("Drug status", true);
        map.put("Targeting", true);
        map.put("Primary Tumor type", true);
        map.put("Metastatic Tumor Type", true);
        map.put("Evidence level", true);
        map.put("Source", true);
        map.put("Curator", true);
        map.put("Assay type", true);
        return new JsonDatamodelChecker("CgiEntry", map);
    }

    @NotNull
    static JsonDatamodelChecker jaxTrialsEntryChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("nctId", true);
        map.put("title", true);
        map.put("variantRequirements", true);
        map.put("variantRequirementDetails", true);
        map.put("indications", true);
        map.put("therapies", true);
        map.put("gender", true);
        map.put("recruitment", true);
        map.put("phase", true);
        map.put("sponsors", true);
        map.put("updateDate", true);
        return new JsonDatamodelChecker("JaxTrialsEntry", map);
    }

    @NotNull
    static JsonDatamodelChecker jaxTrialsVariantRequirementDetailsChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("molecularProfile", true);
        map.put("requirementType", true);
        return new JsonDatamodelChecker("JaxTrialsVariantRequirementDetails", map);
    }

    @NotNull
    static JsonDatamodelChecker jaxTrialsMolecularProfileChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("profileName", true);
        map.put("id", true);
        return new JsonDatamodelChecker("JaxTrialsMolecularProfile", map);
    }

    @NotNull
    static JsonDatamodelChecker jaxTrialsIndicationChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("name", true);
        map.put("source", true);
        map.put("id", true);
        return new JsonDatamodelChecker("JaxTrialsIndication", map);
    }

    @NotNull
    static JsonDatamodelChecker jaxTrialsTherapyChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("therapyName", true);
        map.put("id", true);
        return new JsonDatamodelChecker("JaxTrialsTherapy", map);
    }

    @NotNull
    static JsonDatamodelChecker molecularMatchEntryChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("direction", true);
        map.put("biomarkerClass", true);
        map.put("mutations", true);
        map.put("variantInfo", true);
        map.put("prevalence", true);
        map.put("_score", true);
        map.put("sources", true);
        map.put("clinicalSignificance", true);
        map.put("tier", true);
        map.put("tierExplanation", true);
        map.put("ampcap", true);
        map.put("civic", true);
        map.put("regulatoryBody", true);
        map.put("regulatoryBodyApproved", true);
        map.put("guidelineBody", false);
        map.put("guidelineVersion", false);
        map.put("includeGene1", false);
        map.put("includeFinding1", false);
        map.put("includeCondition1", true);
        map.put("includeMutation1", false);
        map.put("includeDrug1", false);
        map.put("includeDrugclass1", false);
        map.put("includeResistance1", false);
        map.put("includeStage0", false);
        map.put("includeGene0", false);
        map.put("includeCondition0", true);
        map.put("includeMutation0", false);
        map.put("criteriaMet", true);
        map.put("criteriaUnmet", true);
        map.put("ast", true);
        map.put("institution", false);
        map.put("tags", true);
        map.put("classifications", true);
        map.put("noTherapyAvailable", false);
        map.put("therapeuticContext", true);
        map.put("sixtier", true);
        map.put("mvld", true);
        map.put("autoGenerateNarrative", true);
        map.put("narrative", true);
        map.put("expression", true);
        map.put("customer", true);
        map.put("version", true);
        map.put("id", true);
        map.put("external_id", false);
        map.put("uniqueKey", true);
        map.put("hashKey", true);
        return new JsonDatamodelChecker("MolecularMatchEntry", map);
    }

    @NotNull
    static JsonDatamodelChecker molecularMatchMutationChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("geneSymbol", true);
        map.put("name", true);
        map.put("transcriptRecognized", false);
        map.put("transcript", false);
        map.put("longestTranscript", false);
        map.put("uniprotTranscript", false);
        map.put("transcriptConsequence", false);
        map.put("parents", true);
        map.put("wgsaData", false);
        map.put("wgsaMap", false);
        map.put("exonsInfo", false);
        map.put("fusionData", false);
        map.put("mutation_type", true);
        map.put("sources", true);
        map.put("synonyms", true);
        map.put("GRCh37_location", true);
        map.put("pathology", true);
        map.put("cdna", true);
        map.put("description", true);
        map.put("_src", true);
        map.put("id", true);
        return new JsonDatamodelChecker("MolecularMatchMutation", map);
    }

    @NotNull
    static JsonDatamodelChecker molecularMatchTranscriptConsequenceChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("chr", false);
        map.put("start", false);
        map.put("stop", false);
        map.put("ref", false);
        map.put("alt", false);
        map.put("referenceGenome", true);
        map.put("transcript", true);
        map.put("strand", true);
        map.put("cdna", false);
        map.put("amino_acid_change", false);
        map.put("intronNumber", true);
        map.put("exonNumber", false);
        map.put("suppress", true);
        map.put("custom", true);
        map.put("validated", true);
        map.put("compositeKey", true);
        return new JsonDatamodelChecker("MolecularMatchTranscriptConsequence", map);
    }

    @NotNull
    static JsonDatamodelChecker molecularMatchWGSADataChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("locations", true);
        return new JsonDatamodelChecker("MolecularMatchWGSAData", map);
    }

    @NotNull
    static JsonDatamodelChecker molecularMatchWGSALocationChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("Gene", true);
        map.put("Chr", true);
        map.put("Start", true);
        map.put("End", true);
        map.put("Ref", true);
        map.put("Alt", true);
        map.put("Chr_Start_Ref_Alt", true);
        map.put("Transcript", true);
        map.put("NucleotideChange", true);
        map.put("AA", false);
        map.put("FullAA", true);
        map.put("ExonicFunc", false);
        map.put("PopFreqMax", true);
        map.put("ClinVar_DIS", false);
        map.put("ClinVar_SIG", false);
        map.put("ClinVar_STATUS", false);
        map.put("ClinVar_DBID", false);
        map.put("ExAC_AFR", false);
        map.put("ExAC_AMR", false);
        map.put("ExAC_EAS", false);
        map.put("ExAC_FIN", false);
        map.put("ExAC_NFE", false);
        map.put("ExAC_SAS", false);
        map.put("ExAC_Freq", false);
        map.put("1000G_AFR", false);
        map.put("1000G_AMR", false);
        map.put("1000G_EUR", false);
        map.put("1000G_EAS", false);
        map.put("1000G_SAS", false);
        map.put("1000G_ALL", false);
        map.put("FATHMM", true);
        map.put("FATHMM_Pred", true);
        map.put("ESP6500si_AA", false);
        map.put("ESP6500si_EA", false);
        map.put("dbSNP", false);
        map.put("COSMIC_ID", false);
        map.put("phyloP46way_placental", true);
        map.put("phyloP100way_vertebrate", true);
        map.put("SiPhy_29way_logOdds", true);
        map.put("GWAS_SNP", false);
        map.put("GWAS_DIS", false);
        map.put("GWAS_PUBMED", false);
        map.put("GERP++_RS", true);
        map.put("Func", true);
        map.put("wgRna", false);
        map.put("targetScanS", false);
        map.put("_key", true);
        return new JsonDatamodelChecker("MolecularMatchWGSALocation", map);
    }

    @NotNull
    static JsonDatamodelChecker molecularMatchWGSAMapChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("name", true);
        map.put("Gene", true);
        map.put("Transcript", true);
        map.put("Exon", false);
        map.put("GRCh37_Chr_Start_Ref_Alt", true);
        map.put("NucleotideChange", true);
        map.put("AA", false);
        map.put("Synonyms", true);
        map.put("ProtCoords", true);
        return new JsonDatamodelChecker("MolecularMatchWGSAMap", map);
    }

    @NotNull
    static JsonDatamodelChecker molecularMatchExonsInfoChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("chr", true);
        map.put("transcript", true);
        map.put("txStart", false);
        map.put("txEnd", false);
        map.put("cdsStart", false);
        map.put("cdsEnd", false);
        map.put("exonBoundaries", true);
        return new JsonDatamodelChecker("MolecularMatchExonsInfo", map);
    }

    @NotNull
    static JsonDatamodelChecker molecularMatchExonBoundariesChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("1", false);
        map.put("2", false);
        map.put("3", false);
        map.put("4", false);
        map.put("5", false);
        map.put("6", false);
        map.put("7", false);
        map.put("8", false);
        map.put("9", false);
        map.put("10", false);
        map.put("11", false);
        map.put("12", false);
        map.put("13", false);
        map.put("14", false);
        map.put("15", false);
        map.put("16", false);
        map.put("17", false);
        map.put("18", false);
        map.put("19", false);
        map.put("20", false);
        map.put("21", false);
        map.put("22", false);
        map.put("23", false);
        map.put("24", false);
        map.put("25", false);
        map.put("26", false);
        map.put("27", false);
        map.put("28", false);
        map.put("29", false);
        map.put("30", false);
        map.put("31", false);
        map.put("32", false);
        map.put("33", false);
        map.put("34", false);
        map.put("35", false);
        map.put("36", false);
        map.put("37", false);
        map.put("38", false);
        map.put("39", false);
        map.put("40", false);
        map.put("41", false);
        return new JsonDatamodelChecker("MolecularMatchExonBoundaries", map);
    }

    @NotNull
    static JsonDatamodelChecker molecularMatchPositionChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("start", true);
        map.put("stop", true);
        return new JsonDatamodelChecker("MolecularMatchPosition", map);
    }

    @NotNull
    static JsonDatamodelChecker molecularMatchFusionDataChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("source", false);
        map.put("synonym", false);
        map.put("Achr", false);
        map.put("Aband", false);
        map.put("Agene", false);
        map.put("Acoord", false);
        map.put("Atx", false);
        map.put("Aori", false);
        map.put("Agreg", false);
        map.put("Bchr", false);
        map.put("Bband", false);
        map.put("Bgene", false);
        map.put("Bcoord", false);
        map.put("Btx", false);
        map.put("Bori", false);
        map.put("Bgreg", false);
        map.put("ins", false);
        map.put("Paper", false);
        return new JsonDatamodelChecker("MolecularMatchFusionData", map);
    }

    @NotNull
    static JsonDatamodelChecker molecularMatchFusionGenomicRegionChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("num", true);
        map.put("type", true);
        return new JsonDatamodelChecker("MolecularMatchFusionGenomicRegion", map);
    }

    @NotNull
    static JsonDatamodelChecker molecularMatchGRCh37LocationChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("chr", true);
        map.put("start", true);
        map.put("stop", true);
        map.put("ref", true);
        map.put("alt", true);
        map.put("strand", true);
        map.put("transcript_consequences", true);
        map.put("validated", true);
        map.put("compositeKey", true);
        return new JsonDatamodelChecker("MolecularMatchGRCh37Location", map);
    }

    @NotNull
    static JsonDatamodelChecker molecularMatchGRCh37TranscriptConsequenceChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("transcript", true);
        map.put("cdna", true);
        map.put("amino_acid_change", true);
        map.put("txSites", true);
        map.put("intronNumber", true);
        map.put("exonNumber", true);

        return new JsonDatamodelChecker("MolecularMatchGRCh37TranscriptConsequence", map);
    }

    @NotNull
    static JsonDatamodelChecker molecularMatchVariantInfoChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("name", true);
        map.put("gene", true);
        map.put("transcript", true);
        map.put("classification", true);
        map.put("consequences", true);
        map.put("fusions", true);
        map.put("locations", true);
        map.put("geneFusionPartner", true);
        map.put("COSMIC_ID", true);
        map.put("popFreqMax", true);
        return new JsonDatamodelChecker("MolecularMatchVariantInfo", map);
    }

    @NotNull
    static JsonDatamodelChecker molecularMatchFusionChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("chr", true);
        map.put("referenceGenome", true);
        map.put("LBPWREP", true);
        map.put("LBPWLEP", true);
        map.put("RBPWREP", true);
        map.put("RBPWLEP", true);
        map.put("intronNumber", true);
        map.put("exonNumber", true);
        return new JsonDatamodelChecker("MolecularMatchFusion", map);
    }

    @NotNull
    static JsonDatamodelChecker molecularMatchLocationChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("chr", true);
        map.put("start", true);
        map.put("stop", true);
        map.put("ref", false);
        map.put("alt", false);
        map.put("cdna", false);
        map.put("amino_acid_change", false);
        map.put("referenceGenome", false);
        map.put("strand", false);
        map.put("intronNumber", false);
        map.put("exonNumber", false);
        return new JsonDatamodelChecker("MolecularMatchLocation", map);
    }

    @NotNull
    static JsonDatamodelChecker molecularMatchPrevalenceChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("studyId", true);
        map.put("count", true);
        map.put("samples", true);
        map.put("percent", true);
        map.put("molecular", false);
        map.put("condition", false);
        return new JsonDatamodelChecker("MolecularMatchPrevalence", map);
    }

    @NotNull
    static JsonDatamodelChecker molecularMatchSourceChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("name", true);
        map.put("type", true);
        map.put("subType", false);
        map.put("valid", true);
        map.put("pubId", true);
        map.put("link", true);
        map.put("trialId", false);
        map.put("trialPhase", false);
        map.put("year", true);
        map.put("functionalConsequence", false);
        map.put("institution", false);
        map.put("trustRating", false);
        map.put("suppress", true);
        map.put("id", true);
        return new JsonDatamodelChecker("MolecularMatchSource", map);
    }

    @NotNull
    static JsonDatamodelChecker molecularMatchTierExplanationChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("tier", true);
        map.put("step", true);
        map.put("message", true);
        map.put("success", true);
        return new JsonDatamodelChecker("MolecularMatchTierExplanation", map);
    }

    @NotNull
    static JsonDatamodelChecker molecularMatchCriteriaUnmetChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("term", true);
        map.put("filterType", true);
        map.put("priority", true);
        map.put("facet", true);
        map.put("valid", false);
        map.put("transcript", false);
        map.put("isNew", false);
        map.put("generatedBy", false);
        map.put("generatedByTerm", false);
        map.put("suppress", true);
        map.put("manualSuppress", false);
        map.put("primary", false);
        map.put("compositeKey", true);
        map.put("custom", false);
        return new JsonDatamodelChecker("MolecularMatchCriteriaUnmet", map);
    }

    @NotNull
    static JsonDatamodelChecker molecularMatchAstChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("type", true);
        map.put("raw", false);
        map.put("value", false);
        map.put("operator", false);
        map.put("left", false);
        map.put("right", false);
        return new JsonDatamodelChecker("MolecularMatchAst", map);
    }

    @NotNull
    static JsonDatamodelChecker molecularMatchTagChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("term", true);
        map.put("facet", true);
        map.put("filterType", false);
        map.put("priority", true);
        map.put("transcript", false);
        map.put("valid", false);
        map.put("generatedBy", false);
        map.put("generatedByTerm", false);
        map.put("isNew", false);
        map.put("primary", false);
        map.put("custom", false);
        map.put("suppress", false);
        map.put("manualSuppress", false);
        map.put("composite", false);
        map.put("compositeKey", false);
        return new JsonDatamodelChecker("MolecularMatchTag", map);
    }

    @NotNull
    static JsonDatamodelChecker molecularMatchClassificationChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("name", false);
        map.put("geneSymbol", false);
        map.put("expandGeneSearch", false);
        map.put("transcript", false);
        map.put("transcripts", false);
        map.put("Chr", false);
        map.put("Start", false);
        map.put("End", false);
        map.put("Ref", false);
        map.put("Alt", false);
        map.put("NucleotideChange", false);
        map.put("Exon", false);
        map.put("ExonicFunc", false);
        map.put("classification", true);
        map.put("classificationOverride", false);
        map.put("pathology", false);
        map.put("copyNumberType", false);
        map.put("drugsApprovedOnLabelCount", false);
        map.put("drugsApprovedOffLabelCount", false);
        map.put("drugsExperimentalCount", false);
        map.put("trialCount", false);
        map.put("publicationCount", false);
        map.put("sources", false);
        map.put("dbSNP", false);
        map.put("COSMIC_ID", false);
        map.put("PopFreqMax", false);
        map.put("parents", false);
        map.put("rootTerm", false);
        map.put("alias", false);
        map.put("priority", false);
        map.put("description", false);
        return new JsonDatamodelChecker("MolecularMatchClassification", map);
    }

    @NotNull
    static JsonDatamodelChecker molecularMatchParentChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("name", true);
        map.put("type", false);
        map.put("actionableParent", false);
        map.put("transcripts", true);
        return new JsonDatamodelChecker("MolecularMatchParent", map);
    }

    @NotNull
    static JsonDatamodelChecker molecularMatchTherapeuticContextChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("name", true);
        map.put("facet", true);
        map.put("suppress", true);
        map.put("valid", false);
        return new JsonDatamodelChecker("MolecularMatchTherapeuticContext", map);
    }

    @NotNull
    static JsonDatamodelChecker molecularMatchTrialsEntryChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("status", true);
        map.put("startDate", false);
        map.put("title", true);
        map.put("molecularAlterations", true);
        map.put("_score", true);
        map.put("interventions", true);
        map.put("locations", true);
        map.put("briefTitle", true);
        map.put("overallContact", true);
        map.put("link", true);
        map.put("phase", true);
        map.put("tags", true);
        map.put("id", true);
        map.put("studyType", true);
        return new JsonDatamodelChecker("MolecularMatchTrialsEntry", map);
    }

    @NotNull
    static JsonDatamodelChecker molecularMatchTrialsInterventionChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("intervention_name", false);
        map.put("other_name", false);
        map.put("description", false);
        map.put("arm_group_label", false);
        map.put("intervention_type", false);
        return new JsonDatamodelChecker("MolecularMatchTrialsIntervention", map);
    }

    @NotNull
    static JsonDatamodelChecker molecularMatchTrialsLocationChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("status", true);
        map.put("last_name", false);
        map.put("email", false);
        map.put("phone", false);
        map.put("phone_backup", false);
        map.put("email_backup", false);
        map.put("last_name_backup", false);
        map.put("phone_ext_backup", false);
        map.put("phone_ext", false);
        map.put("city", false);
        map.put("_valid", false);
        map.put("zip", false);
        map.put("created", false);
        map.put("country", false);
        map.put("number", false);
        map.put("id", false);
        map.put("lastUpdated", false);
        map.put("contact", false);
        map.put("state", false);
        map.put("street", false);
        map.put("location", false);
        map.put("po_box", false);
        map.put("failedGeocode", false);
        map.put("geo", false);
        map.put("_validMessage", false);
        map.put("name", false);
        return new JsonDatamodelChecker("MolecularMatchTrialsLocation", map);
    }

    @NotNull
    static JsonDatamodelChecker molecularMatchTrialsOverallContactChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("phone", false);
        map.put("last_name", false);
        map.put("email", false);
        map.put("affiliation", false);
        map.put("phone_ext", false);
        map.put("country", false);
        map.put("city", false);
        map.put("name", false);
        map.put("zip", false);
        map.put("url", false);
        map.put("street", false);
        map.put("type", false);
        return new JsonDatamodelChecker("MolecularMatchTrialsOverallContact", map);
    }

    @NotNull
    static JsonDatamodelChecker molecularMatchTrialsGeoChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("lat", true);
        map.put("lon", true);
        return new JsonDatamodelChecker("MolecularMatchTrialsGeo", map);
    }

    @NotNull
    static JsonDatamodelChecker molecularMatchTrialsSubLocationChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("type", true);
        map.put("coordinates", true);
        return new JsonDatamodelChecker("MolecularMatchTrialsSubLocation", map);
    }

    @NotNull
    static JsonDatamodelChecker molecularMatchTrialsContactChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("phone", false);
        map.put("name", false);
        map.put("email", false);
        return new JsonDatamodelChecker("MolecularMatchTrialsContact", map);
    }

    @NotNull
    static JsonDatamodelChecker molecularMatchTrialsTagChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("facet", true);
        map.put("compositeKey", true);
        map.put("suppress", true);
        map.put("filterType", true);
        map.put("term", true);
        map.put("custom", true);
        map.put("priority", true);
        map.put("alias", false);
        map.put("manualSuppress", false);
        map.put("generatedBy", false);
        map.put("generatedByTerm", false);
        map.put("id", false);
        map.put("manualPriority", false);
        return new JsonDatamodelChecker("MolecularMatchTrialsTag", map);
    }

    @NotNull
    static JsonDatamodelChecker sageEntryChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("gene", true);
        map.put("entrez_id", true);
        map.put("clinical_manifestation", true);
        map.put("response_type", true);
        map.put("evidence_label", true);
        map.put("drug_labels", true);
        map.put("germline_or_somatic", true);
        map.put("publication_url", true);
        return new JsonDatamodelChecker("SageEntry", map);
    }
}
