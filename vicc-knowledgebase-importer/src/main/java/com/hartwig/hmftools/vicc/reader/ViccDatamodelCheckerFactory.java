package com.hartwig.hmftools.vicc.reader;

import java.util.Map;

import com.google.common.collect.Maps;

import org.jetbrains.annotations.NotNull;

final class ViccDatamodelCheckerFactory {

    private ViccDatamodelCheckerFactory() {
    }

    @NotNull
    static ViccDatamodelChecker viccEntryChecker() {
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
        return new ViccDatamodelChecker("ViccEntry", map);
    }

    @NotNull
    static ViccDatamodelChecker geneIdentifierChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("symbol", true);
        map.put("entrez_id", true);
        map.put("ensembl_gene_id", true);
        return new ViccDatamodelChecker("GeneIdentifier", map);
    }

    @NotNull
    static ViccDatamodelChecker featureChecker() {
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
        return new ViccDatamodelChecker("Feature", map);
    }

    @NotNull
    static ViccDatamodelChecker featureInfoChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("germline_or_somatic", true);
        return new ViccDatamodelChecker("FeatureInfo", map);
    }

    @NotNull
    static ViccDatamodelChecker featureAttributeChecker() {
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
        return new ViccDatamodelChecker("FeatureAttribute", map);
    }

    @NotNull
    static ViccDatamodelChecker featureAttributeObjectChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("string_value", true);
        return new ViccDatamodelChecker("FeatureAttributeObject", map);
    }

    @NotNull
    static ViccDatamodelChecker sequenceOntologyChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("hierarchy", false);
        map.put("soid", true);
        map.put("parent_soid", true);
        map.put("name", true);
        map.put("parent_name", true);
        return new ViccDatamodelChecker("SequenceOntology", map);
    }

    @NotNull
    static ViccDatamodelChecker associationChecker() {
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
        return new ViccDatamodelChecker("Association", map);
    }

    @NotNull
    static ViccDatamodelChecker evidenceChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("info", true);
        map.put("evidenceType", true);
        map.put("description", true);
        return new ViccDatamodelChecker("Evidence", map);
    }

    @NotNull
    static ViccDatamodelChecker evidenceInfoChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("publications", true);
        return new ViccDatamodelChecker("EvidenceInfo", map);
    }

    @NotNull
    static ViccDatamodelChecker evidenceTypeChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("sourceName", true);
        map.put("id", false);
        return new ViccDatamodelChecker("EvidenceType", map);
    }

    @NotNull
    static ViccDatamodelChecker environmentalContextChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("term", false);
        map.put("description", true);
        map.put("taxonomy", false);
        map.put("source", false);
        map.put("usan_stem", false);
        map.put("approved_countries", false);
        map.put("toxicity", false);
        map.put("id", false);
        return new ViccDatamodelChecker("EnvironmentalContext", map);
    }

    @NotNull
    static ViccDatamodelChecker taxonomyChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("kingdom", true);
        map.put("direct-parent", true);
        map.put("class", true);
        map.put("subclass", false);
        map.put("superclass", true);
        return new ViccDatamodelChecker("Taxonomy", map);
    }

    @NotNull
    static ViccDatamodelChecker phenotypeChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("type", false);
        map.put("description", true);
        map.put("family", true);
        map.put("id", false);
        return new ViccDatamodelChecker("Phenotype", map);
    }

    @NotNull
    static ViccDatamodelChecker phenotypeTypeChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("source", true);
        map.put("term", true);
        map.put("id", true);
        return new ViccDatamodelChecker("PhenotypeType", map);
    }

    @NotNull
    static ViccDatamodelChecker brcaEntryChecker() {
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

        return new ViccDatamodelChecker("BrcaEntry", map);
    }

    @NotNull
    static ViccDatamodelChecker civicEntryChecker() {
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
        return new ViccDatamodelChecker("CivicEntry", map);
    }

    @NotNull
    static ViccDatamodelChecker civicCoordinatesChecker() {
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
        return new ViccDatamodelChecker("CivicCoordinates", map);
    }

    @NotNull
    static ViccDatamodelChecker civicSourceChecker() {
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
        return new ViccDatamodelChecker("CivicSource", map);
    }

    @NotNull
    static ViccDatamodelChecker civicProvisionalValueChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("description", false);
        return new ViccDatamodelChecker("CivicProvisionalValue", map);
    }

    @NotNull
    static ViccDatamodelChecker civicProvisionalValueDescriptionChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("revision_id", false);
        map.put("value", false);
        return new ViccDatamodelChecker("CivicProvisionalValueDescription", map);
    }

    @NotNull
    static ViccDatamodelChecker civicVariantGroupChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("name", true);
        map.put("type", true);
        map.put("description", true);
        map.put("variants", true);
        map.put("id", true);
        return new ViccDatamodelChecker("CivicVariantGroup", map);
    }

    @NotNull
    static ViccDatamodelChecker civicVariantChecker() {
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
        return new ViccDatamodelChecker("CivicVariant", map);
    }

    @NotNull
    static ViccDatamodelChecker civicVariantTypeChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("name", true);
        map.put("display_name", true);
        map.put("description", true);
        map.put("url", true);
        map.put("so_id", true);
        map.put("id", true);
        return new ViccDatamodelChecker("CivicVariantType", map);
    }

    @NotNull
    static ViccDatamodelChecker civicEvidenceItemChecker() {
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
        return new ViccDatamodelChecker("CivicEvidenceItem", map);
    }

    @NotNull
    static ViccDatamodelChecker civicClinicalTrialChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("name", true);
        map.put("nct_id", true);
        map.put("clinical_trial_url", true);
        map.put("description", true);
        return new ViccDatamodelChecker("CivicClinicalTrial", map);
    }

    @NotNull
    static ViccDatamodelChecker civicPublicationDateChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("year", false);
        map.put("month", false);
        map.put("day", false);
        return new ViccDatamodelChecker("CivicPublicationDate", map);
    }

    @NotNull
    static ViccDatamodelChecker civicDiseaseChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("name", true);
        map.put("display_name", true);
        map.put("doid", true);
        map.put("url", true);
        map.put("id", true);
        return new ViccDatamodelChecker("CivicDisease", map);
    }

    @NotNull
    static ViccDatamodelChecker civicDrugChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("name", true);
        map.put("pubchem_id", true);
        map.put("id", true);
        return new ViccDatamodelChecker("CivicDrug", map);
    }

    @NotNull
    static ViccDatamodelChecker civicLifecycleActionsChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("last_commented_on", false);
        map.put("last_modified", false);
        map.put("last_reviewed", false);
        return new ViccDatamodelChecker("CivicLifecycleActions", map);
    }

    @NotNull
    static ViccDatamodelChecker civicLastCommentedOnChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("timestamp", true);
        map.put("user", true);
        return new ViccDatamodelChecker("CivicLastCommentedOn", map);
    }

    @NotNull
    static ViccDatamodelChecker civicLastModifiedChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("timestamp", true);
        map.put("user", true);
        return new ViccDatamodelChecker("CivicLastModified", map);
    }

    @NotNull
    static ViccDatamodelChecker civicLastReviewedChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("timestamp", true);
        map.put("user", true);
        return new ViccDatamodelChecker("CivicLastReviewed", map);
    }

    @NotNull
    static ViccDatamodelChecker civicUserChecker() {
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
        return new ViccDatamodelChecker("CivicUser", map);
    }

    @NotNull
    static ViccDatamodelChecker civicOrganizationChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("name", false);
        map.put("url", false);
        map.put("profile_image", false);
        map.put("id", false);
        map.put("description", false);
        return new ViccDatamodelChecker("CivicOrganization", map);
    }

    @NotNull
    static ViccDatamodelChecker civicProfileImageChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("x14", true);
        map.put("x32", true);
        map.put("x64", true);
        map.put("x128", true);
        map.put("x256", true);
        return new ViccDatamodelChecker("CivicProfileImage", map);
    }

    @NotNull
    static ViccDatamodelChecker civicAvatarsChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("x14", true);
        map.put("x32", true);
        map.put("x64", true);
        map.put("x128", true);
        return new ViccDatamodelChecker("CivicAvatars", map);
    }

    @NotNull
    static ViccDatamodelChecker jaxEntryChecker() {
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
        return new ViccDatamodelChecker("JaxEntry", map);
    }

    @NotNull
    static ViccDatamodelChecker jaxMolecularProfileChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("profileName", true);
        map.put("id", true);
        return new ViccDatamodelChecker("JaxMolecularProfile", map);
    }

    @NotNull
    static ViccDatamodelChecker jaxTherapyChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("therapyName", true);
        map.put("id", true);
        return new ViccDatamodelChecker("JaxTherapy", map);
    }

    @NotNull
    static ViccDatamodelChecker jaxIndicationChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("source", true);
        map.put("id", true);
        map.put("name", true);
        return new ViccDatamodelChecker("JaxIndication", map);
    }

    @NotNull
    static ViccDatamodelChecker jaxReferenceChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("url", true);
        map.put("id", true);
        map.put("pubMedId", true);
        map.put("title", true);
        return new ViccDatamodelChecker("JaxReference", map);
    }

    @NotNull
    static ViccDatamodelChecker oncoKbEntryChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("biological", false);
        map.put("clinical", false);
        return new ViccDatamodelChecker("OncoKbEntry", map);
    }

    @NotNull
    static ViccDatamodelChecker oncoKbBiologicalChecker() {
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
        return new ViccDatamodelChecker("OncoKbBiological", map);
    }

    @NotNull
    static ViccDatamodelChecker oncoKbClinicalChecker() {
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
        return new ViccDatamodelChecker("OncoKbClinical", map);
    }

    @NotNull
    static ViccDatamodelChecker oncoKbDrugsAbstractChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("text", true);
        map.put("link", true);
        return new ViccDatamodelChecker("OncoKbDrugsAbstract", map);
    }

    @NotNull
    static ViccDatamodelChecker oncoKbVariantChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("variantResidues", true);
        map.put("proteinStart", true);
        map.put("name", true);
        map.put("proteinEnd", true);
        map.put("refResidues", true);
        map.put("alteration", true);
        map.put("consequence", true);
        map.put("gene", true);
        return new ViccDatamodelChecker("OncoKbVariant", map);
    }

    @NotNull
    static ViccDatamodelChecker oncoKbConsequenceChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("term", true);
        map.put("description", true);
        map.put("isGenerallyTruncating", true);
        return new ViccDatamodelChecker("OncoKbConsequence", map);
    }

    @NotNull
    static ViccDatamodelChecker oncoKbGeneChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("oncogene", true);
        map.put("name", true);
        map.put("hugoSymbol", true);
        map.put("curatedRefSeq", true);
        map.put("entrezGeneId", true);
        map.put("geneAliases", true);
        map.put("tsg", true);
        map.put("curatedIsoform", true);
        return new ViccDatamodelChecker("OncoKbGene", map);
    }

    @NotNull
    static ViccDatamodelChecker molecularMatchEntryChecker() {
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
        map.put("guidelineBody", true);
        map.put("guidelineVersion", true);
        map.put("includeCondition1", false);
        map.put("includeMutation1", false);
        map.put("includeDrug1", false);
        map.put("includeStage0", false);
        map.put("includeDrug0", false);
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
        return new ViccDatamodelChecker("MolecularMatchEntry", map);
    }

    @NotNull
    static ViccDatamodelChecker molecularMatchMutationChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("geneSymbol", true);
        map.put("name", true);
        map.put("transcriptRecognized", false);
        map.put("transcript", true);
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
        return new ViccDatamodelChecker("MolecularMatchMutation", map);
    }

    @NotNull
    static ViccDatamodelChecker molecularMatchTranscriptConsequenceChecker() {
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
        return new ViccDatamodelChecker("MolecularMatchTranscriptConsequence", map);
    }

    @NotNull
    static ViccDatamodelChecker molecularMatchWGSADataChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("locations", true);
        return new ViccDatamodelChecker("MolecularMatchWGSAData", map);
    }
    @NotNull
    static ViccDatamodelChecker molecularMatchWGSALocationChecker() {
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
        return new ViccDatamodelChecker("MolecularMatchWGSALocation", map);
    }

    @NotNull
    static ViccDatamodelChecker molecularMatchTrialsEntryChecker() {
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
        return new ViccDatamodelChecker("MolecularMatchTrialsEntry", map);
    }

    @NotNull
    static ViccDatamodelChecker molecularMatchTrialsInterventionChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("intervention_name", false);
        map.put("other_name", false);
        map.put("description", false);
        map.put("arm_group_label", false);
        map.put("intervention_type", false);
        return new ViccDatamodelChecker("MolecularMatchTrialsIntervention", map);
    }

    @NotNull
    static ViccDatamodelChecker molecularMatchTrialsLocationChecker() {
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
        return new ViccDatamodelChecker("MolecularMatchTrialsLocation", map);
    }

    @NotNull
    static ViccDatamodelChecker molecularMatchTrialsOverallContactChecker() {
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
        return new ViccDatamodelChecker("MolecularMatchTrialsOverallContact", map);
    }

    @NotNull
    static ViccDatamodelChecker molecularMatchTrialsGeoChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("lat", true);
        map.put("lon", true);
        return new ViccDatamodelChecker("MolecularMatchTrialsGeo", map);
    }

    @NotNull
    static ViccDatamodelChecker molecularMatchTrialsSubLocationChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("type", true);
        map.put("coordinates", true);
        return new ViccDatamodelChecker("MolecularMatchTrialsSubLocation", map);
    }

    @NotNull
    static ViccDatamodelChecker molecularMatchTrialsContactChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("phone", false);
        map.put("name", false);
        map.put("email", false);
        return new ViccDatamodelChecker("MolecularMatchTrialsContact", map);
    }

    @NotNull
    static ViccDatamodelChecker molecularMatchTrialsTagChecker() {
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
        return new ViccDatamodelChecker("MolecularMatchTrialsTag", map);
    }

    @NotNull
    static ViccDatamodelChecker sageEntryChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("gene", true);
        map.put("entrez_id", true);
        map.put("clinical_manifestation", true);
        map.put("response_type", true);
        map.put("evidence_label", true);
        map.put("drug_labels", true);
        map.put("germline_or_somatic", true);
        map.put("publication_url", true);
        return new ViccDatamodelChecker("SageEntry", map);
    }
}
