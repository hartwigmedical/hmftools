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
        return new ViccDatamodelChecker("MolecularMatchTrialsEntryChecker", map);
    }

    @NotNull
    static ViccDatamodelChecker molecularMatchTrialsInterventionChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("intervention_name", false);
        map.put("other_name", false);
        map.put("description", false);
        map.put("arm_group_label", false);
        map.put("intervention_type", false);
        return new ViccDatamodelChecker("MolecularMatchTrialsInterventionChecker", map);
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
        return new ViccDatamodelChecker("MolecularMatchTrialsLocationChecker", map);
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
        return new ViccDatamodelChecker("MolecularMatchTrialsOverallContactChecker", map);
    }

    @NotNull
    static ViccDatamodelChecker molecularMatchTrialsGeoChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("lat", true);
        map.put("lon", true);
        return new ViccDatamodelChecker("MolecularMatchTrialsGeoChecker", map);
    }

    @NotNull
    static ViccDatamodelChecker molecularMatchTrialsSubLocationChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("type", true);
        map.put("coordinates", true);
        return new ViccDatamodelChecker("MolecularMatchTrialsSubLocationChecker", map);
    }

    @NotNull
    static ViccDatamodelChecker molecularMatchTrialsContactChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("phone", false);
        map.put("name", false);
        map.put("email", false);
        return new ViccDatamodelChecker("MolecularMatchTrialsContactChecker", map);
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
        return new ViccDatamodelChecker("MolecularMatchTrialsTagChecker", map);
    }

    @NotNull
    static ViccDatamodelChecker sageChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("gene", true);
        map.put("entrez_id", true);
        map.put("clinical_manifestation", true);
        map.put("response_type", true);
        map.put("evidence_label", true);
        map.put("drug_labels", true);
        map.put("germline_or_somatic", true);
        map.put("publication_url", true);
        return new ViccDatamodelChecker("SageChecker", map);
    }
}
