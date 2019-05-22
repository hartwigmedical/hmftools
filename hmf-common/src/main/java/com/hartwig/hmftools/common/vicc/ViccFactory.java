package com.hartwig.hmftools.common.vicc;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;

import com.google.common.collect.Lists;
import com.google.gson.JsonArray;
import com.google.gson.JsonElement;
import com.google.gson.JsonObject;
import com.google.gson.JsonParser;
import com.google.gson.stream.JsonReader;
import com.google.gson.stream.JsonToken;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public abstract class ViccFactory {
    private static final Logger LOGGER = LogManager.getLogger(ViccFactory.class);

    private ViccFactory() {
    }

    public static void extractAllFileSpecificFields(@NotNull String allJsonPath) throws IOException {
        final String csvFileName = "/Users/liekeschoenmaker/hmf/tmp/all.csv";
        PrintWriter writer = new PrintWriter(new File(csvFileName));
        JsonParser parser = new JsonParser();
        JsonReader reader = new JsonReader(new FileReader(allJsonPath));
        reader.setLenient(true);
        int index = 1;
        StringBuilder headerCSV = new StringBuilder();
        String headerSource = "source;";
        String headerGenes = "genes;";
        String header = "Source;Primary Tumor type;Drug full name;"
                + "Drug family;Alteration;Biomarker;Gene;Evidence level;Association;";

        headerCSV.append(headerSource);
        headerCSV.append(headerGenes);
        headerCSV.append(header);
        writer.append(headerCSV);
        writer.append("\n");

        while (reader.peek() != JsonToken.END_DOCUMENT && index < 1100) {
            JsonObject object = parser.parse(reader).getAsJsonObject();
            StringBuilder stringToCSVAll = new StringBuilder();

            StringBuilder stringToCSVSource = source.readObjectSource(object);
            StringBuilder stringToCSVGenes = genes.readObjectGenes(object);
            StringBuilder stringToCSVCGI = cgi.readObjectCGISpecificFields(object);
            StringBuilder stringToCSVsage = sage.readObjectSageSpecificFields(object);
            stringToCSVAll.append(stringToCSVSource);
            stringToCSVAll.append(stringToCSVGenes);
            stringToCSVAll.append(stringToCSVCGI);
            stringToCSVAll.append(stringToCSVsage);
            writer.append(stringToCSVAll);
            writer.append("\n");
            index++;
        }
        reader.close();
        writer.close();

    }

    public static void extractAllFile(@NotNull String allJsonPath) throws IOException {
        final String csvFileName = "/Users/liekeschoenmaker/hmf/tmp/all.csv";
        PrintWriter writer = new PrintWriter(new File(csvFileName));
        JsonParser parser = new JsonParser();
        JsonReader reader = new JsonReader(new FileReader(allJsonPath));
        reader.setLenient(true);
        int index = 1;
        StringBuilder headerCSV = new StringBuilder();

        String headerIndex = "index;";
        String headerSource = "source;";
        String headerGenes = "genes;";
        String headerTags = "tags;";
        String headerDevTags = "dev_tags;";
        String headerGeneIdentifiers = "gene_identifiers.Symbol;gene_identifiers.entrez_id;gene_identifiers.ensembl_gene_id;";
        String headerSage = "sage.entrez_id;sage.clinical_manifestation;sage.publication_url;sage.germline_or_somatic;sage.evidence_label;"
                + "sage.drug_labels;sage.response_type;sage.gene;";
        String headerPmkb = "pmkb.tumor.id;pmkb.tumor.name;pmkb.tissues.id;pmkb.tissues.name;pmkb.variant.amino_acid_change;"
                + "pmkb.variant.germline;pmkb.variant.partner_gene;pmkb.variant.codons;pmkb.variant.description;pmkb.variant.exons;"
                + "pmkb.variant.notes;pmkb.variant.cosmic;pmkb.variant.effect;pmkb.variant.cnv_type;pmkb.variant.id;"
                + "pmkb.variant.cytoband;pmkb.variant.variant_type;pmkb.variant.dna_change;pmkb.variant.coordinates;"
                + "pmkb.variant.chromosome_based_cnv;pmkb.variant.gene.description;pmkb.variant.gene.created_at;"
                + "pmkb.variant.gene.updated_at;pmkb.variant.gene.active_ind;pmkb.variant.gene.external_id;"
                + "pmkb.variant.gene.id;pmkb.variant.gene.name;pmkb.transcript;pmkb.description_type;pmkb.chromosome;pmkb.name;";
        String headerBRCA = "brca.Variant_frequency_LOVD;brca.Allele_frequency_FIN_ExAC;brca.ClinVarAccession_ENIGMA;"
                + "brca.Homozygous_count_AFR_ExAC;brca.BX_ID_ExAC;brca.Variant_in_LOVD;brca.Allele_frequency_AFR_ExAC;brca.DBID_LOVD;"
                + "brca.Chr;brca.BX_ID_ENIGMA;brca.Co_occurrence_LR_exLOVD;brca.Homozygous_count_EAS_ExAC;brca.Submitter_ClinVar;"
                + "brca.Allele_frequency_EAS_ExAC;brca.Hg37_End;brca.Submitters_LOVD;brca.Clinical_classification_BIC;"
                + "brca.Homozygous_count_NFE_ExAC;brca.Allele_count_SAS_ExAC;brca.Method_ClinVar;brca.Allele_count_NFE_ExAC;"
                + "brca.Pathogenicity_all;brca.Germline_or_Somatic_BIC;brca.Homozygous_count_SAS_ExAC;brca.BIC_Nomenclature;"
                + "brca.Assertion_method_ENIGMA;brca.Literature_source_exLOVD;brca.Change_Type_id;brca.Collection_method_ENIGMA;"
                + "brca.Sum_family_LR_exLOVD;brca.HGVS_cDNA_LOVD;brca.Homozygous_count_FIN_ExAC;brca.EAS_Allele_frequency_1000_Genomes;"
                + "brca.Ethnicity_BIC;brca.Individuals_LOVD;brca.Variant_in_ExAC;brca.URL_ENIGMA;brca.Allele_Origin_ClinVar;"
                + "brca.Allele_frequency_AMR_ExAC;brca.Variant_in_1000_Genomes;brca.AFR_Allele_frequency_1000_Genomes;"
                + "brca.BX_ID_exLOVD;brca.Source;brca.Condition_ID_value_ENIGMA;brca.HGVS_Protein;brca.Ref;brca.Allele_number_AFR_ExAC;"
                + "brca.Allele_count_AFR_ExAC;brca.BX_ID_LOVD;brca.Synonyms;brca.Gene_Symbol;brca.Comment_on_clinical_significance_ENIGMA;"
                + "brca.Missense_analysis_prior_probability_exLOVD;brca.Allele_number_FIN_ExAC;brca.Posterior_probability_exLOVD;"
                + "brca.Polyphen_Score;brca.Reference_Sequence;brca.Allele_count_EAS_ExAC;brca.Hg38_End;brca.HGVS_cDNA;"
                + "brca.Functional_analysis_technique_LOVD;brca.SAS_Allele_frequency_1000_Genomes;brca.RNA_LOVD;"
                + "brca.Combined_prior_probablility_exLOVD;brca.BX_ID_ClinVar;brca.IARC_class_exLOVD;brca.BX_ID_BIC;brca.Sift_Prediction;"
                + "brca.Allele_number_NFE_ExAC;brca.Allele_origin_ENIGMA;brca.Allele_number_OTH_ExAC;brca.Hg36_End;"
                + "brca.Allele_frequency_SAS_ExAC;brca.Date_Last_Updated_ClinVar;brca.Allele_number_EAS_ExAC;"
                + "brca.Allele_frequency_OTH_ExAC;brca.Source_URL;brca.SCV_ClinVar;brca.Pathogenicity_expert;"
                + "brca.Allele_frequency_1000_Genomes;brca.Functional_analysis_result_LOVD;brca.AMR_Allele_frequency_1000_Genomes;"
                + "brca.Variant_in_ESP;brca.Variant_in_BIC;brca.Clinical_significance_ENIGMA;brca.Max_Allele_Frequency;"
                + "brca.Allele_count_AMR_ExAC;brca.Variant_in_ENIGMA;brca.BX_ID_ESP;brca.Patient_nationality_BIC;brca.BX_ID_1000_Genomes;"
                + "brca.Genomic_Coordinate_hg37;brca.Genomic_Coordinate_hg36;brca.EUR_Allele_frequency_1000_Genomes;"
                + "brca.Number_of_family_member_carrying_mutation_BIC;brca.Segregation_LR_exLOVD;brca.Allele_Frequency;"
                + "brca.Minor_allele_frequency_percent_ESP;brca.Allele_frequency_ExAC;brca.Mutation_type_BIC;"
                + "brca.Assertion_method_citation_ENIGMA;brca.Condition_ID_type_ENIGMA;brca.Allele_count_OTH_ExAC;brca.HGVS_protein_LOVD;"
                + "brca.Variant_in_ClinVar;brca.Clinical_importance_BIC;brca.Discordant;brca.Allele_count_FIN_ExAC;"
                + "brca.Condition_category_ENIGMA;brca.Allele_Frequency_ESP;brca.Homozygous_count_OTH_ExAC;brca.Genetic_origin_LOVD;"
                + "brca.id;brca.Homozygous_count_AMR_ExAC;brca.Clinical_Significance_ClinVar;brca.AA_Allele_Frequency_ESP;"
                + "brca.Protein_Change;brca.Variant_in_exLOVD;brca.EA_Allele_Frequency_ESP;brca.HGVS_RNA;"
                + "brca.Clinical_significance_citations_ENIGMA;brca.Variant_effect_LOVD;brca.Polyphen_Prediction;brca.Data_Release_id;"
                + "brca.Hg37_Start;brca.Hg36_Start;brca.Sift_Score;brca.Genomic_Coordinate_hg38;brca.Alt;brca.Literature_citation_BIC;"
                + "brca.Variant_haplotype_LOVD;brca.Allele_frequency_NFE_ExAC;brca.Hg38_Start;brca.Pos;brca.Date_last_evaluated_ENIGMA;"
                + "brca.Allele_number_SAS_ExAC;brca.Allele_number_AMR_ExAC;";
        String headerCGI = "cgi.Targeting;cgi.Source”;cgi.cDNA;cgi.Primary Tumor type;cgi.individual_mutation;cgi.Drug full name;"
                + "cgi.Curator;cgi.Drug family;cgi.Alteration;cgi.Drug;cgi.Biomarker;cgi.gDNA;cgi.Drug status;cgi.Gene;cgi.transcript;"
                + "cgi.strand;cgi.info;cgi.Assay type;cgi.Alteration type;cgi.region;cgi.Evidence level;cgi.Association;"
                + "cgi.Metastatic Tumor Type;";
        String headerOncokb = "oncokb.biological.mutationEffectPmids;oncokb.biological.Isoform;oncokb.clinical.RefSeq;"
                + "oncokb.clinical.level;oncokb.clinical.Isoform;oncokb.biological.variant.variantResidues;"
                + "oncokb.biological.variant.proteinStart;oncokb.biological.variant.name;oncokb.biological.variant.proteinEnd;"
                + "oncokb.biological.variant.refResidues;oncokb.biological.variant.alteration;oncokb.biological.variant.consequence.term;"
                + "oncokb.biological.variant.consequence.description;oncokb.biological.variant.consequence.isGenerallyTruncating;"
                + "oncokb.biological.variant.gene.oncogene;oncokb.biological.variant.gene.name;oncokb.biological.variant.gene.hugoSymbol;"
                + "oncokb.biological.variant.gene.curatedRefSeq;oncokb.biological.variant.gene.entrezGeneId;"
                + "oncokb.biological.variant.gene.geneAliases;oncokb.biological.variant.gene.tsg;"
                + "oncokb.biological.variant.gene.curatedIsoform;oncokb.biological.Entrez Gene ID;oncokb.biological.oncogenic;"
                + "oncokb.biological.mutationEffect;oncokb.biological.RefSeq;oncokb.biological.gene;"
                + "oncokb.biological.mutationEffectAbstracts;oncokb.clinical.Entrez Gene ID;oncokb.clinical.drugPmids;"
                + "oncokb.clinical.cancerType;oncokb.clinical.drug;oncokb.clinical.drugAbstracts.text;oncokb.clinical.drugAbstracts.link;"
                + "oncokb.clinical.gene;oncokb.clinical.level_label;";
        String headerJax = "jax.responseType;jax.approvalStatus;jax.molecularProfile.profileName;jax.molecularProfile.id;jax.therapy.id;"
                + "jax.therapy.therapyName;jax.evidenceType;jax.indication.source;jax.indication.id;jax.indication.name;"
                + "jax.efficacyEvidence;jax.references.url;jax.references.id;jax.references.pubMedId”;jax.references.title;jax.id;";
        String headerJaxTrials = "jax_trials.indications.source;jax_trials.indications.id;jax_trials.indications.name;jax_trials.title;"
                + "jax_trials.gender;jax_trials.nctId;jax_trials.sponsors;jax_trials.recruitment;jax_trials.variantRequirements;"
                + "jax_trials.updateDate;jax_trials.phase;jax_trials.variantRequirementDetails.molecularProfile.profileName;"
                + "jax_trials.variantRequirementDetails.molecularProfile.id;jax_trials.variantRequirementDetails.requirementType;"
                + "jax_trials.therapies.id;jax_trials.therapies.therapyName;";
        String headerAssociation = "association.drug_labels;association.description;association.publication_url;association.source_link;"
                + "association.variant_name;association.evidence.info.publications;association.evidence.evidenceType.sourceName;"
                + "association.evidence.evidenceType.id;association.evidence.description;" + "association.evidence_label;"
                + "association.phenotype.type.source;association.phenotype.type.term;association.phenotype.type.id;"
                + "association.phenotype.description;association.phenotype.family;association.phenotype.id;association.evidence_level;"
                + "association..response_type;association.environmentalContexts.term;"
                + "association.environmentalContexts.description;association.environmentalContexts.taxonomy.kingdom;"
                + "association.environmentalContexts.taxonomy.direct-parent;association.environmentalContexts.taxonomy.class;"
                + "association.environmentalContexts.taxonomy.subclass;association.environmentalContexts.taxonomy.superclass;"
                + "association.environmentalContexts.source;association.environmentalContexts.usan_stem;association.environmentalContexts.toxicity;"
                + "association.environmentalContexts.approved_countries;association.environmentalContexts.id;";
        String headerFeaturesNames = "feature_names;";
        String headerFeatures = "features.provenance_rule;features.entrez_id;features.end;features.name;features.links;"
                + "features.sequence_ontology.hierarchy;features.sequence_ontology.soid;features.sequence_ontology.parent_soid;"
                + "features.sequence_ontology.name;features.sequence_ontology.parent_name;features.provenance;features.start;"
                + "features.synonyms;features.biomarker_type;features.referenceName;features.geneSymbol;features.alt;"
                + "features.ref;features.chromosome;features.description;";
        String headerCIVIC = "civic.variant_groups;civic.entrez_name;civic.variant_types.display_name;civic.variant_types.description;"
                + "civic.variant_types.url;civic.variant_types.so_id;civic.variant_types.d;civic.variant_types.name;"
                + "civic.civic_actionability_score;civic.clinvar_entries;civic.lifecycle_actions.last_commented_on.timestamp;"
                + "civic.lifecycle_actions.last_commented_on.user.username;"
                + "civic.lifecycle_actions.last_commented_on.user.area_of_expertise;"
                + "civic.lifecycle_actions.last_commented_on.user.organization.url;"
                + "civic.lifecycle_actions.last_commented_on.user.organization.id;"
                + "civic.lifecycle_actions.last_commented_on.user.organization.profile_image;"
                + "civic.lifecycle_actions.last_commented_on.user.organization.x32;"
                + "civic.lifecycle_actions.last_commented_on.user.organization.x256;"
                + "civic.lifecycle_actions.last_commented_on.user.organization.x14;"
                + "civic.lifecycle_actions.last_commented_on.user.organization.x64;"
                + "civic.lifecycle_actions.last_commented_on.user.organization.x128;"
                + "civic.lifecycle_actions.last_commented_on.user.organization.description;"
                + "civic.lifecycle_actions.last_commented_on.user.organization.name;"
                + "civic.lifecycle_actions.last_commented_on.user.twitter_handle;civic.lifecycle_actions.last_commented_on.user.name;"
                + "civic.lifecycle_actions.last_commented_on.user.bio;civic.lifecycle_actions.last_commented_on.user.url;"
                + "civic.lifecycle_actions.last_commented_on.user.created_at;civic.lifecycle_actions.last_commented_on.user.avatarsx32;"
                + "civic.lifecycle_actions.last_commented_on.user.avatarsx14;civic.lifecycle_actions.last_commented_on.user.avatarsx64;"
                + "civic.lifecycle_actions.last_commented_on.user.avatarsx128; "
                + "civic.lifecycle_actions.last_commented_on.user.accepted_license;"
                + "civic.lifecycle_actions.last_commented_on.user.affiliation;civic.lifecycle_actions.last_commented_on.user.avatar_url;"
                + "civic.lifecycle_actions.last_commented_on.user.role;civic.lifecycle_actions.last_commented_on.user.facebook_profile;"
                + "civic.lifecycle_actions.last_commented_on.user.linkedin_profile;civic.lifecycle_actions.last_commented_on.user.orcid;"
                + "civic.lifecycle_actions.last_commented_on.user.display_name;civic.lifecycle_actions.last_commented_on.user.last_seen_at;"
                + "civic.lifecycle_actions.last_commented_on.user.featured_expert;civic.lifecycle_actions.last_commented_on.user.id;"
                + "civic.lifecycle_actions.last_commented_on.user.signup_complete;civic.last_modified.timestamp;"
                + "civic.last_modified.user.username;civic.last_modified.user.area_of_expertise;civic.last_modified.user.organization.url;"
                + "civic.last_modified.user.organization.id;civic.last_modified.user.organization.profile_image.x32;"
                + "civic.last_modified.user.organization.profile_image.x256;civic.last_modified.user.organization.profile_image.x14;"
                + "civic.last_modified.user.organization.profile_image.x64;civic.last_modified.user.organization.profile_image.x128;"
                + "civic.last_modified.user.organization.description;"
                + "civic.last_modified.user.organization.namecivic.last_modified.user.twitter_handle;civic.last_modified.user.name;"
                + "civic.last_modified.user.bio;civic.last_modified.user.url;civic.last_modified.user.created_at;"
                + "civic.last_modified.user.avatars.x32;civic.last_modified.user.avatars.x14;,civic.last_modified.user.avatars.x64;"
                + "civic.last_modified.user.avatars.x128;civic.last_modified.user.accepted_license;"
                + "civic.last_modified.user.affiliation;civic.last_modified.user.avatar_url;civic.last_modified.user.role;"
                + "civic.last_modified.user.facebook_profile;civic.last_modified.user.linkedin_profile;civic.last_modified.user.orcid;"
                + "civic.last_modified.user.display_name;civic.last_modified.user.last_seen_at;civic.last_modified.user.featured_expert;"
                + "civic.last_modified.user.id;civic.last_modified.user.signup_complete;civic.last_reviewed.timestamp;"
                + "civic.last_reviewed.user.username;civic.last_reviewed.user.area_of_expertise;civic.last_reviewed.user.organization.url;"
                + "civic.last_reviewed.user.organization.id;civic.last_reviewed.user.organization.profile_image.x32;"
                + "civic.last_reviewed.user.organization.profile_image.x256;civic.last_reviewed.user.organization.profile_image.x14”;"
                + "civic.last_reviewed.user.organization.profile_image.x64;civic.last_reviewed.user.organization.profile_image.x128;"
                + "civic.last_reviewed.user.organization.description;"
                + "civic.last_reviewed.user.organization..name;civic.last_reviewed.user..twitter_handle;"
                + "civic.last_reviewed.user..name;civic.last_reviewed.user..bio;civic.last_reviewed.user.url;"
                + ",civic.last_reviewed.user.created_at;civic.last_reviewed.user.avatars.x32;civic.last_reviewed.user.avatars.x14;"
                + "civic.last_reviewed.user.avatars.x64;civic.last_reviewed.user.avatars.x128;"
                + "civic.last_reviewed.user..accepted_license;civic.last_reviewed.user..affiliation;civic.last_reviewed.user..avatar_url;"
                + "civic.last_reviewed.user..role;civic.last_reviewed.user..facebook_profile;civic.last_reviewed.user..linkedin_profile;"
                + "civic.last_reviewed.user..orcid;civic.last_reviewed.user..display_name;civic.last_reviewed.user..last_seen_at;"
                + "civic.last_reviewed.user..featured_expert;civic.last_reviewed.user.id;civic.last_reviewed.user..signup_complete;"
                + "civic.variant_aliases;civic.allele_registry_id;civic.provisional_values;civic.gene_id;civic.name;"
                + "civic.evidence_items.status;civic.evidence_items.rating;civic.evidence_items.drug_interaction_type;"
                + "civic.evidence_items.description;civic.evidence_items.open_change_count;civic.evidence_items.evidence_type;"
                + "civic.evidence_items.drugs.pubchem_id;civic.evidence_items.drugs.id;civic.evidence_items.drugs.name;"
                + "civic.evidence_items.variant_origin;civic.evidence_items.disease.doid;civic.evidence_items.disease.url;"
                + "civic.evidence_items.disease.display_name;civic.evidence_items.disease.id;civic.evidence_items.disease.name;"
                + "civic.evidence_items.source.status;civic.evidence_items.source.open_access;civic.evidence_items.source.name;"
                + "civic.evidence_items.source.journal;civic.evidence_items.source.citation;civic.evidence_items.source.pmc_id;"
                + "civic.evidence_items.source.full_journal_title;civic.evidence_items.source.source_url;"
                + "civic.evidence_items.source.clinical_trials;civic.evidence_items.source.pubmed_id;"
                + "civic.evidence_items.source.is_review;civic.evidence_items.source.publication_date.year;"
                + "civic.evidence_items.source.publication_date.day;civic.evidence_items.source.publication_date.month;"
                + "civic.evidence_items.source.id;civic.evidence_items.evidence_direction;civic.evidence_items.variant_id;"
                + "civic.evidence_items.clinical_significance;civic.evidence_items.evidence_level;civic.evidence_items.type;"
                + "civic.evidence_items.id;civic.evidence_items.name;civic.sources;civic.entrez_id;civic.assertions;civic.hgvs_expressions;"
                + "civic.errors;civic.coordinates.chromosome2;civic.coordinates.reference_bases;civic.coordinates.start2;"
                + "civic.coordinates.variant_bases;civic.coordinates.stop;civic.coordinates.stop2;"
                + "civic.coordinates.representative_transcript2;civic.coordinates.start;"
                + "civic.coordinates.representative_transcript;civic.coordinates.ensembl_version;"
                + "civic.coordinates.chromosome;civic.coordinates.reference_build;civic.type;civic.id;civic.description;";
        String headerMolecularMatch = "molecularmatch.criteriaUnmet.priority;molecularmatch.criteriaUnmet.compositeKey;"
                + "molecularmatch.criteriaUnmet.suppress;molecularmatch.criteriaUnmet.filterType;molecularmatch.criteriaUnmet.term;"
                + "molecularmatch.criteriaUnmet.primary;molecularmatch.criteriaUnmet.facet;molecularmatch.criteriaUnmet.valid;"
                + "molecularmatch.criteriaUnmet.custom;molecularmatch.prevalence.count;molecularmatch.prevalence.percent;"
                + "molecularmatch.prevalence.samples;molecularmatch.prevalence.studyId;molecularmatch.prevalence.molecular;"
                + "molecularmatch.prevalence.condition;molecularmatch._score;molecularmatch.autoGenerateNarrative;"
                + "molecularmatch.mutations.transcriptConsequence.amino_acid_change;"
                + "molecularmatch.mutations.transcriptConsequence.compositeKey;molecularmatch.mutations.transcriptConsequence.intronNumber;"
                + "molecularmatch.mutations.transcriptConsequence.exonNumber;molecularmatch.mutations.transcriptConsequence.suppress;"
                + "molecularmatch.mutations.transcriptConsequence.stop;molecularmatch.mutations.transcriptConsequence.custom;"
                + "molecularmatch.mutations.transcriptConsequence.start;molecularmatch.mutations.transcriptConsequence.chr;"
                + "molecularmatch.mutations.transcriptConsequence.strand;molecularmatch.mutations.transcriptConsequence.validated;"
                + "molecularmatch.mutations.transcriptConsequence.transcript;molecularmatch.mutations.transcriptConsequence.cdna;"
                + "molecularmatch.mutations.transcriptConsequence.referenceGenome;molecularmatch.mutations.longestTranscript;"
                + "molecularmatch.mutations.description;molecularmatch.mutations.mutation_type;molecularmatch.mutations._src;"
                + "molecularmatch.mutations.ources;molecularmatch.mutations.synonyms;molecularmatch.mutations.parents;"
                + "molecularmatch.mutations.GRCh37_location.compositeKey;molecularmatch.mutations.GRCh37_location.ref;"
                + "molecularmatch.mutations.GRCh37_location.stop;molecularmatch.mutations.GRCh37_location.start;"
                + "molecularmatch.mutations.GRCh37_location.chr”;molecularmatch.mutations.GRCh37_location.alt;"
                + "molecularmatch.mutations.GRCh37_location.validated;"
                + "molecularmatch.mutations.GRCh37_location.transcript_consequences.amino_acid_change;"
                + "molecularmatch.mutations.GRCh37_location.transcript_consequences.txSites;"
                + "molecularmatch.mutations.GRCh37_location.transcript_consequences.exonNumber;"
                + "molecularmatch.mutations.GRCh37_location.transcript_consequences.intronNumber;"
                + "molecularmatch.mutations.GRCh37_location.transcript_consequences.transcript;"
                + "molecularmatch.mutations.GRCh37_location.transcript_consequences.cdna;molecularmatch.mutations.GRCh37_location.strand;"
                + "molecularmatch.mutations.uniprotTranscript;molecularmatch.mutations.geneSymbol;molecularmatch.mutations.pathology;"
                + "molecularmatch.mutations.transcript;molecularmatch.mutations.id;molecularmatch.mutations.cdna;"
                + "molecularmatch.mutations.name;molecularmatch.sources.name;molecularmatch.sources.suppress;molecularmatch.sources.pubId;"
                + "molecularmatch.sources.subType;molecularmatch.sources.valid;molecularmatch.sources.link;molecularmatch.sources.year;"
                + "molecularmatch.sources.trialId;molecularmatch.sources.type;molecularmatch.sources.id;"
                + "molecularmatch.clinicalSignificance;molecularmatch.id;molecularmatch.includeCondition0;"
                + "molecularmatch.includeCondition1;molecularmatch.uniqueKey;molecularmatch.civic;molecularmatch.hashKey;"
                + "molecularmatch.regulatoryBodyApproved;molecularmatch.version;molecularmatch.includeMutation1;"
                + "molecularmatch.guidelineBody;molecularmatch.regulatoryBody;molecularmatch.customer;molecularmatch.direction;"
                + "molecularmatch.ampcap;molecularmatch.ast.operator;molecularmatch.ast.right.raw;molecularmatch.ast.right.type;"
                + "molecularmatch.ast.right.value;molecularmatch.ast.type;molecularmatch.ast.left.raw;molecularmatch.ast.left.type;"
                + "molecularmatch.ast.left.value;molecularmatch.variantInfo.classification;molecularmatch.variantInfo.name;"
                + "molecularmatch.variantInfo.consequences;molecularmatch.variantInfo.fusions;molecularmatch.variantInfo.locations;"
                + "molecularmatch.variantInfo.locations.amino_acid_change;molecularmatch.variantInfo.locations.intronNumber;"
                + "molecularmatch.variantInfo.locations.exonNumber;molecularmatch.variantInfo.locations.stop;"
                + "molecularmatch.variantInfo.locations.start;molecularmatch.variantInfo.locations.chr;"
                + "molecularmatch.variantInfo.locations.strand;molecularmatch.variantInfo.locations.alt;"
                + "molecularmatch.variantInfo.locations.ref;molecularmatch.variantInfo.locations.cdna;"
                + "molecularmatch.variantInfo.geneFusionPartner;molecularmatch.variantInfo.COSMIC_ID;,molecularmatch.variantInfo.gene;"
                + "molecularmatch.variantInfo.transcript;molecularmatch.variantInfo.popFreqMax;molecularmatch.tier;"
                + "molecularmatch.tierExplanation.tier;molecularmatch.tierExplanation.step;molecularmatch.tierExplanation.message;"
                + "molecularmatch.tierExplanation.success;molecularmatch.mvld;molecularmatch.tags.priority;"
                + "molecularmatch.tags.compositeKey;molecularmatch.tags.suppress;molecularmatch.tags.filterType;"
                + "molecularmatch.tags.term;molecularmatch.tags.primary;molecularmatch.tags.facet;molecularmatch.tags.valid;"
                + "molecularmatch.tags.custom;molecularmatch.criteriaMet;molecularmatch.biomarkerClass;molecularmatch.classifications.End;"
                + "molecularmatch.classifications.classification;molecularmatch.classifications.classificationOverride;"
                + "molecularmatch.classifications.Start;molecularmatch.classifications.Chr;molecularmatch.classifications.geneSymbol;"
                + "molecularmatch.classifications.pathology;molecularmatch.classifications.Ref;molecularmatch.classifications.description;"
                + "molecularmatch.classifications.priority;molecularmatch.classifications.NucleotideChange;"
                + "molecularmatch.classifications.parents;molecularmatch.classifications.drugsExperimentalCount;"
                + "molecularmatch.classifications.Exon;molecularmatch.classifications.drugsApprovedOffLabelCount;"
                + "molecularmatch.classifications.ExonicFunc;molecularmatch.classifications.PopFreqMax;"
                + "molecularmatch.classifications.copyNumberType;molecularmatch.classifications.publicationCount;"
                + "molecularmatch.classifications.transcript;molecularmatch.classifications.dbSNP;molecularmatch.classifications.Alt;"
                + "molecularmatch.classifications.name;molecularmatch.classifications.rootTerm;molecularmatch.classifications.sources;"
                + "molecularmatch.classifications.drugsApprovedOnLabelCount;molecularmatch.classifications.trialCount;"
                + "molecularmatch.classifications.alias;molecularmatch.includeDrug1;molecularmatch.therapeuticContext;"
                + "molecularmatch.sixtier;molecularmatch.narrative;molecularmatch.expression;molecularmatch.includeGene0;";
        String headerMolecularMatchTrials = "molecularmatch_trials.status;molecularmatch_trials.startDate;molecularmatch_trials.title;"
                + "molecularmatch_trials.molecularAlterations;molecularmatch_trials._score;"
                + "molecularmatch_trials.interventions.intervention_name;"
                + "molecularmatch_trials.interventions.other_name;molecularmatch_trials.interventions.description;"
                + "molecularmatch_trials.interventions.arm_group_label;molecularmatch_trials.interventions.iintervention_type;"
                + "molecularmatch_trials.locations.status;"
                + "molecularmatch_trials.locations.city;molecularmatch_trials.locations._valid;molecularmatch_trials.locations.zip;"
                + "molecularmatch_trials.locations.created;molecularmatch_trials.locations.country;molecularmatch_trials.locations.id;"
                + "molecularmatch_trials.locations.lastUpdated;molecularmatch_trials.locations.contact.phone;"
                + "molecularmatch_trials.locations.contact.name;molecularmatch_trials.locations.contact.email;"
                + "molecularmatch_trials.locations.state;molecularmatch_trials.locations.street;"
                + "molecularmatch_trials.locations.location.type;molecularmatch_trials.locations.location.coordinates;"
                + "molecularmatch_trials.locations.po_box;molecularmatch_trials.locations.failedGeocode;"
                + "molecularmatch_trials.locations.geo.lat;molecularmatch_trials.locations.geo.on;"
                + "molecularmatch_trials.locations._validMessage;molecularmatch_trials.locations.name;molecularmatch_trials.briefTitle;"
                + "molecularmatch_trials.overallContact.phone;molecularmatch_trials.overallContact.last_name;"
                + "molecularmatch_trials.overallContact.email;molecularmatch_trials.overallContact.affiliation;"
                + "molecularmatch_trials.link”;molecularmatch_trials.phase;molecularmatch_trials.tags.facet;"
                + "molecularmatch_trials.tags.compositeKey;molecularmatch_trials.tags.suppress;molecularmatch_trials.tags.filterType;"
                + "molecularmatch_trials.tags.term;molecularmatch_trials.tags.custom;molecularmatch_trials.tags.priority;"
                + "molecularmatch_trials.tags.alias;molecularmatch_trials.tags.generatedByTerm;molecularmatch_trials.tags.generatedBy;"
                + "molecularmatch_trials.id;molecularmatch_trials.studyType;";

        headerCSV.append(headerIndex);
        headerCSV.append(headerSource);
        //        headerCSV.append(headerGenes);
        //        headerCSV.append(headerTags);
        //        headerCSV.append(headerDevTags);
        //        headerCSV.append(headerGeneIdentifiers);
        //        headerCSV.append(headerFeaturesNames);
        //        headerCSV.append(headerSage);
        //        headerCSV.append(headerPmkb);
        //        headerCSV.append(headerBRCA);
        //        headerCSV.append(headerCGI);
        //        headerCSV.append(headerOncokb);
        //        headerCSV.append(headerJax);
        //        headerCSV.append(headerJaxTrials);

        //                headerCSV.append(headerAssociation); //TODO check fields for every db
        //        headerCSV.append(headerFeatures); //TODO check fields for every db
        //        headerCSV.append(headerCIVIC); //TODO check fields for every db
        //        headerCSV.append(headerMolecularMatch); //TODO check fields for every db
        headerCSV.append(headerMolecularMatchTrials);//TODO check fields for every db

        writer.append(headerCSV);
        writer.append("\n");

        while (reader.peek() != JsonToken.END_DOCUMENT && index < 5) {
            JsonObject object = parser.parse(reader).getAsJsonObject();
            StringBuilder stringToCSVSource = source.readObjectSource(object);
            LOGGER.info(index + " " + stringToCSVSource);
            StringBuilder stringToCSVAll = new StringBuilder();
            headerCSV.append("index").append(";");
            StringBuilder stringToCSVGenes = genes.readObjectGenes(object);
            StringBuilder stringToCSVTags = tags.readObjectTags(object);
            StringBuilder stringToCSVDevTags = devTags.readObjectDevTags(object);
            StringBuilder stringToCSVGeneIdentifiers = geneIdentifiers.readObjectGeneIdentifiers(object);
            StringBuilder stringToCSVFeaturesNames = featuresNames.readObjectFeaturesNames(object);

            StringBuilder stringToCSVSage = sage.readObjectSage(object);
            StringBuilder stringToCSVPmkb = pmkb.readObjectPmkb(object);
            StringBuilder stringToCSVBrca = brca.readObjectBRCA(object);
            StringBuilder stringToCSVCGI = cgi.readObjectCGI(object);
            StringBuilder stringToCSVOncokb = oncokb.readObjectOncokb(object);
            StringBuilder stringToCSVJax = jax.readObjectJax(object);
            StringBuilder stringToCSVJaxTrials = jaxTrials.readObjectJaxTrials(object);

            StringBuilder stringToCSVAssociation = association.readObjectAssociation(object); //TODO check fields for every db
            StringBuilder stringToCSVFeatures = features.readObjectFeatures(object); //TODO check fields for every db
            StringBuilder stringToCSVCIVIC = civic.readObjectCIVIC(object); //TODO check fields for every db
            StringBuilder stringToCSVMolecularMatch = molecularMatch.readObjectMolecularMatch(object); //TODO check fields for every db
            StringBuilder stringToCSVMolecularMatchTrials =
                    molecularMatchTrial.readObjectMolecularMatchTrials(object); //TODO check fields for every db

            stringToCSVAll.append(index).append(";");
            stringToCSVAll.append(stringToCSVSource);
            //            stringToCSVAll.append(stringToCSVGenes);
            //            stringToCSVAll.append(stringToCSVTags);
            //            stringToCSVAll.append(stringToCSVDevTags);
            //            stringToCSVAll.append(stringToCSVGeneIdentifiers);
            //            stringToCSVAll.append(stringToCSVFeaturesNames);
            //            stringToCSVAll.append(stringToCSVSage);
            //            stringToCSVAll.append(stringToCSVPmkb);
            //            stringToCSVAll.append(stringToCSVBrca);
            //            stringToCSVAll.append(stringToCSVCGI);
            //            stringToCSVAll.append(stringToCSVOncokb);
            //            stringToCSVAll.append(stringToCSVJax);
            //            stringToCSVAll.append(stringToCSVJaxTrials);
            //                        stringToCSVAll.append(stringToCSVAssociation);
            //            stringToCSVAll.append(stringToCSVFeatures);
            //            stringToCSVAll.append(stringToCSVCIVIC);
            //            stringToCSVAll.append(stringToCSVMolecularMatch);
            stringToCSVAll.append(stringToCSVMolecularMatchTrials);
            stringToCSVAll.append(stringToCSVSource);
            writer.append(stringToCSVAll);
            writer.append("\n");
            index++;

        }
        reader.close();
        writer.close();
    }

    public static void extractBRCAFile(@NotNull String brcaJsonPath) throws IOException {
    }

    public static void extractCgiFile(@NotNull String cgiJsonPath) throws IOException {
    }

    public static void extractCivicFile(@NotNull String civicJsonPath) throws IOException {
    }

    public static void extractJaxFile(@NotNull String jaxJsonPath) throws IOException {
    }

    public static void extractJaxTrialsFile(@NotNull String jaxTrialsJsonPath) throws IOException {
    }

    public static void extractMolecularMatchFile(@NotNull String molecularMatchJsonPath) throws IOException {
    }

    public static void extractMolecularMatchTrailsFile(@NotNull String molecularMatchTrialsJsonPath) throws IOException {
    }

    public static void extractOncokbFile(@NotNull String oncokbJsonPath) throws IOException {
    }

    public static void extractPmkbFile(@NotNull String pmkbJsonPath) throws IOException {
    }

    public static void extractSageFile(@NotNull String sageJsonPath) throws IOException {
    }

}
