package com.hartwig.hmftools.vicc.dao;

import static com.hartwig.hmftools.vicc.database.Tables.APPROVEDCOUNTRY;
import static com.hartwig.hmftools.vicc.database.Tables.ASSOCIATION;
import static com.hartwig.hmftools.vicc.database.Tables.DEVTAG;
import static com.hartwig.hmftools.vicc.database.Tables.ENVIRONMENTALCONTEXT;
import static com.hartwig.hmftools.vicc.database.Tables.EVIDENCE;
import static com.hartwig.hmftools.vicc.database.Tables.EVIDENCEINFO;
import static com.hartwig.hmftools.vicc.database.Tables.EVIDENCETYPE;
import static com.hartwig.hmftools.vicc.database.Tables.FEATURE;
import static com.hartwig.hmftools.vicc.database.Tables.FEATURENAME;
import static com.hartwig.hmftools.vicc.database.Tables.GENE;
import static com.hartwig.hmftools.vicc.database.Tables.GENEIDENTIFIER;
import static com.hartwig.hmftools.vicc.database.Tables.HIERARCHY;
import static com.hartwig.hmftools.vicc.database.Tables.LINK;
import static com.hartwig.hmftools.vicc.database.Tables.PHENOTYPE;
import static com.hartwig.hmftools.vicc.database.Tables.PHENOTYPETYPE;
import static com.hartwig.hmftools.vicc.database.Tables.PROVENANCE;
import static com.hartwig.hmftools.vicc.database.Tables.PUBLICATIONURL;
import static com.hartwig.hmftools.vicc.database.Tables.SAGE;
import static com.hartwig.hmftools.vicc.database.Tables.SEQUENCEONTOLOGY;
import static com.hartwig.hmftools.vicc.database.Tables.SYNONYM;
import static com.hartwig.hmftools.vicc.database.Tables.TAG;
import static com.hartwig.hmftools.vicc.database.Tables.TAXONOMY;
import static com.hartwig.hmftools.vicc.database.Tables.VICCENTRY;
import static com.hartwig.hmftools.vicc.database.Tables.BRCA;

import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.SQLException;
import java.util.List;

import com.hartwig.hmftools.vicc.ViccJsonToSQLImporter;
import com.hartwig.hmftools.vicc.database.Tables;
import com.hartwig.hmftools.vicc.datamodel.Association;
import com.hartwig.hmftools.vicc.datamodel.BRCA;
import com.hartwig.hmftools.vicc.datamodel.EnvironmentalContext;
import com.hartwig.hmftools.vicc.datamodel.Evidence;
import com.hartwig.hmftools.vicc.datamodel.EvidenceInfo;
import com.hartwig.hmftools.vicc.datamodel.EvidenceType;
import com.hartwig.hmftools.vicc.datamodel.Feature;
import com.hartwig.hmftools.vicc.datamodel.GeneIdentifier;
import com.hartwig.hmftools.vicc.datamodel.KbSpecificObject;
import com.hartwig.hmftools.vicc.datamodel.Phenotype;
import com.hartwig.hmftools.vicc.datamodel.PhenotypeType;
import com.hartwig.hmftools.vicc.datamodel.Sage;
import com.hartwig.hmftools.vicc.datamodel.SequenceOntology;
import com.hartwig.hmftools.vicc.datamodel.Taxonomy;
import com.hartwig.hmftools.vicc.datamodel.ViccEntry;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;
import org.jooq.DSLContext;
import org.jooq.SQLDialect;
import org.jooq.impl.DSL;

public class ViccDAO {

    private static final Logger LOGGER = LogManager.getLogger(ViccJsonToSQLImporter.class);

    @NotNull
    private final DSLContext context;

    public static ViccDAO connectToViccDAO(@NotNull final String userName, @NotNull final String password, @NotNull final String url)
            throws SQLException {
        final Connection conn = DriverManager.getConnection(url, userName, password);
        final String catalog = conn.getCatalog();

        LOGGER.debug("Connecting to database {}", catalog);
        return new ViccDAO(DSL.using(conn, SQLDialect.MYSQL));
    }

    private ViccDAO(@NotNull final DSLContext context) {
        this.context = context;
    }

    public void writeViccEntry(@NotNull ViccEntry viccEntry) {
        int id = context.insertInto(VICCENTRY, VICCENTRY.SOURCE)
                .values(viccEntry.source())
                .returning(VICCENTRY.ID)
                .fetchOne()
                .getValue(VICCENTRY.ID);
        LOGGER.info("vicc entry id: " + id);
        writeTags(id, viccEntry.tags());
        writeDevTags(id, viccEntry.devTags());
        writeGeneIdentifiers(id, viccEntry.geneIdentifiers());
        writeGenes(id, viccEntry.genes());
        writeFeatureNames(id, viccEntry.featureNames());
        writeFeatures(id, viccEntry.features());
        writeAssociation(id, viccEntry.association());
        writeKbSpecificObject(id, viccEntry.KbSpecificObject());
    }

    private void writeKbSpecificObject(int viccEntryId, @NotNull KbSpecificObject object) {
        if (object instanceof Sage) {
            Sage sage = (Sage) object;
            context.insertInto(SAGE,
                    SAGE.ENTREZID,
                    SAGE.CLINICALMANIFESTATION,
                    SAGE.PUBLICATIONURL,
                    SAGE.GERMLINEORSOMATIC,
                    SAGE.EVIDENCELABEL,
                    SAGE.DRUGLABEL,
                    SAGE.RESPONSETYPE,
                    SAGE.GENE,
                    SAGE.VICCENTRYID)
                    .values(sage.entrezId(),
                            sage.clinicalManifestation(),
                            sage.publicationUrl(),
                            sage.germlineOrSomatic(),
                            sage.evidenceLabel(),
                            sage.drugLabels(),
                            sage.responseType(),
                            sage.gene(),
                            viccEntryId)
                    .execute();
        }
        if (object instanceof BRCA) {
            BRCA brca = (BRCA) object;
            context.insertInto(BRCA,
                    BRCA.VARIANT_FREQUENCY_LOVD,
                    BRCA.ALLELE_FREQUENCY_FIN_EXAC,
                    BRCA.CLINVARACCESSION_ENIGMA,
                    BRCA.HOMOZYGOUS_COUNT_AFR_EXAC,
                    BRCA.BX_ID_EXAC,
                    BRCA.VARIANT_IN_LOVD,
                    BRCA.ALLELE_FREQUENCY_AFR_EXAC,
                    BRCA.DBID_LOVD,
                    BRCA.CHR,
                    BRCA.BX_ID_ENIGMA,
                    BRCA.CO_OCCURRENCE_LR_EXLOVD,
                    BRCA.HOMOZYGOUS_COUNT_EAS_EXAC,
                    BRCA.SUBMITTER_CLINVAR,
                    BRCA.ALLELE_FREQUENCY_EAS_EXAC,
                    BRCA.HG37_END,
                    BRCA.SUBMITTERS_LOVD,
                    BRCA.CLINICAL_CLASSIFICATION_BIC,
                    BRCA.HOMOZYGOUS_COUNT_NFE_EXAC,
                    BRCA.ALLELE_COUNT_SAS_EXAC,
                    BRCA.METHOD_CLINVAR,
                    BRCA.ALLELE_COUNT_NFE_EXAC,
                    BRCA.PATHOGENICITY_ALL,
                    BRCA.GERMLINE_OR_SOMATIC_BIC,
                    BRCA.HOMOZYGOUS_COUNT_SAS_EXAC,
                    BRCA.BIC_NOMENCLATURE,
                    BRCA.ASSERTION_METHOD_ENIGMA,
                    BRCA.LITERATURE_SOURCE_EXLOVD,
                    BRCA.CHANGE_TYPE_ID,
                    BRCA.COLLECTION_METHOD_ENIGMA,
                    BRCA.SUM_FAMILY_LR_EXLOVD,
                    BRCA.HGVS_CDNA_LOVD,
                    BRCA.HOMOZYGOUS_COUNT_FIN_EXAC,
                    BRCA.EAS_ALLELE_FREQUENCY_1000_GENOMES,
                    BRCA.ETHNICITY_BIC,
                    BRCA.INDIVIDUALS_LOVD,
                    BRCA.VARIANT_IN_EXAC,
                    BRCA.URL_ENIGMA,
                    BRCA.ALLELE_ORIGIN_CLINVAR,
                    BRCA.ALLELE_FREQUENCY_AMR_EXAC,
                    BRCA.VARIANT_IN_1000_GENOMES,
                    BRCA.AFR_ALLELE_FREQUENCY_1000_GENOMES,
                    BRCA.BX_ID_EXLOVD,
                    BRCA.SOURCE,
                    BRCA.CONDITION_ID_VALUE_ENIGMA,
                    BRCA.HGVS_PROTEIN,
                    BRCA.REF,
                    BRCA.ALLELE_NUMBER_AFR_EXAC,
                    BRCA.ALLELE_COUNT_AFR_EXAC,
                    BRCA.BX_ID_LOVD,
                    BRCA.SYNONYMS,
                    BRCA.GENE_SYMBOL,
                    BRCA.COMMENT_ON_CLINICAL_SIGNIFICANCE_ENIGMA,
                    BRCA.MISSENSE_ANALYSIS_PRIOR_PROBABILITY_EXLOVD,
                    BRCA.ALLELE_NUMBER_FIN_EXAC,
                    BRCA.POSTERIOR_PROBABILITY_EXLOVD,
                    BRCA.POLYPHEN_SCORE,
                    BRCA.REFERENCE_SEQUENCE,
                    BRCA.ALLELE_COUNT_EAS_EXAC,
                    BRCA.HG38_END,
                    BRCA.HGVS_CDNA,
                    BRCA.FUNCTIONAL_ANALYSIS_TECHNIQUE_LOVD,
                    BRCA.SAS_ALLELE_FREQUENCY_1000_GENOMES,
                    BRCA.RNA_LOVD,
                    BRCA.COMBINED_PRIOR_PROBABLILITY_EXLOVD,
                    BRCA.BX_ID_CLINVAR,
                    BRCA.IARC_CLASS_EXLOVD,
                    BRCA.BX_ID_BIC,
                    BRCA.SIFT_PREDICTION,
                    BRCA.ALLELE_NUMBER_NFE_EXAC,
                    BRCA.ALLELE_ORIGIN_ENIGMA,
                    BRCA.ALLELE_NUMBER_OTH_EXAC,
                    BRCA.HG36_END,
                    BRCA.ALLELE_FREQUENCY_SAS_EXAC,
                    BRCA.DATE_LAST_UPDATED_CLINVAR,
                    BRCA.ALLELE_NUMBER_EAS_EXAC,
                    BRCA.ALLELE_FREQUENCY_OTH_EXAC,
                    BRCA.SOURCE_URL,
                    BRCA.SCV_CLINVAR,
                    BRCA.PATHOGENICITY_EXPERT,
                    BRCA.ALLELE_FREQUENCY_1000_GENOMES,
                    BRCA.FUNCTIONAL_ANALYSIS_RESULT_LOVD,
                    BRCA.AMR_ALLELE_FREQUENCY_1000_GENOMES,
                    BRCA.VARIANT_IN_ESP,
                    BRCA.VARIANT_IN_BIC,
                    BRCA.CLINICAL_SIGNIFICANCE_ENIGMA,
                    BRCA.MAX_ALLELE_FREQUENCY,
                    BRCA.ALLELE_COUNT_AMR_EXAC,
                    BRCA.VARIANT_IN_ENIGMA,
                    BRCA.BX_ID_ESP,
                    BRCA.PATIENT_NATIONALITY_BIC,
                    BRCA.BX_ID_1000_GENOMES,
                    BRCA.GENOMIC_COORDINATE_HG37,
                    BRCA.GENOMIC_COORDINATE_HG36,
                    BRCA.EUR_ALLELE_FREQUENCY_1000_GENOMES,
                    BRCA.NUMBER_OF_FAMILY_MEMBER_CARRYING_MUTATION_BIC,
                    BRCA.SEGREGATION_LR_EXLOVD,
                    BRCA.ALLELE_FREQUENCY,
                    BRCA.MINOR_ALLELE_FREQUENCY_PERCENT_ESP,
                    BRCA.ALLELE_FREQUENCY_EXAC,
                    BRCA.MUTATION_TYPE_BIC,
                    BRCA.ASSERTION_METHOD_CITATION_ENIGMA,
                    BRCA.CONDITION_ID_TYPE_ENIGMA,
                    BRCA.ALLELE_COUNT_OTH_EXAC,
                    BRCA.HGVS_PROTEIN_LOVD,
                    BRCA.VARIANT_IN_CLINVAR,
                    BRCA.CLINICAL_IMPORTANCE_BIC,
                    BRCA.DISCORDANT,
                    BRCA.ALLELE_COUNT_FIN_EXAC,
                    BRCA.CONDITION_CATEGORY_ENIGMA,
                    BRCA.ALLELE_FREQUENCY_ESP,
                    BRCA.HOMOZYGOUS_COUNT_OTH_EXAC,
                    BRCA.GENETIC_ORIGIN_LOVD,
                    BRCA.HOMOZYGOUS_COUNT_AMR_EXAC,
                    BRCA.CLINICAL_SIGNIFICANCE_CLINVAR,
                    BRCA.AA_ALLELE_FREQUENCY_ESP,
                    BRCA.PROTEIN_CHANGE,
                    BRCA.VARIANT_IN_EXLOVD,
                    BRCA.EA_ALLELE_FREQUENCY_ESP,
                    BRCA.HGVS_RNA,
                    BRCA.CLINICAL_SIGNIFICANCE_CITATIONS_ENIGMA,
                    BRCA.VARIANT_EFFECT_LOVD,
                    BRCA.POLYPHEN_PREDICTION,
                    BRCA.DATA_RELEASE_ID,
                    BRCA.HG37_START,
                    BRCA.HG36_START,
                    BRCA.SIFT_SCORE,
                    BRCA.GENOMIC_COORDINATE_HG38,
                    BRCA.ALT,
                    BRCA.LITERATURE_CITATION_BIC,
                    BRCA.VARIANT_HAPLOTYPE_LOVD,
                    BRCA.ALLELE_FREQUENCY_NFE_EXAC,
                    BRCA.HG38_START,
                    BRCA.POS,
                    BRCA.DATE_LAST_EVALUATED_ENIGMA,
                    BRCA.ALLELE_NUMBER_SAS_EXAC,
                    BRCA.ALLELE_NUMBER_AMR_EXAC,
                    BRCA.VICCENTRYID)
                    .values(brca.brcApart1().Variant_frequency_LOVD(),
                            brca.brcApart1().Allele_frequency_FIN_ExAC(),
                            brca.brcApart1().ClinVarAccession_ENIGMA(),
                            brca.brcApart1().Homozygous_count_AFR_ExAC(),
                            brca.brcApart1().BX_ID_ExAC(),
                            brca.brcApart1().Variant_in_LOVD(),
                            brca.brcApart1().Allele_frequency_AFR_ExAC(),
                            brca.brcApart2().DBID_LOVD(),
                            brca.brcApart1().Chr(),
                            brca.brcApart1().BX_ID_ENIGMA(),
                            brca.brcApart1().Co_occurrence_LR_exLOVD(),
                            brca.brcApart1().Homozygous_count_EAS_ExAC(),
                            brca.brcApart1().Submitter_ClinVar(),
                            brca.brcApart1().Allele_frequency_EAS_ExAC(),
                            brca.brcApart1().Hg37_End(),
                            brca.brcApart1().Submitters_LOVD(),
                            brca.brcApart1().Clinical_classification_BIC(),
                            brca.brcApart1().Homozygous_count_NFE_ExAC(),
                            brca.brcApart1().Allele_count_SAS_ExAC(),
                            brca.brcApart1().Method_ClinVar(),
                            brca.brcApart1().Allele_count_NFE_ExAC(),
                            brca.brcApart1().Pathogenicity_all(),
                            brca.brcApart1().Germline_or_Somatic_BIC(),
                            brca.brcApart1().Homozygous_count_SAS_ExAC(),
                            brca.brcApart1().BIC_Nomenclature(),
                            brca.brcApart1().Assertion_method_ENIGMA(),
                            brca.brcApart1().Literature_source_exLOVD(),
                            brca.brcApart1().Change_Type_id(),
                            brca.brcApart1().Collection_method_ENIGMA(),
                            brca.brcApart1().Sum_family_LR_exLOVD(),
                            brca.brcApart1().HGVS_cDNA_LOVD(),
                            brca.brcApart1().Homozygous_count_FIN_ExAC(),
                            brca.brcApart1().EAS_Allele_frequency_1000_Genomes(),
                            brca.brcApart1().Ethnicity_BIC(),
                            brca.brcApart1().Individuals_LOVD(),
                            brca.brcApart1().Variant_in_ExAC(),
                            brca.brcApart1().URL_ENIGMA(),
                            brca.brcApart1().Allele_Origin_ClinVar(),
                            brca.brcApart1().Allele_frequency_AMR_ExAC(),
                            brca.brcApart1().Variant_in_1000_Genomes(),
                            brca.brcApart1().AFR_Allele_frequency_1000_Genomes(),
                            brca.brcApart1().BX_ID_exLOVD(),
                            brca.brcApart1().Source(),
                            brca.brcApart1().Condition_ID_value_ENIGMA(),
                            brca.brcApart1().HGVS_Protein(),
                            brca.brcApart1().Ref(),
                            brca.brcApart1().Allele_number_AFR_ExAC(),
                            brca.brcApart1().Allele_count_AFR_ExAC(),
                            brca.brcApart1().BX_ID_LOVD(),
                            brca.brcApart1().Synonyms(),
                            brca.brcApart1().Gene_Symbol(),
                            brca.brcApart1().Comment_on_clinical_significance_ENIGMA(),
                            brca.brcApart1().Missense_analysis_prior_probability_exLOVD(),
                            brca.brcApart1().Allele_number_FIN_ExAC(),
                            brca.brcApart1().Posterior_probability_exLOVD(),
                            brca.brcApart1().Polyphen_Score(),
                            brca.brcApart1().Reference_Sequence(),
                            brca.brcApart1().Allele_count_EAS_ExAC(),
                            brca.brcApart1().Hg38_End(),
                            brca.brcApart1().HGVS_cDNA(),
                            brca.brcApart1().Functional_analysis_technique_LOVD(),
                            brca.brcApart1().SAS_Allele_frequency_1000_Genomes(),
                            brca.brcApart1().RNA_LOVD(),
                            brca.brcApart1().Combined_prior_probablility_exLOVD(),
                            brca.brcApart1().BX_ID_ClinVar(),
                            brca.brcApart1().IARC_class_exLOVD(),
                            brca.brcApart1().BX_ID_BIC(),
                            brca.brcApart1().Sift_Prediction(),
                            brca.brcApart1().Allele_number_NFE_ExAC(),
                            brca.brcApart1().Allele_origin_ENIGMA(),
                            brca.brcApart1().Allele_number_OTH_ExAC(),
                            brca.brcApart1().Hg36_End(),
                            brca.brcApart1().Allele_frequency_SAS_ExAC(),
                            brca.brcApart1().Date_Last_Updated_ClinVar(),
                            brca.brcApart1().Allele_number_EAS_ExAC(),
                            brca.brcApart1().Allele_frequency_OTH_ExAC(),
                            brca.brcApart1().Source_URL(),
                            brca.brcApart1().SCV_ClinVar(),
                            brca.brcApart1().Pathogenicity_expert(),
                            brca.brcApart1().Allele_frequency_1000_Genomes(),
                            brca.brcApart1().Functional_analysis_result_LOVD(),
                            brca.brcApart1().AMR_Allele_frequency_1000_Genomes(),
                            brca.brcApart1().Variant_in_ESP(),
                            brca.brcApart1().Variant_in_BIC(),
                            brca.brcApart1().Clinical_significance_ENIGMA(),
                            brca.brcApart1().Max_Allele_Frequency(),
                            brca.brcApart1().Allele_count_AMR_ExAC(),
                            brca.brcApart1().Variant_in_ENIGMA(),
                            brca.brcApart1().BX_ID_ESP(),
                            brca.brcApart1().Patient_nationality_BIC(),
                            brca.brcApart1().BX_ID_1000_Genomes(),
                            brca.brcApart1().Genomic_Coordinate_hg37(),
                            brca.brcApart1().Genomic_Coordinate_hg36(),
                            brca.brcApart1().EUR_Allele_frequency_1000_Genomes(),
                            brca.brcApart1().Number_of_family_member_carrying_mutation_BIC(),
                            brca.brcApart1().Segregation_LR_exLOVD(),
                            brca.brcApart1().Allele_Frequency(),
                            brca.brcApart1().Minor_allele_frequency_percent_ESP(),
                            brca.brcApart1().Allele_frequency_ExAC(),
                            brca.brcApart1().Mutation_type_BIC(),
                            brca.brcApart1().Assertion_method_citation_ENIGMA(),
                            brca.brcApart1().Condition_ID_type_ENIGMA(),
                            brca.brcApart1().Allele_count_OTH_ExAC(),
                            brca.brcApart1().HGVS_protein_LOVD(),
                            brca.brcApart1().Variant_in_ClinVar(),
                            brca.brcApart1().Clinical_importance_BIC(),
                            brca.brcApart1().Discordant(),
                            brca.brcApart2().Allele_count_FIN_ExAC(),
                            brca.brcApart2().Condition_category_ENIGMA(),
                            brca.brcApart2().Allele_Frequency_ESP(),
                            brca.brcApart2().Homozygous_count_OTH_ExAC(),
                            brca.brcApart2().Genetic_origin_LOVD(),
                            brca.brcApart2().Homozygous_count_AMR_ExAC(),
                            brca.brcApart2().Clinical_Significance_ClinVar(),
                            brca.brcApart2().AA_Allele_Frequency_ESP(),
                            brca.brcApart2().Protein_Change(),
                            brca.brcApart2().Variant_in_exLOVD(),
                            brca.brcApart2().EA_Allele_Frequency_ESP(),
                            brca.brcApart2().HGVS_RNA(),
                            brca.brcApart2().Clinical_significance_citations_ENIGMA(),
                            brca.brcApart2().Variant_effect_LOVD(),
                            brca.brcApart2().Polyphen_Prediction(),
                            brca.brcApart2().Data_Release_id(),
                            brca.brcApart2().Hg37_Start(),
                            brca.brcApart2().Hg36_Start(),
                            brca.brcApart2().Sift_Score(),
                            brca.brcApart2().Genomic_Coordinate_hg38(),
                            brca.brcApart2().Alt(),
                            brca.brcApart2().Literature_citation_BIC(),
                            brca.brcApart2().Variant_haplotype_LOVD(),
                            brca.brcApart2().Allele_frequency_NFE_ExAC(),
                            brca.brcApart2().Hg38_Start(),
                            brca.brcApart2().Pos(),
                            brca.brcApart2().Date_last_evaluated_ENIGMA(),
                            brca.brcApart2().Allele_number_SAS_ExAC(),
                            brca.brcApart2().Allele_number_AMR_ExAC(),
                            viccEntryId)
                    .execute();
        }
    }

    private void writeAssociation(int viccEntryId, @NotNull Association association) {
        int id = context.insertInto(ASSOCIATION,
                ASSOCIATION.VARIANTNAME,
                ASSOCIATION.EVIDENCELEVEL,
                ASSOCIATION.EVIDENCELABEL,
                ASSOCIATION.RESPONSETYPE,
                ASSOCIATION.DRUGLABEL,
                ASSOCIATION.SOURCELINK,
                ASSOCIATION.DESCRIPTION,
                ASSOCIATION.ONCOGENIC,
                ASSOCIATION.VICCENTRYID)
                .values(association.variantName(),
                        association.evidenceLevel(),
                        association.evidenceLabel(),
                        association.responseType(),
                        association.drugLabels(),
                        association.sourceLink(),
                        association.description(),
                        association.oncogenic(),
                        viccEntryId)
                .returning(ASSOCIATION.ID)
                .fetchOne()
                .getValue(ASSOCIATION.ID);
        writeEvidence(id, association.evidence());
        writePublicationsUrls(id, association.publicationUrls());
        writePhenotype(id, association.phenotype());
        writeEnvironmentalContexts(id, association.environmentalContexts());
    }

    private void writeEnvironmentalContexts(int associationId, @Nullable List<EnvironmentalContext> environmentalContexts) {
        if (environmentalContexts != null) {
            for (EnvironmentalContext environmentalContext : environmentalContexts) {
                int id = context.insertInto(ENVIRONMENTALCONTEXT,
                        ENVIRONMENTALCONTEXT.TERM,
                        ENVIRONMENTALCONTEXT.DESCRIPTION,
                        ENVIRONMENTALCONTEXT.SOURCE,
                        ENVIRONMENTALCONTEXT.USANSTEM,
                        ENVIRONMENTALCONTEXT.IDENVIRONMENTALCONTEXT,
                        ENVIRONMENTALCONTEXT.ASSOCIATIONID)
                        .values(environmentalContext.term(),
                                environmentalContext.description(),
                                environmentalContext.source(),
                                environmentalContext.usanStem(),
                                environmentalContext.id(),
                                associationId)
                        .returning(ENVIRONMENTALCONTEXT.ID)
                        .fetchOne()
                        .getValue(ENVIRONMENTALCONTEXT.ID);
                writeApprovedCountries(id, environmentalContext.approvedCountries());
                writeTaxonomy(id, environmentalContext.taxonomy());
            }
        }
    }

    private void writeTaxonomy(int environmentalContextsId, @Nullable Taxonomy taxonomy) {
        if (taxonomy != null) {
            context.insertInto(TAXONOMY,
                    TAXONOMY.KINGDOM,
                    TAXONOMY.DIRECTPARENT,
                    TAXONOMY.CLASS,
                    TAXONOMY.SUBCLASS,
                    TAXONOMY.SUPERCLASS,
                    TAXONOMY.ENVIRONMENTALCONTEXTSID)
                    .values(taxonomy.kingdom(),
                            taxonomy.directParent(),
                            taxonomy.classs(),
                            taxonomy.subClass(),
                            taxonomy.superClass(),
                            environmentalContextsId)
                    .execute();
        }
    }

    private void writeApprovedCountries(int environmentalContextsId, @NotNull List<String> approvedCountries) {
        for (String approvesCountry : approvedCountries) {
            context.insertInto(APPROVEDCOUNTRY, APPROVEDCOUNTRY.APPROVEDCOUNTRYNAME, APPROVEDCOUNTRY.ENVIRONMENTALCONTEXTSID)
                    .values(approvesCountry, environmentalContextsId)
                    .execute();
        }
    }

    private void writePhenotype(int associationId, @Nullable Phenotype phenotype) {
        if (phenotype != null) {
            int id = context.insertInto(PHENOTYPE, PHENOTYPE.DESCRIPTION, PHENOTYPE.FAMILY, PHENOTYPE.IDPHENOTYPE, PHENOTYPE.ASSOCIATIONID)
                    .values(phenotype.description(), phenotype.family(), phenotype.id(), associationId)
                    .returning(PHENOTYPE.ID)
                    .fetchOne()
                    .getValue(PHENOTYPE.ID);
            writePhenotypeType(id, phenotype.type());
        }
    }

    private void writePhenotypeType(int phenotypeId, @Nullable PhenotypeType phenotypeType) {
        if (phenotypeType != null) {
            context.insertInto(PHENOTYPETYPE,
                    PHENOTYPETYPE.SOURCE,
                    PHENOTYPETYPE.TERM,
                    PHENOTYPETYPE.IDPHENOTYPETYPE,
                    PHENOTYPETYPE.PHENOTYPEID)
                    .values(phenotypeType.source(), phenotypeType.term(), phenotypeType.id(), phenotypeId)
                    .execute();
        }
    }

    private void writePublicationsUrls(int associationId, @Nullable List<String> publicationsUrls) {
        if (publicationsUrls != null) {
            for (String publicationUrl : publicationsUrls) {
                context.insertInto(PUBLICATIONURL, PUBLICATIONURL.URLOFPUBLICATION, PUBLICATIONURL.ASSOCIATIONID)
                        .values(publicationUrl, associationId)
                        .execute();
            }
        }
    }

    private void writeEvidence(int associationId, @NotNull List<Evidence> evidences) {
        for (Evidence evidence : evidences) {
            int id = context.insertInto(EVIDENCE, EVIDENCE.DESCRIPTION, EVIDENCE.ASSOCIATIONID)
                    .values(evidence.description(), associationId)
                    .returning(EVIDENCE.ID)
                    .fetchOne()
                    .getValue(EVIDENCE.ID);
            writeEvidenceInfo(id, evidence.info());
            writeEvidenceType(id, evidence.evidenceType());
        }
    }

    private void writeEvidenceInfo(int evidenceId, @Nullable EvidenceInfo evidenceInfo) {
        if (evidenceInfo != null) {
            for (String publication : evidenceInfo.publications()) {
                context.insertInto(EVIDENCEINFO, EVIDENCEINFO.PUBLICATION, EVIDENCEINFO.EVIDENCEID)
                        .values(publication, evidenceId)
                        .execute();
            }
        }
    }

    private void writeEvidenceType(int evidenceId, @NotNull EvidenceType evidenceType) {
        context.insertInto(EVIDENCETYPE, EVIDENCETYPE.SOURCENAME, EVIDENCETYPE.IDEVIDENCETYPE, EVIDENCETYPE.EVIDENCEID)
                .values(evidenceType.sourceName(), evidenceType.id(), evidenceId)
                .execute();
    }

    private void writeFeatures(int viccEntryId, @NotNull List<Feature> features) {
        for (Feature feature : features) {
            int id = context.insertInto(FEATURE,
                    FEATURE.NAME,
                    FEATURE.BIOMARKERTYPE,
                    FEATURE.REFERENCENAME,
                    FEATURE.CHROMOSOME,
                    FEATURE.START,
                    FEATURE.END,
                    FEATURE.REF,
                    FEATURE.ALT,
                    FEATURE.PROVENANCERULE,
                    FEATURE.GENESYMBOL,
                    FEATURE.ENTREZID,
                    FEATURE.DESCRIPTION,
                    FEATURE.VICCENTRYID)
                    .values(feature.name(),
                            feature.biomarkerType(),
                            feature.referenceName(),
                            feature.chromosome(),
                            feature.start(),
                            feature.end(),
                            feature.ref(),
                            feature.alt(),
                            feature.provenanceRule(),
                            feature.geneSymbol(),
                            feature.entrezId(),
                            feature.description(),
                            viccEntryId)
                    .returning(FEATURE.ID)
                    .fetchOne()
                    .getValue(FEATURE.ID);
            writeProvenance(id, feature.provenance());
            writeSynonyms(id, feature.synonyms());
            writeLinks(id, feature.links());
            writeSequenceOntology(id, feature.sequenceOntology());
        }
    }

    private void writeProvenance(int featureId, @NotNull List<String> provenances) {
        for (String provenance : provenances) {
            context.insertInto(PROVENANCE, PROVENANCE.PROVENANCENAME, PROVENANCE.FEATUREID).values(provenance, featureId).execute();
        }
    }

    private void writeSynonyms(int featureId, @Nullable List<String> synonyms) {
        if (synonyms != null) {
            for (String synonym : synonyms) {
                context.insertInto(SYNONYM, SYNONYM.SYNONYMNAME, SYNONYM.FEATUREID).values(synonym, featureId).execute();
            }
        }
    }

    private void writeLinks(int featureId, @Nullable List<String> links) {
        if (links != null) {
            for (String link : links) {
                context.insertInto(LINK, LINK.LINKNAME, LINK.FEATUREID).values(link, featureId).execute();
            }
        }
    }

    private void writeSequenceOntology(int featureId, @Nullable SequenceOntology sequenceOntologies) {
        if (sequenceOntologies != null) {
            int id = context.insertInto(SEQUENCEONTOLOGY,
                    SEQUENCEONTOLOGY.SOID,
                    SEQUENCEONTOLOGY.PARENTSOID,
                    SEQUENCEONTOLOGY.NAME,
                    SEQUENCEONTOLOGY.PARENTNAME,
                    SEQUENCEONTOLOGY.FEATUREID)
                    .values(sequenceOntologies.soid(),
                            sequenceOntologies.parentSoid(),
                            sequenceOntologies.name(),
                            sequenceOntologies.parentName(),
                            featureId)
                    .returning(SEQUENCEONTOLOGY.ID)
                    .fetchOne()
                    .getValue(SEQUENCEONTOLOGY.ID);
            writeHierarchy(id, sequenceOntologies.hierarchy());
        }
    }

    private void writeHierarchy(int sequenceOntologyId, @Nullable List<String> hierarchies) {
        if (hierarchies != null) {
            for (String hierarchy : hierarchies) {
                context.insertInto(HIERARCHY, HIERARCHY.HIERARCHYNAME, HIERARCHY.SEQUENCEONTOLOGYID)
                        .values(hierarchy, sequenceOntologyId)
                        .execute();
            }
        }
    }

    private void writeTags(int viccEntryId, @NotNull List<String> tags) {
        for (String tag : tags) {
            context.insertInto(TAG, TAG.TAGNAME, TAG.VICCENTRYID).values(tag, viccEntryId).execute();
        }
    }

    private void writeDevTags(int viccEntryId, @NotNull List<String> devTags) {
        for (String devTag : devTags) {
            context.insertInto(DEVTAG, DEVTAG.DEVTAGNAME, DEVTAG.VICCENTRYID).values(devTag, viccEntryId).execute();
        }
    }

    private void writeGeneIdentifiers(int viccEntryId, @NotNull List<GeneIdentifier> geneIdentifiers) {
        for (GeneIdentifier geneIdentifier : geneIdentifiers) {
            context.insertInto(GENEIDENTIFIER,
                    GENEIDENTIFIER.SYMBOL,
                    GENEIDENTIFIER.ENTREZID,
                    GENEIDENTIFIER.ENSEMBLGENEID,
                    GENEIDENTIFIER.VICCENTRYID)
                    .values(geneIdentifier.symbol(), geneIdentifier.entrezId(), geneIdentifier.ensemblGeneId(), viccEntryId)
                    .execute();
        }
    }

    private void writeGenes(int viccEntryId, @NotNull List<String> genes) {
        for (String gene : genes) {
            context.insertInto(GENE, GENE.GENENAME, GENE.VICCENTRYID).values(gene, viccEntryId).execute();
        }
    }

    private void writeFeatureNames(int viccEntryId, @Nullable List<String> featureNames) {
        if (featureNames != null) {
            for (String featureName : featureNames) {
                context.insertInto(FEATURENAME, FEATURENAME.NAMEOFFEATURE, FEATURENAME.VICCENTRYID)
                        .values(featureName, viccEntryId)
                        .execute();
            }
        }
    }

    public void deleteAll() {
        LOGGER.info("Deleting all from vicc db");

        context.deleteFrom(TAG).execute();
        context.deleteFrom(DEVTAG).execute();
        context.deleteFrom(GENEIDENTIFIER).execute();
        context.deleteFrom(GENE).execute();
        context.deleteFrom(FEATURENAME).execute();
        context.deleteFrom(FEATURE).execute();
        context.deleteFrom(PROVENANCE).execute();
        context.deleteFrom(SYNONYM).execute();
        context.deleteFrom(LINK).execute();
        context.deleteFrom(SEQUENCEONTOLOGY).execute();
        context.deleteFrom(HIERARCHY).execute();
        context.deleteFrom(ASSOCIATION).execute();
        context.deleteFrom(EVIDENCE).execute();
        context.deleteFrom(EVIDENCEINFO).execute();
        context.deleteFrom(EVIDENCETYPE).execute();
        context.deleteFrom(PUBLICATIONURL).execute();
        context.deleteFrom(PHENOTYPE).execute();
        context.deleteFrom(PHENOTYPETYPE).execute();
        context.deleteFrom(ENVIRONMENTALCONTEXT).execute();
        context.deleteFrom(APPROVEDCOUNTRY).execute();
        context.deleteFrom(TAXONOMY).execute();
        context.deleteFrom(VICCENTRY).execute();
        context.deleteFrom(SAGE).execute();
    }
}
