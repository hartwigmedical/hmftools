package com.hartwig.hmftools.vicc.dao;

import static com.hartwig.hmftools.vicc.database.Tables.APPROVEDCOUNTRY;
import static com.hartwig.hmftools.vicc.database.Tables.ASSOCIATION;
import static com.hartwig.hmftools.vicc.database.Tables.BRCAPART1;
import static com.hartwig.hmftools.vicc.database.Tables.BRCAPART2;
import static com.hartwig.hmftools.vicc.database.Tables.CGI;
import static com.hartwig.hmftools.vicc.database.Tables.CGICDNA;
import static com.hartwig.hmftools.vicc.database.Tables.CGIGDNA;
import static com.hartwig.hmftools.vicc.database.Tables.CGIINDIVIDUALMUTATION;
import static com.hartwig.hmftools.vicc.database.Tables.CGIINFO;
import static com.hartwig.hmftools.vicc.database.Tables.CGIREGION;
import static com.hartwig.hmftools.vicc.database.Tables.CGISTRAND;
import static com.hartwig.hmftools.vicc.database.Tables.CGITRANSCRIPT;
import static com.hartwig.hmftools.vicc.database.Tables.CIVIC;
import static com.hartwig.hmftools.vicc.database.Tables.CIVICASSERTIONS;
import static com.hartwig.hmftools.vicc.database.Tables.CIVICCLINICALTRIAL;
import static com.hartwig.hmftools.vicc.database.Tables.CIVICCLINVARENTRIES;
import static com.hartwig.hmftools.vicc.database.Tables.CIVICCOORDINATES;
import static com.hartwig.hmftools.vicc.database.Tables.CIVICDESCRIPTION;
import static com.hartwig.hmftools.vicc.database.Tables.CIVICDISEASE;
import static com.hartwig.hmftools.vicc.database.Tables.CIVICDRUGS;
import static com.hartwig.hmftools.vicc.database.Tables.CIVICERROR;
import static com.hartwig.hmftools.vicc.database.Tables.CIVICEVIDENCEITEMS;
import static com.hartwig.hmftools.vicc.database.Tables.CIVICEVIDENCEITEMSCLINICALTRIAL;
import static com.hartwig.hmftools.vicc.database.Tables.CIVICEVIDENCEITEMSPUBLICATION;
import static com.hartwig.hmftools.vicc.database.Tables.CIVICEVIDENCEITEMSSOURCE;
import static com.hartwig.hmftools.vicc.database.Tables.CIVICHGVSEXPRESSIONS;
import static com.hartwig.hmftools.vicc.database.Tables.CIVICLASTCOMMENTEDON;
import static com.hartwig.hmftools.vicc.database.Tables.CIVICLASTCOMMENTEDONAVATARS;
import static com.hartwig.hmftools.vicc.database.Tables.CIVICLASTCOMMENTEDONORGANIZATION;
import static com.hartwig.hmftools.vicc.database.Tables.CIVICLASTCOMMENTEDONPROFILEIMAGE;
import static com.hartwig.hmftools.vicc.database.Tables.CIVICLASTCOMMENTEDONUSER;
import static com.hartwig.hmftools.vicc.database.Tables.CIVICLASTMODIFIED;
import static com.hartwig.hmftools.vicc.database.Tables.CIVICLASTMODIFIEDAVATARS;
import static com.hartwig.hmftools.vicc.database.Tables.CIVICLASTMODIFIEDORGANIZATION;
import static com.hartwig.hmftools.vicc.database.Tables.CIVICLASTMODIFIEDPROFILEIMAGE;
import static com.hartwig.hmftools.vicc.database.Tables.CIVICLASTMODIFIEDUSER;
import static com.hartwig.hmftools.vicc.database.Tables.CIVICLASTREVIEWED;
import static com.hartwig.hmftools.vicc.database.Tables.CIVICLASTREVIEWEDAVATARS;
import static com.hartwig.hmftools.vicc.database.Tables.CIVICLASTREVIEWEDORGANIZATION;
import static com.hartwig.hmftools.vicc.database.Tables.CIVICLASTREVIEWEDPROFILEIMAGE;
import static com.hartwig.hmftools.vicc.database.Tables.CIVICLASTREVIEWEDUSER;
import static com.hartwig.hmftools.vicc.database.Tables.CIVICLIFECYCLEACTIONS;
import static com.hartwig.hmftools.vicc.database.Tables.CIVICPUBLICATION;
import static com.hartwig.hmftools.vicc.database.Tables.CIVICSOURCE;
import static com.hartwig.hmftools.vicc.database.Tables.CIVICVARIANTALIASES;
import static com.hartwig.hmftools.vicc.database.Tables.CIVICVARIANTSGROUPS;
import static com.hartwig.hmftools.vicc.database.Tables.CIVICVARIANTSGROUPSCOORDINATES;
import static com.hartwig.hmftools.vicc.database.Tables.CIVICVARIANTSGROUPSTYPES;
import static com.hartwig.hmftools.vicc.database.Tables.CIVICVARIANTSGROUPSVARIANTS;
import static com.hartwig.hmftools.vicc.database.Tables.CIVICVARIANTTYPES;
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
import static com.hartwig.hmftools.vicc.database.Tables.JAX;
import static com.hartwig.hmftools.vicc.database.Tables.JAXINDICATIONS;
import static com.hartwig.hmftools.vicc.database.Tables.JAXMOLECULARPROFILE;
import static com.hartwig.hmftools.vicc.database.Tables.JAXREFERENCES;
import static com.hartwig.hmftools.vicc.database.Tables.JAXTHERAPY;
import static com.hartwig.hmftools.vicc.database.Tables.JAXTRIALS;
import static com.hartwig.hmftools.vicc.database.Tables.JAXTRIALSINDICATIONS;
import static com.hartwig.hmftools.vicc.database.Tables.JAXTRIALSMOLECULARPROFILE;
import static com.hartwig.hmftools.vicc.database.Tables.JAXTRIALSTHERAPIES;
import static com.hartwig.hmftools.vicc.database.Tables.JAXTRIALSVARIANTREQUIREMENTDETAILS;
import static com.hartwig.hmftools.vicc.database.Tables.LINK;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHTRIALS;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHTRIALSALTERATIONS;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHTRIALSCONTACT;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHTRIALSCOORDINATES;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHTRIALSGEO;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHTRIALSINTERVATIONS;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHTRIALSLOCATION;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHTRIALSLOCATIONS;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHTRIALSOTHERGROUPLABEL;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHTRIALSOTHERNAME;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHTRIALSOVERALLCONTACT;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHTRIALSTAGS;
import static com.hartwig.hmftools.vicc.database.Tables.ONCOKB;
import static com.hartwig.hmftools.vicc.database.Tables.ONCOKBBIOLOGICAL;
import static com.hartwig.hmftools.vicc.database.Tables.ONCOKBCLINICAL;
import static com.hartwig.hmftools.vicc.database.Tables.ONCOKBCONSEQUENCESBIOLOGICAL;
import static com.hartwig.hmftools.vicc.database.Tables.ONCOKBCONSEQUENCESCLINICAL;
import static com.hartwig.hmftools.vicc.database.Tables.ONCOKBDRUGABSTRACTSCLINICAL;
import static com.hartwig.hmftools.vicc.database.Tables.ONCOKBGENEALIASESBIOLOGICAL;
import static com.hartwig.hmftools.vicc.database.Tables.ONCOKBGENEALIASESCLINICAL;
import static com.hartwig.hmftools.vicc.database.Tables.ONCOKBGENEBIOLOGICAL;
import static com.hartwig.hmftools.vicc.database.Tables.ONCOKBGENECLINICAL;
import static com.hartwig.hmftools.vicc.database.Tables.ONCOKBVARIANTBIOLOGICAL;
import static com.hartwig.hmftools.vicc.database.Tables.ONCOKBVARIANTCLINICAL;
import static com.hartwig.hmftools.vicc.database.Tables.PHENOTYPE;
import static com.hartwig.hmftools.vicc.database.Tables.PHENOTYPETYPE;
import static com.hartwig.hmftools.vicc.database.Tables.PMKB;
import static com.hartwig.hmftools.vicc.database.Tables.PMKBGENE;
import static com.hartwig.hmftools.vicc.database.Tables.PMKBTISSUE;
import static com.hartwig.hmftools.vicc.database.Tables.PMKBTUMOR;
import static com.hartwig.hmftools.vicc.database.Tables.PMKBVARIANT;
import static com.hartwig.hmftools.vicc.database.Tables.PROVENANCE;
import static com.hartwig.hmftools.vicc.database.Tables.PUBLICATIONURL;
import static com.hartwig.hmftools.vicc.database.Tables.SAGE;
import static com.hartwig.hmftools.vicc.database.Tables.SEQUENCEONTOLOGY;
import static com.hartwig.hmftools.vicc.database.Tables.SYNONYM;
import static com.hartwig.hmftools.vicc.database.Tables.TAG;
import static com.hartwig.hmftools.vicc.database.Tables.TAXONOMY;
import static com.hartwig.hmftools.vicc.database.Tables.VICCENTRY;

import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.SQLException;
import java.util.List;
import java.util.Optional;

import com.hartwig.hmftools.vicc.ViccJsonToSQLImporter;
import com.hartwig.hmftools.vicc.database.Tables;
import com.hartwig.hmftools.vicc.database.tables.Oncokbbiological;
import com.hartwig.hmftools.vicc.datamodel.Association;
import com.hartwig.hmftools.vicc.datamodel.BRCA;
import com.hartwig.hmftools.vicc.datamodel.Cgi;
import com.hartwig.hmftools.vicc.datamodel.Civic;
import com.hartwig.hmftools.vicc.datamodel.CivicClinicalTrial;
import com.hartwig.hmftools.vicc.datamodel.CivicDescription;
import com.hartwig.hmftools.vicc.datamodel.CivicDrugs;
import com.hartwig.hmftools.vicc.datamodel.CivicEvidenceItems;
import com.hartwig.hmftools.vicc.datamodel.CivicPublicationDate;
import com.hartwig.hmftools.vicc.datamodel.CivicSource;
import com.hartwig.hmftools.vicc.datamodel.CivicVariantGroup;
import com.hartwig.hmftools.vicc.datamodel.CivicVariantTypes;
import com.hartwig.hmftools.vicc.datamodel.CivicVariants;
import com.hartwig.hmftools.vicc.datamodel.EnvironmentalContext;
import com.hartwig.hmftools.vicc.datamodel.Evidence;
import com.hartwig.hmftools.vicc.datamodel.EvidenceInfo;
import com.hartwig.hmftools.vicc.datamodel.EvidenceType;
import com.hartwig.hmftools.vicc.datamodel.Feature;
import com.hartwig.hmftools.vicc.datamodel.GeneIdentifier;
import com.hartwig.hmftools.vicc.datamodel.Jax;
import com.hartwig.hmftools.vicc.datamodel.JaxReferences;
import com.hartwig.hmftools.vicc.datamodel.JaxTrials;
import com.hartwig.hmftools.vicc.datamodel.JaxTrialsIndications;
import com.hartwig.hmftools.vicc.datamodel.JaxTrialsMolecularProfile;
import com.hartwig.hmftools.vicc.datamodel.JaxTrialsTherapies;
import com.hartwig.hmftools.vicc.datamodel.JaxTrialsVariantRequirementDetails;
import com.hartwig.hmftools.vicc.datamodel.KbSpecificObject;
import com.hartwig.hmftools.vicc.datamodel.MolecularMatchTrials;
import com.hartwig.hmftools.vicc.datamodel.MolecularMatchTrialsContact;
import com.hartwig.hmftools.vicc.datamodel.MolecularMatchTrialsIntervation;
import com.hartwig.hmftools.vicc.datamodel.MolecularMatchTrialsLocations;
import com.hartwig.hmftools.vicc.datamodel.MolecularMatchTrialsTags;
import com.hartwig.hmftools.vicc.datamodel.OncoKbBiological;
import com.hartwig.hmftools.vicc.datamodel.OncoKbClinical;
import com.hartwig.hmftools.vicc.datamodel.OncoKbConsequence;
import com.hartwig.hmftools.vicc.datamodel.OncoKbDrugAbstracts;
import com.hartwig.hmftools.vicc.datamodel.Oncokb;
import com.hartwig.hmftools.vicc.datamodel.OncokbGene;
import com.hartwig.hmftools.vicc.datamodel.OncokbVariant;
import com.hartwig.hmftools.vicc.datamodel.Phenotype;
import com.hartwig.hmftools.vicc.datamodel.PhenotypeType;
import com.hartwig.hmftools.vicc.datamodel.Pmkb;
import com.hartwig.hmftools.vicc.datamodel.PmkbGene;
import com.hartwig.hmftools.vicc.datamodel.PmkbTissue;
import com.hartwig.hmftools.vicc.datamodel.PmkbTumor;
import com.hartwig.hmftools.vicc.datamodel.PmkbVariant;
import com.hartwig.hmftools.vicc.datamodel.Sage;
import com.hartwig.hmftools.vicc.datamodel.SequenceOntology;
import com.hartwig.hmftools.vicc.datamodel.Taxonomy;
import com.hartwig.hmftools.vicc.datamodel.ViccEntry;
import com.sun.org.apache.bcel.internal.generic.IF_ACMPEQ;

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

    private void importSageinSQL(int viccEntryId, @NotNull KbSpecificObject object) {
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

    private void importBRCAinSQL(int viccEntryId, @NotNull KbSpecificObject object) {
        BRCA brca = (BRCA) object;
        context.insertInto(BRCAPART1,
                BRCAPART1.VARIANT_FREQUENCY_LOVD,
                BRCAPART1.ALLELE_FREQUENCY_FIN_EXAC,
                BRCAPART1.CLINVARACCESSION_ENIGMA,
                BRCAPART1.HOMOZYGOUS_COUNT_AFR_EXAC,
                BRCAPART1.BX_ID_EXAC,
                BRCAPART1.VARIANT_IN_LOVD,
                BRCAPART1.ALLELE_FREQUENCY_AFR_EXAC,
                BRCAPART1.DBID_LOVD,
                BRCAPART1.CHR,
                BRCAPART1.BX_ID_ENIGMA,
                BRCAPART1.CO_OCCURRENCE_LR_EXLOVD,
                BRCAPART1.HOMOZYGOUS_COUNT_EAS_EXAC,
                BRCAPART1.SUBMITTER_CLINVAR,
                BRCAPART1.ALLELE_FREQUENCY_EAS_EXAC,
                BRCAPART1.HG37_END,
                BRCAPART1.SUBMITTERS_LOVD,
                BRCAPART1.CLINICAL_CLASSIFICATION_BIC,
                BRCAPART1.HOMOZYGOUS_COUNT_NFE_EXAC,
                BRCAPART1.ALLELE_COUNT_SAS_EXAC,
                BRCAPART1.METHOD_CLINVAR,
                BRCAPART1.ALLELE_COUNT_NFE_EXAC,
                BRCAPART1.PATHOGENICITY_ALL,
                BRCAPART1.GERMLINE_OR_SOMATIC_BIC,
                BRCAPART1.HOMOZYGOUS_COUNT_SAS_EXAC,
                BRCAPART1.BIC_NOMENCLATURE,
                BRCAPART1.ASSERTION_METHOD_ENIGMA,
                BRCAPART1.LITERATURE_SOURCE_EXLOVD,
                BRCAPART1.CHANGE_TYPE_ID,
                BRCAPART1.COLLECTION_METHOD_ENIGMA,
                BRCAPART1.SUM_FAMILY_LR_EXLOVD,
                BRCAPART1.HGVS_CDNA_LOVD,
                BRCAPART1.HOMOZYGOUS_COUNT_FIN_EXAC,
                BRCAPART1.EAS_ALLELE_FREQUENCY_1000_GENOMES,
                BRCAPART1.ETHNICITY_BIC,
                BRCAPART1.INDIVIDUALS_LOVD,
                BRCAPART1.VARIANT_IN_EXAC,
                BRCAPART1.URL_ENIGMA,
                BRCAPART1.ALLELE_ORIGIN_CLINVAR,
                BRCAPART1.ALLELE_FREQUENCY_AMR_EXAC,
                BRCAPART1.VARIANT_IN_1000_GENOMES,
                BRCAPART1.AFR_ALLELE_FREQUENCY_1000_GENOMES,
                BRCAPART1.BX_ID_EXLOVD,
                BRCAPART1.SOURCE,
                BRCAPART1.CONDITION_ID_VALUE_ENIGMA,
                BRCAPART1.HGVS_PROTEIN,
                BRCAPART1.REF,
                BRCAPART1.ALLELE_NUMBER_AFR_EXAC,
                BRCAPART1.ALLELE_COUNT_AFR_EXAC,
                BRCAPART1.BX_ID_LOVD,
                BRCAPART1.SYNONYMS,
                BRCAPART1.GENE_SYMBOL,
                BRCAPART1.COMMENT_ON_CLINICAL_SIGNIFICANCE_ENIGMA,
                BRCAPART1.MISSENSE_ANALYSIS_PRIOR_PROBABILITY_EXLOVD,
                BRCAPART1.ALLELE_NUMBER_FIN_EXAC,
                BRCAPART1.POSTERIOR_PROBABILITY_EXLOVD,
                BRCAPART1.POLYPHEN_SCORE,
                BRCAPART1.REFERENCE_SEQUENCE,
                BRCAPART1.ALLELE_COUNT_EAS_EXAC,
                BRCAPART1.HG38_END,
                BRCAPART1.HGVS_CDNA,
                BRCAPART1.FUNCTIONAL_ANALYSIS_TECHNIQUE_LOVD,
                BRCAPART1.SAS_ALLELE_FREQUENCY_1000_GENOMES,
                BRCAPART1.RNA_LOVD,
                BRCAPART1.COMBINED_PRIOR_PROBABLILITY_EXLOVD,
                BRCAPART1.BX_ID_CLINVAR,
                BRCAPART1.IARC_CLASS_EXLOVD,
                BRCAPART1.BX_ID_BIC,
                BRCAPART1.VICCENTRYID)
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
                        viccEntryId)
                .execute();

        context.insertInto(BRCAPART2,
                BRCAPART2.SIFT_PREDICTION,
                BRCAPART2.ALLELE_NUMBER_NFE_EXAC,
                BRCAPART2.ALLELE_ORIGIN_ENIGMA,
                BRCAPART2.ALLELE_NUMBER_OTH_EXAC,
                BRCAPART2.HG36_END,
                BRCAPART2.ALLELE_FREQUENCY_SAS_EXAC,
                BRCAPART2.DATE_LAST_UPDATED_CLINVAR,
                BRCAPART2.ALLELE_NUMBER_EAS_EXAC,
                BRCAPART2.ALLELE_FREQUENCY_OTH_EXAC,
                BRCAPART2.SOURCE_URL,
                BRCAPART2.SCV_CLINVAR,
                BRCAPART2.PATHOGENICITY_EXPERT,
                BRCAPART2.ALLELE_FREQUENCY_1000_GENOMES,
                BRCAPART2.FUNCTIONAL_ANALYSIS_RESULT_LOVD,
                BRCAPART2.AMR_ALLELE_FREQUENCY_1000_GENOMES,
                BRCAPART2.VARIANT_IN_ESP,
                BRCAPART2.VARIANT_IN_BIC,
                BRCAPART2.CLINICAL_SIGNIFICANCE_ENIGMA,
                BRCAPART2.MAX_ALLELE_FREQUENCY,
                BRCAPART2.ALLELE_COUNT_AMR_EXAC,
                BRCAPART2.VARIANT_IN_ENIGMA,
                BRCAPART2.BX_ID_ESP,
                BRCAPART2.PATIENT_NATIONALITY_BIC,
                BRCAPART2.BX_ID_1000_GENOMES,
                BRCAPART2.GENOMIC_COORDINATE_HG37,
                BRCAPART2.GENOMIC_COORDINATE_HG36,
                BRCAPART2.EUR_ALLELE_FREQUENCY_1000_GENOMES,
                BRCAPART2.NUMBER_OF_FAMILY_MEMBER_CARRYING_MUTATION_BIC,
                BRCAPART2.SEGREGATION_LR_EXLOVD,
                BRCAPART2.ALLELE_FREQUENCY,
                BRCAPART2.MINOR_ALLELE_FREQUENCY_PERCENT_ESP,
                BRCAPART2.ALLELE_FREQUENCY_EXAC,
                BRCAPART2.MUTATION_TYPE_BIC,
                BRCAPART2.ASSERTION_METHOD_CITATION_ENIGMA,
                BRCAPART2.CONDITION_ID_TYPE_ENIGMA,
                BRCAPART2.ALLELE_COUNT_OTH_EXAC,
                BRCAPART2.HGVS_PROTEIN_LOVD,
                BRCAPART2.VARIANT_IN_CLINVAR,
                BRCAPART2.CLINICAL_IMPORTANCE_BIC,
                BRCAPART2.DISCORDANT,
                BRCAPART2.ALLELE_COUNT_FIN_EXAC,
                BRCAPART2.CONDITION_CATEGORY_ENIGMA,
                BRCAPART2.ALLELE_FREQUENCY_ESP,
                BRCAPART2.HOMOZYGOUS_COUNT_OTH_EXAC,
                BRCAPART2.GENETIC_ORIGIN_LOVD,
                BRCAPART2.HOMOZYGOUS_COUNT_AMR_EXAC,
                BRCAPART2.CLINICAL_SIGNIFICANCE_CLINVAR,
                BRCAPART2.AA_ALLELE_FREQUENCY_ESP,
                BRCAPART2.PROTEIN_CHANGE,
                BRCAPART2.VARIANT_IN_EXLOVD,
                BRCAPART2.EA_ALLELE_FREQUENCY_ESP,
                BRCAPART2.HGVS_RNA,
                BRCAPART2.CLINICAL_SIGNIFICANCE_CITATIONS_ENIGMA,
                BRCAPART2.VARIANT_EFFECT_LOVD,
                BRCAPART2.POLYPHEN_PREDICTION,
                BRCAPART2.DATA_RELEASE_ID,
                BRCAPART2.HG37_START,
                BRCAPART2.HG36_START,
                BRCAPART2.SIFT_SCORE,
                BRCAPART2.GENOMIC_COORDINATE_HG38,
                BRCAPART2.ALT,
                BRCAPART2.LITERATURE_CITATION_BIC,
                BRCAPART2.VARIANT_HAPLOTYPE_LOVD,
                BRCAPART2.ALLELE_FREQUENCY_NFE_EXAC,
                BRCAPART2.HG38_START,
                BRCAPART2.POS,
                BRCAPART2.DATE_LAST_EVALUATED_ENIGMA,
                BRCAPART2.ALLELE_NUMBER_SAS_EXAC,
                BRCAPART2.ALLELE_NUMBER_AMR_EXAC,
                BRCAPART2.VICCENTRYID)
                .values(brca.brcApart1().Sift_Prediction(),
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

    private void importCGIinSQL(int viccEntryId, @NotNull KbSpecificObject object) {
        Cgi cgi = (Cgi) object;
        int id = context.insertInto(CGI,
                CGI.TARGETING,
                CGI.SOURCE,
                CGI.PRIMARYTUMORTYPE,
                CGI.DRUGSFULLNAME,
                CGI.CURATOR,
                CGI.DRUGFAMILY,
                CGI.ALTERATION,
                CGI.DRUG,
                CGI.BIOMARKER,
                CGI.DRUGSTATUS,
                CGI.GENE,
                CGI.ASSAYTYPE,
                CGI.ALTERATIONTYPE,
                CGI.EVIDENCELEVEL,
                CGI.ASSOCIATION,
                CGI.METASTATICTUMORTYPE,
                CGI.VICCENTRYID)
                .values(cgi.targeting(),
                        cgi.source(),
                        cgi.primary_tumor_type(),
                        cgi.drugsFullName(),
                        cgi.curator(),
                        cgi.drug_family(),
                        cgi.alteration(),
                        cgi.drug(),
                        cgi.biomarker(),
                        cgi.drug_status(),
                        cgi.gene(),
                        cgi.assay_type(),
                        cgi.alteration_type(),
                        cgi.evidence_level(),
                        cgi.association(),
                        cgi.metastatic_Tumor_Type(),
                        viccEntryId)
                .returning(CGI.ID)
                .fetchOne()
                .getValue(CGI.ID);

        for (String cDNA : cgi.cDNA()) {
            context.insertInto(CGICDNA, CGICDNA.CDNA, CGICDNA.CGIID).values(cDNA, id).execute();
        }

        for (String individualMutation : cgi.individual_mutation()) {
            context.insertInto(CGIINDIVIDUALMUTATION, CGIINDIVIDUALMUTATION.INDIVIDUALMUTATION, CGIINDIVIDUALMUTATION.CGIID)
                    .values(individualMutation, id)
                    .execute();
        }

        for (String gDNA : cgi.gDNA()) {
            context.insertInto(CGIGDNA, CGIGDNA.GDNA, CGIGDNA.CGIID).values(gDNA, id).execute();
        }

        for (String transcript : cgi.transcript()) {
            context.insertInto(CGITRANSCRIPT, CGITRANSCRIPT.TRANSCRIPT, CGITRANSCRIPT.CGIID).values(transcript, id).execute();
        }

        for (String strand : cgi.strand()) {
            context.insertInto(CGISTRAND, CGISTRAND.STRAND, CGISTRAND.CGIID).values(strand, id).execute();
        }

        for (String info : cgi.info()) {
            context.insertInto(CGIINFO, CGIINFO.INFO, CGIINFO.CGIID).values(info, id).execute();
        }

        for (String region : cgi.region()) {
            context.insertInto(CGIREGION, CGIREGION.REGION, CGIREGION.CGIID).values(region, id).execute();
        }
    }

    private void importJAXinSQL(int viccEntryId, @NotNull KbSpecificObject object) {
        Jax jax = (Jax) object;

        int id = context.insertInto(JAX,
                JAX.RESPONSETYPE,
                JAX.APPROVALSTATUS,
                JAX.EVIDENCETYPE,
                JAX.EFFICACYEVIDENCE,
                JAX.IDJAXSOURCE,
                JAX.VICCENTRYID)
                .values(jax.responseType(), jax.approvalStatus(), jax.evidenceType(), jax.efficacyEvidence(), jax.id(), viccEntryId)
                .returning(JAX.ID)
                .fetchOne()
                .getValue(JAX.ID);

        context.insertInto(JAXMOLECULARPROFILE,
                JAXMOLECULARPROFILE.PROFILENAME,
                JAXMOLECULARPROFILE.IDMOLECULARPROFILE,
                JAXMOLECULARPROFILE.JAXID).values(jax.molecularProfile().profileName(), jax.molecularProfile().id(), id).execute();

        context.insertInto(JAXTHERAPY, JAXTHERAPY.THERAPYNAME, JAXTHERAPY.IDTHERAPY, JAXTHERAPY.JAXID)
                .values(jax.therapy().therapyName(), jax.therapy().id(), id)
                .execute();

        context.insertInto(JAXINDICATIONS, JAXINDICATIONS.SOURCE, JAXINDICATIONS.IDINDICATIONS, JAXINDICATIONS.NAME, JAXINDICATIONS.JAXID)
                .values(jax.indications().source(), jax.indications().id(), jax.indications().name(), id)
                .execute();

        for (JaxReferences references : jax.references()) {
            context.insertInto(JAXREFERENCES,
                    JAXREFERENCES.URL,
                    JAXREFERENCES.IDREFERENCES,
                    JAXREFERENCES.PUBMEDID,
                    JAXREFERENCES.TITLE,
                    JAXREFERENCES.JAXID).values(references.url(), references.id(), references.pubMedId(), references.title(), id).execute();
        }
    }

    private void importJaxTrialsinSQL(int viccEntryId, @NotNull KbSpecificObject object) {
        JaxTrials jaxTrials = (JaxTrials) object;

        int id = context.insertInto(JAXTRIALS,
                JAXTRIALS.TITLE,
                JAXTRIALS.GENDER,
                JAXTRIALS.NCTID,
                JAXTRIALS.SPONSORS,
                JAXTRIALS.RECRUITMENT,
                JAXTRIALS.VARIANTREQUIREMENTS,
                JAXTRIALS.UPDATEDATE,
                JAXTRIALS.PHASE,
                JAXTRIALS.VICCENTRYID)
                .values(jaxTrials.title(),
                        jaxTrials.gender(),
                        jaxTrials.nctId(),
                        jaxTrials.sponsors(),
                        jaxTrials.recruitment(),
                        jaxTrials.variantRequirements(),
                        jaxTrials.updateDate(),
                        jaxTrials.updateDate(),
                        viccEntryId)
                .returning(JAXTRIALS.ID)
                .fetchOne()
                .getValue(JAXTRIALS.ID);

        for (JaxTrialsIndications indications : jaxTrials.indications()) {
            context.insertInto(JAXTRIALSINDICATIONS,
                    JAXTRIALSINDICATIONS.SOURCE,
                    JAXTRIALSINDICATIONS.IDINDICATIONS,
                    JAXTRIALSINDICATIONS.NAME,
                    JAXTRIALSINDICATIONS.JAXTRIALSID).values(indications.source(), indications.id(), indications.name(), id).execute();
        }

        for (JaxTrialsVariantRequirementDetails variantRequirementDetails : jaxTrials.variantRequirementDetails()) {
            int id1 = context.insertInto(JAXTRIALSVARIANTREQUIREMENTDETAILS,
                    JAXTRIALSVARIANTREQUIREMENTDETAILS.REQUIREMENTTYPE,
                    JAXTRIALSVARIANTREQUIREMENTDETAILS.JAXTRIALSID)
                    .values(variantRequirementDetails.requirementType(), id)
                    .returning(JAXTRIALSVARIANTREQUIREMENTDETAILS.ID)
                    .fetchOne()
                    .getValue(JAXTRIALSVARIANTREQUIREMENTDETAILS.ID);

            for (JaxTrialsMolecularProfile molecularProfile : variantRequirementDetails.molecularProfiles()) {
                context.insertInto(JAXTRIALSMOLECULARPROFILE,
                        JAXTRIALSMOLECULARPROFILE.PROFILENAME,
                        JAXTRIALSMOLECULARPROFILE.IDMOLECULARPROFILE,
                        JAXTRIALSMOLECULARPROFILE.JAXTRIALSVARIANTREQUIREMENTDETAILSID)
                        .values(molecularProfile.profileName(), molecularProfile.id(), id1)
                        .execute();
            }
        }

        for (JaxTrialsTherapies therapies : jaxTrials.therapies()) {
            context.insertInto(JAXTRIALSTHERAPIES,
                    JAXTRIALSTHERAPIES.IDTHERAPIES,
                    JAXTRIALSTHERAPIES.THERAPYNAME,
                    JAXTRIALSTHERAPIES.JAXTRIALSID).values(therapies.id(), therapies.therapyName(), id).execute();
        }

    }

    private void importPmkbinSQL(int viccEntryId, @NotNull KbSpecificObject object) {
        Pmkb pmkb = (Pmkb) object;

        int id = context.insertInto(PMKB, PMKB.VICCENTRYID).values(viccEntryId).returning(PMKB.ID).fetchOne().getValue(PMKB.ID);

        for (PmkbTumor tumor : pmkb.tumor()) {
            context.insertInto(PMKBTUMOR, PMKBTUMOR.IDTUMOR, PMKBTUMOR.NAME, PMKBTUMOR.VICCENTRYID)
                    .values(tumor.id(), tumor.name(), id)
                    .execute();
        }

        for (PmkbTissue tissue : pmkb.tissue()) {
            context.insertInto(PMKBTISSUE, PMKBTISSUE.IDTISSUE, PMKBTISSUE.NAME, PMKBTISSUE.VICCENTRYID)
                    .values(tissue.id(), tissue.name(), id)
                    .execute();
        }

        for (PmkbVariant variant : pmkb.variant()) {
            int idVariant = context.insertInto(PMKBVARIANT,
                    PMKBVARIANT.AMINOACIDCHANGE,
                    PMKBVARIANT.GERMLINE,
                    PMKBVARIANT.PARTNERGENE,
                    PMKBVARIANT.CODONS,
                    PMKBVARIANT.DESCRIPTION,
                    PMKBVARIANT.EXONS,
                    PMKBVARIANT.NOTES,
                    PMKBVARIANT.COSMIC,
                    PMKBVARIANT.EFFECT,
                    PMKBVARIANT.CNVTYPE,
                    PMKBVARIANT.IDVARIANT,
                    PMKBVARIANT.CYTOBAND,
                    PMKBVARIANT.VARIANTTYPE,
                    PMKBVARIANT.DNACHANGE,
                    PMKBVARIANT.COORDINATES,
                    PMKBVARIANT.CHROMOSOMEBASEDCNV,
                    PMKBVARIANT.TRANSCRIPT,
                    PMKBVARIANT.DESCRIPTIONTYPE,
                    PMKBVARIANT.CHROMOSOME,
                    PMKBVARIANT.NAME,
                    PMKBVARIANT.VICCENTRYID)
                    .values(variant.aminoAcidChange(),
                            variant.germline(),
                            variant.partnerGene(),
                            variant.codons(),
                            variant.description(),
                            variant.exons(),
                            variant.notes(),
                            variant.cosmic(),
                            variant.effect(),
                            variant.cnvType(),
                            variant.id(),
                            variant.cytoband(),
                            variant.variantType(),
                            variant.dnaChange(),
                            variant.coordinates(),
                            variant.chromosomeBasedCnv(),
                            variant.transcript(),
                            variant.descriptionType(),
                            variant.chromosome(),
                            variant.name(),
                            viccEntryId)
                    .returning(PMKBVARIANT.ID)
                    .fetchOne()
                    .getValue(PMKBVARIANT.ID);

            for (PmkbGene gene : variant.gene()) {
                context.insertInto(PMKBGENE,
                        PMKBGENE.DESCRIPTION,
                        PMKBGENE.CREATEDAT,
                        PMKBGENE.UPDATEDAT,
                        PMKBGENE.ACTIVEIND,
                        PMKBGENE.EXTERNALID,
                        PMKBGENE.IDGENE,
                        PMKBGENE.NAME,
                        PMKBGENE.PMKBVARIANTID)
                        .values(gene.description(),
                                gene.createdAt(),
                                gene.updatedAt(),
                                gene.activeInd(),
                                gene.externalId(),
                                gene.id(),
                                gene.name(),
                                idVariant)
                        .execute();
            }
        }
    }

    private void importOncokbinSQL(int viccEntryId, @NotNull KbSpecificObject object) {
        Oncokb oncokb = (Oncokb) object;
        int id = context.insertInto(ONCOKB, ONCOKB.VICCENTRYID).values(viccEntryId).returning(ONCOKB.ID).fetchOne().getValue(ONCOKB.ID);

        OncoKbBiological oncokbBiological = oncokb.oncoKbBiological();
        if (oncokbBiological != null) {
            int idBiological = context.insertInto(ONCOKBBIOLOGICAL,
                    ONCOKBBIOLOGICAL.MUTATIONEFFECTPMIDS,
                    ONCOKBBIOLOGICAL.ISOFORM,
                    ONCOKBBIOLOGICAL.ENTREZGENEID,
                    ONCOKBBIOLOGICAL.ONCOGENIC,
                    ONCOKBBIOLOGICAL.MUTATIONEFFECT,
                    ONCOKBBIOLOGICAL.REFSEQ,
                    ONCOKBBIOLOGICAL.GENE,
                    ONCOKBBIOLOGICAL.MUTATIONEFFECTABSTRACTS,
                    ONCOKBBIOLOGICAL.VICCENTRYID)
                    .values(oncokbBiological.mutationEffectPmids(),
                            oncokbBiological.Isoform(),
                            oncokbBiological.entrezGeneID(),
                            oncokbBiological.oncogenic(),
                            oncokbBiological.mutationEffect(),
                            oncokbBiological.RefSeq(),
                            oncokbBiological.gene(),
                            oncokbBiological.mutationEffectAbstracts(),
                            id)
                    .returning(ONCOKBBIOLOGICAL.ID)
                    .fetchOne()
                    .getValue(ONCOKBBIOLOGICAL.ID);

            OncokbVariant oncokbVariant = oncokb.oncoKbBiological().oncokbVariant();

            int idVariant = context.insertInto(ONCOKBVARIANTBIOLOGICAL,
                    ONCOKBVARIANTBIOLOGICAL.VARIANTRESIDUES,
                    ONCOKBVARIANTBIOLOGICAL.PROTEINSTART,
                    ONCOKBVARIANTBIOLOGICAL.NAME,
                    ONCOKBVARIANTBIOLOGICAL.PROTEINEND,
                    ONCOKBVARIANTBIOLOGICAL.REFRESIDUES,
                    ONCOKBVARIANTBIOLOGICAL.ALTERATION,
                    ONCOKBVARIANTBIOLOGICAL.ONCOKBBIOLOGICALID)
                    .values(oncokbVariant.variantResidues(),
                            oncokbVariant.proteinStart(),
                            oncokbVariant.name(),
                            oncokbVariant.proteinEnd(),
                            oncokbVariant.refResidues(),
                            oncokbVariant.alteration(),
                            idBiological)
                    .returning(ONCOKBVARIANTBIOLOGICAL.ID)
                    .fetchOne()
                    .getValue(ONCOKBVARIANTBIOLOGICAL.ID);

            OncoKbConsequence oncoKbConsequence = oncokb.oncoKbBiological().oncokbVariant().oncoKbConsequence();

            context.insertInto(ONCOKBCONSEQUENCESBIOLOGICAL,
                    ONCOKBCONSEQUENCESBIOLOGICAL.TERM,
                    ONCOKBCONSEQUENCESBIOLOGICAL.DESCRIPTION,
                    ONCOKBCONSEQUENCESBIOLOGICAL.ISGENERALLYTRUNCATING,
                    ONCOKBCONSEQUENCESBIOLOGICAL.ONCOKBVARIANTBIOLOGICALID)
                    .values(oncoKbConsequence.term(), oncoKbConsequence.description(), oncoKbConsequence.isGenerallyTruncating(), idVariant)
                    .execute();

            OncokbGene oncoKbGene = oncokb.oncoKbBiological().oncokbVariant().oncokbGene();

            int idGene = context.insertInto(ONCOKBGENEBIOLOGICAL,
                    ONCOKBGENEBIOLOGICAL.ONCOGENE,
                    ONCOKBGENEBIOLOGICAL.NAME,
                    ONCOKBGENEBIOLOGICAL.HUGOSYMBOL,
                    ONCOKBGENEBIOLOGICAL.CURATEDREFSEQ,
                    ONCOKBGENEBIOLOGICAL.ENTREZGENEID,
                    ONCOKBGENEBIOLOGICAL.TSG,
                    ONCOKBGENEBIOLOGICAL.CURATEDISOFORM,
                    ONCOKBGENEBIOLOGICAL.ONCOKBBIOLOGICALID)
                    .values(oncoKbGene.oncogene(),
                            oncoKbGene.name(),
                            oncoKbGene.hugoSymbol(),
                            oncoKbGene.curatedRefSeq(),
                            oncoKbGene.entrezGeneId(),
                            oncoKbGene.tsg(),
                            oncoKbGene.curatedIsoform(),
                            idVariant)
                    .returning(ONCOKBGENEBIOLOGICAL.ID)
                    .fetchOne()
                    .getValue(ONCOKBGENEBIOLOGICAL.ID);

            for (String geneAliases : oncoKbGene.geneAliases()) {
                context.insertInto(ONCOKBGENEALIASESBIOLOGICAL,
                        ONCOKBGENEALIASESBIOLOGICAL.GENEALIASES,
                        ONCOKBGENEALIASESBIOLOGICAL.ONCOKBGENEBIOLOGICALID).values(geneAliases, idGene);
            }
        }

        OncoKbClinical oncokbClinical = oncokb.oncoKbClinical();
        if (oncokbClinical != null) {

            int idClinical = context.insertInto(ONCOKBCLINICAL,
                    ONCOKBCLINICAL.REFSEQ,
                    ONCOKBCLINICAL.LEVEL,
                    ONCOKBCLINICAL.ENTREZGENEID,
                    ONCOKBCLINICAL.DRUGPMIDS,
                    ONCOKBCLINICAL.CANCERTYPE,
                    ONCOKBCLINICAL.DRUG,
                    ONCOKBCLINICAL.GENE,
                    ONCOKBCLINICAL.LEVELLABEL,
                    ONCOKBCLINICAL.VICCENTRYID)
                    .values(oncokbClinical.RefSeq(),
                            oncokbClinical.level(),
                            oncokbClinical.entrezGeneID(),
                            oncokbClinical.drugPmids(),
                            oncokbClinical.cancerType(),
                            oncokbClinical.drug(),
                            oncokbClinical.gene(),
                            oncokbClinical.levelLabel(),
                            id)
                    .returning(ONCOKBCLINICAL.ID)
                    .fetchOne()
                    .getValue(ONCOKBCLINICAL.ID);

            for (OncoKbDrugAbstracts drugAbstracts : oncokbClinical.oncoKbDrugAbstracts()) {
                context.insertInto(ONCOKBDRUGABSTRACTSCLINICAL,
                        ONCOKBDRUGABSTRACTSCLINICAL.TEXT,
                        ONCOKBDRUGABSTRACTSCLINICAL.LINK,
                        ONCOKBDRUGABSTRACTSCLINICAL.ONCOKBCLINICALID)
                        .values(drugAbstracts.text(), drugAbstracts.link(), idClinical)
                        .execute();
            }

            OncokbVariant variantClinical = oncokbClinical.oncokbVariant();
            int idClinicalVariant = context.insertInto(ONCOKBVARIANTCLINICAL,
                    ONCOKBVARIANTCLINICAL.VARIANTRESIDUES,
                    ONCOKBVARIANTCLINICAL.PROTEINSTART,
                    ONCOKBVARIANTCLINICAL.NAME,
                    ONCOKBVARIANTCLINICAL.PROTEINEND,
                    ONCOKBVARIANTCLINICAL.REFRESIDUES,
                    ONCOKBVARIANTCLINICAL.ALTERATION,
                    ONCOKBVARIANTCLINICAL.ONCOKBCLINICALID)
                    .values(variantClinical.variantResidues(),
                            variantClinical.proteinStart(),
                            variantClinical.name(),
                            variantClinical.proteinEnd(),
                            variantClinical.refResidues(),
                            variantClinical.alteration(),
                            idClinical)
                    .returning(ONCOKBVARIANTCLINICAL.ID)
                    .fetchOne()
                    .getValue(ONCOKBVARIANTCLINICAL.ID);

            OncoKbConsequence consequenceaClinical = oncokbClinical.oncokbVariant().oncoKbConsequence();

            context.insertInto(ONCOKBCONSEQUENCESCLINICAL,
                    ONCOKBCONSEQUENCESCLINICAL.TERM,
                    ONCOKBCONSEQUENCESCLINICAL.DESCRIPTION,
                    ONCOKBCONSEQUENCESCLINICAL.ISGENERALLYTRUNCATING,
                    ONCOKBCONSEQUENCESCLINICAL.ONCOKBVARIANTCLINICALID)
                    .values(consequenceaClinical.term(),
                            consequenceaClinical.description(),
                            consequenceaClinical.isGenerallyTruncating(),
                            idClinicalVariant)
                    .execute();

            OncokbGene geneClinical = oncokbClinical.oncokbVariant().oncokbGene();

            int idGeneClinical = context.insertInto(ONCOKBGENECLINICAL,
                    ONCOKBGENECLINICAL.ONCOGENE,
                    ONCOKBGENECLINICAL.NAME,
                    ONCOKBGENECLINICAL.HUGOSYMBOL,
                    ONCOKBGENECLINICAL.CURATEDREFSEQ,
                    ONCOKBGENECLINICAL.ENTREZGENEID,
                    ONCOKBGENECLINICAL.TSG,
                    ONCOKBGENECLINICAL.CURATEDISOFORM,
                    ONCOKBGENECLINICAL.ONCOKBCLINICALID)
                    .values(geneClinical.oncogene(),
                            geneClinical.name(),
                            geneClinical.hugoSymbol(),
                            geneClinical.curatedRefSeq(),
                            geneClinical.entrezGeneId(),
                            geneClinical.tsg(),
                            geneClinical.curatedIsoform(),
                            idClinicalVariant)
                    .returning(ONCOKBGENECLINICAL.ID)
                    .fetchOne()
                    .getValue(ONCOKBGENECLINICAL.ID);

            for (String geneAliasesClinical : geneClinical.geneAliases()) {
                context.insertInto(ONCOKBGENEALIASESCLINICAL,
                        ONCOKBGENEALIASESCLINICAL.GENEALIASES,
                        ONCOKBGENEALIASESCLINICAL.ONCOKBGENECLINICALID).values(geneAliasesClinical, idGeneClinical).execute();
            }

        }

    }

    private void importMolecularMatchTrials(int viccEntryId, @NotNull KbSpecificObject object) {
        MolecularMatchTrials molecularMatchTrials = (MolecularMatchTrials) object;
        int id = context.insertInto(MOLECULARMATCHTRIALS,
                MOLECULARMATCHTRIALS.STATUS,
                MOLECULARMATCHTRIALS.STARTDATE,
                MOLECULARMATCHTRIALS.TITLE,
                MOLECULARMATCHTRIALS.SCORE,
                MOLECULARMATCHTRIALS.BRIEFTITLE,
                MOLECULARMATCHTRIALS.LINK,
                MOLECULARMATCHTRIALS.PHASE,
                MOLECULARMATCHTRIALS.IDMOLECULARMATCHTRIALS,
                MOLECULARMATCHTRIALS.STUDYTYPE,
                MOLECULARMATCHTRIALS.VICCENTRYID)
                .values(molecularMatchTrials.status(),
                        molecularMatchTrials.startDate(),
                        molecularMatchTrials.title(),
                        molecularMatchTrials.score(),
                        molecularMatchTrials.briefTitle(),
                        molecularMatchTrials.link(),
                        molecularMatchTrials.phase(),
                        molecularMatchTrials.id(),
                        molecularMatchTrials.studyType(),
                        viccEntryId)
                .returning(MOLECULARMATCHTRIALS.ID)
                .fetchOne()
                .getValue(MOLECULARMATCHTRIALS.ID);

        for (String molecularAlterations : molecularMatchTrials.molecularAlterations()) {
            context.insertInto(MOLECULARMATCHTRIALSALTERATIONS,
                    MOLECULARMATCHTRIALSALTERATIONS.MOLECULARALTERATIONS,
                    MOLECULARMATCHTRIALSALTERATIONS.MOLECULARMATCHTRIALSID).values(molecularAlterations, id).execute();
        }

        for (MolecularMatchTrialsIntervation intervation : molecularMatchTrials.intervation()) {
            int idIntervation = context.insertInto(MOLECULARMATCHTRIALSINTERVATIONS,
                    MOLECULARMATCHTRIALSINTERVATIONS.INTERVENTION_NAME,
                    MOLECULARMATCHTRIALSINTERVATIONS.DESCRIPTION,
                    MOLECULARMATCHTRIALSINTERVATIONS.INTERVENTION_TYPE,
                    MOLECULARMATCHTRIALSINTERVATIONS.MOLECULARMATCHTRIALSID)
                    .values(intervation.intervention_name(), intervation.description(), intervation.intervention_type(), viccEntryId)
                    .returning(MOLECULARMATCHTRIALSINTERVATIONS.ID)
                    .fetchOne()
                    .getValue(MOLECULARMATCHTRIALSINTERVATIONS.ID);

            if (intervation.other_name() != null) {
                for (String otherName : intervation.other_name()) {
                    context.insertInto(MOLECULARMATCHTRIALSOTHERNAME,
                            MOLECULARMATCHTRIALSOTHERNAME.OTHER_NAME,
                            MOLECULARMATCHTRIALSOTHERNAME.MOLECULARMATCHTRIALSINTERVATIONSID).values(otherName, idIntervation).execute();
                }
            }

            if (intervation.arm_group_label() != null) {
                for (String armGroupLabel : intervation.arm_group_label()) {
                    context.insertInto(MOLECULARMATCHTRIALSOTHERGROUPLABEL,
                            MOLECULARMATCHTRIALSOTHERGROUPLABEL.ARM_GROUP_LABEL,
                            MOLECULARMATCHTRIALSOTHERGROUPLABEL.MOLECULARMATCHTRIALSINTERVATIONSID)
                            .values(armGroupLabel, idIntervation)
                            .execute();
                }
            }
        }

        for (MolecularMatchTrialsLocations locations : molecularMatchTrials.locations()) {
            int idLocations = context.insertInto(MOLECULARMATCHTRIALSLOCATIONS,
                    MOLECULARMATCHTRIALSLOCATIONS.STATUS,
                    MOLECULARMATCHTRIALSLOCATIONS.LAST_NAME,
                    MOLECULARMATCHTRIALSLOCATIONS.EMAIL,
                    MOLECULARMATCHTRIALSLOCATIONS.PHONE,
                    MOLECULARMATCHTRIALSLOCATIONS.PHONE_BACKUP,
                    MOLECULARMATCHTRIALSLOCATIONS.EMAIL_BACKUP,
                    MOLECULARMATCHTRIALSLOCATIONS.LAST_NAME_BACKUP,
                    MOLECULARMATCHTRIALSLOCATIONS.PHONE_EXT_BACKUP,
                    MOLECULARMATCHTRIALSLOCATIONS.PHONE_EXT,
                    MOLECULARMATCHTRIALSLOCATIONS.CITY,
                    MOLECULARMATCHTRIALSLOCATIONS.VALID,
                    MOLECULARMATCHTRIALSLOCATIONS.ZIP,
                    MOLECULARMATCHTRIALSLOCATIONS.CREATED,
                    MOLECULARMATCHTRIALSLOCATIONS.COUNTRY,
                    MOLECULARMATCHTRIALSLOCATIONS.NUMBER,
                    MOLECULARMATCHTRIALSLOCATIONS.IDLOCATIONS,
                    MOLECULARMATCHTRIALSLOCATIONS.LASTUPDATED,
                    MOLECULARMATCHTRIALSLOCATIONS.STATE,
                    MOLECULARMATCHTRIALSLOCATIONS.STREET,
                    MOLECULARMATCHTRIALSLOCATIONS.PO_BOX,
                    MOLECULARMATCHTRIALSLOCATIONS.FAILEDGEOCODE,
                    MOLECULARMATCHTRIALSLOCATIONS.VALIDMESSAGE,
                    MOLECULARMATCHTRIALSLOCATIONS.NAME,
                    MOLECULARMATCHTRIALSLOCATIONS.MOLECULARMATCHTRIALSID)
                    .values(locations.status(),
                            locations.last_name(),
                            locations.email(),
                            locations.phone(),
                            locations.phone_backup(),
                            locations.email_backup(),
                            locations.last_name_backup(),
                            locations.phone_ext_backup(),
                            locations.phone_ext(),
                            locations.city(),
                            locations.valid(),
                            locations.zip(),
                            locations.created(),
                            locations.country(),
                            locations.number(),
                            locations.id(),
                            locations.lastUpdated(),
                            locations.state(),
                            locations.street(),
                            locations.po_box(),
                            locations.failedGeocode(),
                            locations.validMessage(),
                            locations.name(),
                            id)
                    .returning(MOLECULARMATCHTRIALSLOCATIONS.ID)
                    .fetchOne()
                    .getValue(MOLECULARMATCHTRIALSLOCATIONS.ID);

            if (locations.contact() != null) {
                context.insertInto(MOLECULARMATCHTRIALSCONTACT,
                        MOLECULARMATCHTRIALSCONTACT.PHONE,
                        MOLECULARMATCHTRIALSCONTACT.NAME,
                        MOLECULARMATCHTRIALSCONTACT.EMAIL,
                        MOLECULARMATCHTRIALSCONTACT.MOLECULARMATCHTRIALSLOCATIONSID)
                        .values(locations.contact().phone(), locations.contact().name(), locations.contact().email(), idLocations)
                        .execute();
            }

            if (locations.location() != null) {
                int idLocation = context.insertInto(MOLECULARMATCHTRIALSLOCATION,
                        MOLECULARMATCHTRIALSLOCATION.TYPE,
                        MOLECULARMATCHTRIALSLOCATION.MOLECULARMATCHTRIALSLOCATIONSID)
                        .values(locations.location().type(), idLocations)
                        .returning(MOLECULARMATCHTRIALSLOCATION.ID)
                        .fetchOne()
                        .getValue(MOLECULARMATCHTRIALSLOCATION.ID);

                for (String coordinates : locations.location().coordinates()) {
                    context.insertInto(MOLECULARMATCHTRIALSCOORDINATES,
                            MOLECULARMATCHTRIALSCOORDINATES.COORDINATES,
                            MOLECULARMATCHTRIALSCOORDINATES.MOLECULARMATCHTRIALSLOCATIONID).values(coordinates, idLocation).execute();
                }
            }

            if (locations.geo() != null) {
                context.insertInto(MOLECULARMATCHTRIALSGEO,
                        MOLECULARMATCHTRIALSGEO.LAT,
                        MOLECULARMATCHTRIALSGEO.LON,
                        MOLECULARMATCHTRIALSGEO.MOLECULARMATCHTRIALSLOCATIONSID)
                        .values(locations.geo().lat(), locations.geo().lon(), idLocations)
                        .execute();
            }
        }

        if (molecularMatchTrials.overallContact() != null) {
            context.insertInto(MOLECULARMATCHTRIALSOVERALLCONTACT,
                    MOLECULARMATCHTRIALSOVERALLCONTACT.PHONE,
                    MOLECULARMATCHTRIALSOVERALLCONTACT.LAST_NAME,
                    MOLECULARMATCHTRIALSOVERALLCONTACT.PHONE_EXT,
                    MOLECULARMATCHTRIALSOVERALLCONTACT.COUNTRY,
                    MOLECULARMATCHTRIALSOVERALLCONTACT.EMAIL,
                    MOLECULARMATCHTRIALSOVERALLCONTACT.AFFILIATION,
                    MOLECULARMATCHTRIALSOVERALLCONTACT.CITY,
                    MOLECULARMATCHTRIALSOVERALLCONTACT.NAME,
                    MOLECULARMATCHTRIALSOVERALLCONTACT.ZIP,
                    MOLECULARMATCHTRIALSOVERALLCONTACT.URL,
                    MOLECULARMATCHTRIALSOVERALLCONTACT.STREET,
                    MOLECULARMATCHTRIALSOVERALLCONTACT.TYPE,
                    MOLECULARMATCHTRIALSOVERALLCONTACT.MOLECULARMATCHTRIALSID)
                    .values(molecularMatchTrials.overallContact().phone(),
                            molecularMatchTrials.overallContact().last_name(),
                            molecularMatchTrials.overallContact().phone_ext(),
                            molecularMatchTrials.overallContact().country(),
                            molecularMatchTrials.overallContact().email(),
                            molecularMatchTrials.overallContact().affiliation(),
                            molecularMatchTrials.overallContact().city(),
                            molecularMatchTrials.overallContact().name(),
                            molecularMatchTrials.overallContact().zip(),
                            molecularMatchTrials.overallContact().url(),
                            molecularMatchTrials.overallContact().street(),
                            molecularMatchTrials.overallContact().type(),
                            id)
                    .execute();
        }

        for (MolecularMatchTrialsTags tags : molecularMatchTrials.tags()) {
            context.insertInto(MOLECULARMATCHTRIALSTAGS,
                    MOLECULARMATCHTRIALSTAGS.FACET,
                    MOLECULARMATCHTRIALSTAGS.COMPOSITEKEY,
                    MOLECULARMATCHTRIALSTAGS.SUPPRESS,
                    MOLECULARMATCHTRIALSTAGS.FILTERTYPE,
                    MOLECULARMATCHTRIALSTAGS.TERM,
                    MOLECULARMATCHTRIALSTAGS.CUSTOM,
                    MOLECULARMATCHTRIALSTAGS.PRIORITY,
                    MOLECULARMATCHTRIALSTAGS.ALIAS,
                    MOLECULARMATCHTRIALSTAGS.MANUALSUPPRESS,
                    MOLECULARMATCHTRIALSTAGS.GENERATEDBY,
                    MOLECULARMATCHTRIALSTAGS.GENERATEDBYTERM,
                    MOLECULARMATCHTRIALSTAGS.IDTAGS,
                    MOLECULARMATCHTRIALSTAGS.MANUALPRIORITY,
                    MOLECULARMATCHTRIALSTAGS.MOLECULARMATCHTRIALSID)
                    .values(tags.facet(),
                            tags.compositeKey(),
                            tags.suppress(),
                            tags.filterType(),
                            tags.term(),
                            tags.custom(),
                            tags.priority(),
                            tags.alias(),
                            tags.manualSuppress(),
                            tags.generatedBy(),
                            tags.generatedByTerm(),
                            tags.id(),
                            tags.manualPriority(),
                            id)
                    .execute();
        }

    }

    private void importCivic(int viccEntryId, @NotNull KbSpecificObject object) {
        Civic civic = (Civic) object;

        int id = context.insertInto(CIVIC,
                CIVIC.ENTREZNAME,
                CIVIC.CIVICACTIONABILITYSCORE,
                CIVIC.ALLELEREGISTRYID,
                CIVIC.GENEID,
                CIVIC.NAME,
                CIVIC.ENTREZID,
                CIVIC.TYPE,
                CIVIC.IDCIVIC,
                CIVIC.DESCRIPTION,
                CIVIC.VICCENTRYID)
                .values(civic.entrezName(),
                        civic.civicActionabilityScore(),
                        civic.alleleRegistryId(),
                        civic.geneId(),
                        civic.name(),
                        civic.entrezId(),
                        civic.type(),
                        civic.id(),
                        civic.description(),
                        viccEntryId)
                .returning(CIVIC.ID)
                .fetchOne()
                .getValue(CIVIC.ID);

        for (String assertions : civic.assertions()) {
            context.insertInto(CIVICASSERTIONS, CIVICASSERTIONS.ASSERTIONS, CIVICASSERTIONS.CIVICID).values(assertions, id).execute();
        }

        for (String hgvsExpression : civic.hgvs_expressions()) {
            context.insertInto(CIVICHGVSEXPRESSIONS, CIVICHGVSEXPRESSIONS.HGVS_EXPRESSIONS, CIVICHGVSEXPRESSIONS.CIVICID)
                    .values(hgvsExpression, id)
                    .execute();
        }

        for (String clinvarEntries : civic.clinvarEntries()) {
            context.insertInto(CIVICCLINVARENTRIES, CIVICCLINVARENTRIES.CLINVARENTRIES, CIVICCLINVARENTRIES.CIVICID)
                    .values(clinvarEntries, id)
                    .execute();
        }

        for (String variantAliases : civic.variantAliases()) {
            context.insertInto(CIVICVARIANTALIASES, CIVICVARIANTALIASES.VARIANTALIASES, CIVICVARIANTALIASES.CIVICID)
                    .values(variantAliases, id)
                    .execute();
        }

        for (CivicVariantTypes variantTypes : civic.variantTypes()) {
            context.insertInto(CIVICVARIANTTYPES,
                    CIVICVARIANTTYPES.DISPLAYNAME,
                    CIVICVARIANTTYPES.DESCRIPTION,
                    CIVICVARIANTTYPES.URL,
                    CIVICVARIANTTYPES.SOID,
                    CIVICVARIANTTYPES.IDVARIANTTYPES,
                    CIVICVARIANTTYPES.NAME,
                    CIVICVARIANTTYPES.CIVICID)
                    .values(variantTypes.displayName(),
                            variantTypes.description(),
                            variantTypes.url(),
                            variantTypes.soId(),
                            variantTypes.id(),
                            variantTypes.name(),
                            id)
                    .execute();
        }

        if (civic.provisional_values() != null) {
            context.insertInto(CIVICDESCRIPTION, CIVICDESCRIPTION.REVISIONID, CIVICDESCRIPTION.VALUE, CIVICDESCRIPTION.CIVICID)
                    .values(civic.provisional_values().revision_id(), civic.provisional_values().value(), id)
                    .execute();
        }

        context.insertInto(CIVICCOORDINATES,
                CIVICCOORDINATES.CHROMOSOME2,
                CIVICCOORDINATES.REFERENCEBASES,
                CIVICCOORDINATES.START2,
                CIVICCOORDINATES.VARIANTBASES,
                CIVICCOORDINATES.STOP,
                CIVICCOORDINATES.STOP2,
                CIVICCOORDINATES.REPRESENTATIVETRANSCRIPT2,
                CIVICCOORDINATES.START,
                CIVICCOORDINATES.REPRESENTATIVETRANSCRIPT,
                CIVICCOORDINATES.ENSEMBLVERSION,
                CIVICCOORDINATES.CHROMOSOME,
                CIVICCOORDINATES.REFERENCEBUILD,
                CIVICCOORDINATES.CIVICID)
                .values(civic.coordinates().chromosome2(),
                        civic.coordinates().referenceBases(),
                        civic.coordinates().start2(),
                        civic.coordinates().variantBases(),
                        civic.coordinates().stop(),
                        civic.coordinates().stop2(),
                        civic.coordinates().representativeTranscript2(),
                        civic.coordinates().start(),
                        civic.coordinates().representativeTranscript(),
                        civic.coordinates().ensemblVersion(),
                        civic.coordinates().chromosome(),
                        civic.coordinates().referenceBuild(),
                        id)
                .execute();

        for (CivicVariantGroup variantGroup : civic.variantGroups()) {
            int idVariantGroup = context.insertInto(CIVICVARIANTSGROUPS,
                    CIVICVARIANTSGROUPS.IDVARIANTSGROUPS,
                    CIVICVARIANTSGROUPS.TYPE,
                    CIVICVARIANTSGROUPS.DESCRIPTION,
                    CIVICVARIANTSGROUPS.NAME,
                    CIVICVARIANTSGROUPS.CIVICID)
                    .values(variantGroup.id(), variantGroup.type(), variantGroup.description(), variantGroup.name(), id)
                    .returning(CIVICVARIANTSGROUPS.ID)
                    .fetchOne()
                    .getValue(CIVICVARIANTSGROUPS.ID);

            for (CivicVariants variants : variantGroup.variants()) {
                int idVariantGroupVariants = context.insertInto(CIVICVARIANTSGROUPSVARIANTS,
                        CIVICVARIANTSGROUPSVARIANTS.ENTREZ_NAME,
                        CIVICVARIANTSGROUPSVARIANTS.DESCRIPTION,
                        CIVICVARIANTSGROUPSVARIANTS.CIVIC_ACTIONABILITY_SCORE,
                        CIVICVARIANTSGROUPSVARIANTS.GENE_ID,
                        CIVICVARIANTSGROUPSVARIANTS.ENTREZ_ID,
                        CIVICVARIANTSGROUPSVARIANTS.TYPE,
                        CIVICVARIANTSGROUPSVARIANTS.IDVARIANTS,
                        CIVICVARIANTSGROUPSVARIANTS.NAME,
                        CIVICVARIANTSGROUPSVARIANTS.CIVICVARIANTSGROUPSID)
                        .values(variants.entrez_name(),
                                variants.description(),
                                variants.civic_actionability_score(),
                                variants.gene_id(),
                                variants.entrez_id(),
                                variants.type(),
                                variants.id(),
                                variants.name(),
                                idVariantGroup)
                        .returning(CIVICVARIANTSGROUPSVARIANTS.ID)
                        .fetchOne()
                        .getValue(CIVICVARIANTSGROUPSVARIANTS.ID);

                if (variants.coordinates() != null) {
                    context.insertInto(CIVICVARIANTSGROUPSCOORDINATES,
                            CIVICVARIANTSGROUPSCOORDINATES.CHROMOSOME2,
                            CIVICVARIANTSGROUPSCOORDINATES.REFERENCEBASES,
                            CIVICVARIANTSGROUPSCOORDINATES.START2,
                            CIVICVARIANTSGROUPSCOORDINATES.VARIANTBASES,
                            CIVICVARIANTSGROUPSCOORDINATES.STOP,
                            CIVICVARIANTSGROUPSCOORDINATES.STOP2,
                            CIVICVARIANTSGROUPSCOORDINATES.REPRESENTATIVETRANSCRIPT2,
                            CIVICVARIANTSGROUPSCOORDINATES.START,
                            CIVICVARIANTSGROUPSCOORDINATES.REPRESENTATIVETRANSCRIPT,
                            CIVICVARIANTSGROUPSCOORDINATES.ENSEMBLVERSION,
                            CIVICVARIANTSGROUPSCOORDINATES.CHROMOSOME,
                            CIVICVARIANTSGROUPSCOORDINATES.REFERENCEBUILD,
                            CIVICVARIANTSGROUPSCOORDINATES.CIVICVARIANTSGROUPSVARIANTSID)
                            .values(variants.coordinates().chromosome2(),
                                    variants.coordinates().referenceBases(),
                                    variants.coordinates().start2(),
                                    variants.coordinates().variantBases(),
                                    variants.coordinates().stop(),
                                    variants.coordinates().stop2(),
                                    variants.coordinates().representativeTranscript2(),
                                    variants.coordinates().start(),
                                    variants.coordinates().representativeTranscript(),
                                    variants.coordinates().ensemblVersion(),
                                    variants.coordinates().chromosome(),
                                    variants.coordinates().referenceBuild(),
                                    idVariantGroupVariants)
                            .execute();
                }

                for (CivicVariantTypes variantTypesGroup : variants.variant_types()) {
                    context.insertInto(CIVICVARIANTSGROUPSTYPES,
                            CIVICVARIANTSGROUPSTYPES.DISPLAYNAME,
                            CIVICVARIANTSGROUPSTYPES.DESCRIPTION,
                            CIVICVARIANTSGROUPSTYPES.URL,
                            CIVICVARIANTSGROUPSTYPES.SOID,
                            CIVICVARIANTSGROUPSTYPES.IDVARIANTTYPES,
                            CIVICVARIANTSGROUPSTYPES.NAME,
                            CIVICVARIANTSGROUPSTYPES.CIVICVARIANTSGROUPSVARIANTSID)
                            .values(variantTypesGroup.displayName(),
                                    variantTypesGroup.description(),
                                    variantTypesGroup.url(),
                                    variantTypesGroup.soId(),
                                    variantTypesGroup.id(),
                                    variantTypesGroup.name(),
                                    idVariantGroupVariants)
                            .execute();
                }
            }
        }

        for (CivicEvidenceItems evidenceItems : civic.evidenceItem()) {
            int idEvidenceItems = context.insertInto(CIVICEVIDENCEITEMS,
                    CIVICEVIDENCEITEMS.STATUS,
                    CIVICEVIDENCEITEMS.RATING,
                    CIVICEVIDENCEITEMS.DRUGINTERACTIONTYPE,
                    CIVICEVIDENCEITEMS.DESCRIPTION,
                    CIVICEVIDENCEITEMS.OPENCHANGECOUNT,
                    CIVICEVIDENCEITEMS.EVIDENCETYPE,
                    CIVICEVIDENCEITEMS.VARIANTORIGIN,
                    CIVICEVIDENCEITEMS.EVIDENCEDIRECTION,
                    CIVICEVIDENCEITEMS.VARIANTID,
                    CIVICEVIDENCEITEMS.CLINICALSIGNIFICANCE,
                    CIVICEVIDENCEITEMS.EVIDENCELEVEL,
                    CIVICEVIDENCEITEMS.TYPE,
                    CIVICEVIDENCEITEMS.IDEVIDENCEITEMS,
                    CIVICEVIDENCEITEMS.NAME,
                    CIVICEVIDENCEITEMS.CIVICID)
                    .values(evidenceItems.status(),
                            evidenceItems.rating(),
                            evidenceItems.drugInteractionType(),
                            evidenceItems.description(),
                            evidenceItems.openChangeCount(),
                            evidenceItems.evidenceType(),
                            evidenceItems.variantOrigin(),
                            evidenceItems.evidenceDirection(),
                            evidenceItems.variantId(),
                            evidenceItems.clinicalSignificance(),
                            evidenceItems.evidenceLevel(),
                            evidenceItems.type(),
                            evidenceItems.id(),
                            evidenceItems.name(),
                            id)
                    .returning(CIVICEVIDENCEITEMS.ID)
                    .fetchOne()
                    .getValue(CIVICEVIDENCEITEMS.ID);

            for (CivicDrugs drugs : evidenceItems.drugs()) {
                context.insertInto(CIVICDRUGS, CIVICDRUGS.PUBCHEMID, CIVICDRUGS.IDDRUGS, CIVICDRUGS.NAME, CIVICDRUGS.CIVICEVIDENCEITEMSID)
                        .values(drugs.pubchemId(), drugs.id(), drugs.name(), idEvidenceItems)
                        .execute();
            }
            context.insertInto(CIVICDISEASE,
                    CIVICDISEASE.DOID,
                    CIVICDISEASE.URL,
                    CIVICDISEASE.DISPLAYNAME,
                    CIVICDISEASE.IDDISEASE,
                    CIVICDISEASE.NAME,
                    CIVICDISEASE.CIVICEVIDENCEITEMSID)
                    .values(evidenceItems.disease().doid(),
                            evidenceItems.disease().url(),
                            evidenceItems.disease().displayName(),
                            evidenceItems.disease().id(),
                            evidenceItems.disease().name(),
                            idEvidenceItems)
                    .execute();

            int idEvidenceItemsSource = context.insertInto(CIVICEVIDENCEITEMSSOURCE,
                    CIVICEVIDENCEITEMSSOURCE.STATUS,
                    CIVICEVIDENCEITEMSSOURCE.OPENACCESS,
                    CIVICEVIDENCEITEMSSOURCE.NAME,
                    CIVICEVIDENCEITEMSSOURCE.JOURNAL,
                    CIVICEVIDENCEITEMSSOURCE.CITATION,
                    CIVICEVIDENCEITEMSSOURCE.PMC_ID,
                    CIVICEVIDENCEITEMSSOURCE.FULLJOURNALTITLE,
                    CIVICEVIDENCEITEMSSOURCE.SOURCEURL,
                    CIVICEVIDENCEITEMSSOURCE.PUBMEDID,
                    CIVICEVIDENCEITEMSSOURCE.ISREVIEW,
                    CIVICEVIDENCEITEMSSOURCE.IDSOURCE,
                    CIVICEVIDENCEITEMSSOURCE.CIVICEVIDENCEITEMSID)
                    .values(evidenceItems.source().status(),
                            evidenceItems.source().openAccess(),
                            evidenceItems.source().name(),
                            evidenceItems.source().journal(),
                            evidenceItems.source().citation(),
                            evidenceItems.source().pmc_Id(),
                            evidenceItems.source().fullJournalTitle(),
                            evidenceItems.source().sourceUrl(),
                            evidenceItems.source().pubmedId(),
                            evidenceItems.source().isReview(),
                            evidenceItems.source().id(),
                            idEvidenceItems)
                    .returning(CIVICEVIDENCEITEMSSOURCE.CIVICEVIDENCEITEMSID)
                    .fetchOne()
                    .getValue(CIVICEVIDENCEITEMSSOURCE.CIVICEVIDENCEITEMSID);

            context.insertInto(CIVICEVIDENCEITEMSPUBLICATION,
                    CIVICEVIDENCEITEMSPUBLICATION.YEAR,
                    CIVICEVIDENCEITEMSPUBLICATION.DAY,
                    CIVICEVIDENCEITEMSPUBLICATION.MONTH,
                    CIVICEVIDENCEITEMSPUBLICATION.CIVICEVIDENCEITEMSSOURCEID)
                    .values(evidenceItems.source().publicationDate().year(),
                            evidenceItems.source().publicationDate().day(),
                            evidenceItems.source().publicationDate().month(),
                            idEvidenceItemsSource)
                    .execute();

            for (CivicClinicalTrial clinicalTrial : evidenceItems.source().clinicalTrials()) {
                context.insertInto(CIVICEVIDENCEITEMSCLINICALTRIAL,
                        CIVICEVIDENCEITEMSCLINICALTRIAL.NCT_ID,
                        CIVICEVIDENCEITEMSCLINICALTRIAL.DESCRIPTION,
                        CIVICEVIDENCEITEMSCLINICALTRIAL.CLINICAL_TRIAL_URL,
                        CIVICEVIDENCEITEMSCLINICALTRIAL.NAME,
                        CIVICEVIDENCEITEMSCLINICALTRIAL.CIVICEVIDENCEITEMSSOURCEID)
                        .values(clinicalTrial.nct_id(),
                                clinicalTrial.description(),
                                clinicalTrial.clinical_trial_url(),
                                clinicalTrial.name(),
                                idEvidenceItemsSource)
                        .execute();
            }

        }

        if (civic.sources() != null) {
            for (CivicSource source : civic.sources()) {
                int idSource = context.insertInto(CIVICSOURCE,
                        CIVICSOURCE.STATUS,
                        CIVICSOURCE.OPENACCESS,
                        CIVICSOURCE.NAME,
                        CIVICSOURCE.JOURNAL,
                        CIVICSOURCE.CITATION,
                        CIVICSOURCE.PMC_ID,
                        CIVICSOURCE.FULLJOURNALTITLE,
                        CIVICSOURCE.SOURCEURL,
                        CIVICSOURCE.PUBMEDID,
                        CIVICSOURCE.ISREVIEW,
                        CIVICSOURCE.IDSOURCE,
                        CIVICSOURCE.CIVICID)
                        .values(source.status(),
                                source.openAccess(),
                                source.name(),
                                source.journal(),
                                source.citation(),
                                source.pmc_Id(),
                                source.fullJournalTitle(),
                                source.sourceUrl(),
                                source.pubmedId(),
                                source.isReview(),
                                source.id(),
                                id)
                        .returning(CIVICSOURCE.CIVICID)
                        .fetchOne()
                        .getValue(CIVICSOURCE.CIVICID);

                //TODO: fix issue with foreign key
//                context.insertInto(CIVICPUBLICATION,
//                        CIVICPUBLICATION.YEAR,
//                        CIVICPUBLICATION.DAY,
//                        CIVICPUBLICATION.MONTH,
//                        CIVICPUBLICATION.CIVICSOURCEID)
//                        .values(source.publicationDate().year(), source.publicationDate().day(), source.publicationDate().month(), idSource)
//                        .execute();

                for (CivicClinicalTrial clinicalTrial : source.clinicalTrials()) {
                    context.insertInto(CIVICCLINICALTRIAL,
                            CIVICCLINICALTRIAL.NCT_ID,
                            CIVICCLINICALTRIAL.DESCRIPTION,
                            CIVICCLINICALTRIAL.CLINICAL_TRIAL_URL,
                            CIVICCLINICALTRIAL.NAME,
                            CIVICCLINICALTRIAL.CIVICSOURCEID)
                            .values(clinicalTrial.nct_id(),
                                    clinicalTrial.description(),
                                    clinicalTrial.clinical_trial_url(),
                                    clinicalTrial.name(),
                                    idSource)
                            .execute();
                }
            }
        }


        context.insertInto(CIVICERROR, CIVICERROR.CIVICID).values(id).execute();

    }

    private void writeKbSpecificObject(int viccEntryId, @NotNull KbSpecificObject object) {
        if (object instanceof Sage) {
            importSageinSQL(viccEntryId, object);
        }
        if (object instanceof BRCA) {
            importBRCAinSQL(viccEntryId, object);
        }
        if (object instanceof Cgi) {
            importCGIinSQL(viccEntryId, object);
        }
        if (object instanceof Jax) {
            importJAXinSQL(viccEntryId, object);
        }
        if (object instanceof JaxTrials) {
            importJaxTrialsinSQL(viccEntryId, object);
        }
        if (object instanceof Pmkb) {
            importPmkbinSQL(viccEntryId, object);
        }
        if (object instanceof Oncokb) {
            importOncokbinSQL(viccEntryId, object);
        }
        if (object instanceof MolecularMatchTrials) {
            importMolecularMatchTrials(viccEntryId, object);
        }
        if (object instanceof Civic) {
            importCivic(viccEntryId, object);
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
                    TAXONOMY.ENVIRONMENTCONTEXTID)
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
            context.insertInto(APPROVEDCOUNTRY, APPROVEDCOUNTRY.APPROVEDCOUNTRYNAME, APPROVEDCOUNTRY.ENVIRONMENTCONTEXTID)
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
        context.deleteFrom(BRCAPART1).execute();
        context.deleteFrom(BRCAPART2).execute();
        context.deleteFrom(CGI).execute();
        context.deleteFrom(CGICDNA).execute();
        context.deleteFrom(CGIINDIVIDUALMUTATION).execute();
        context.deleteFrom(CGIGDNA).execute();
        context.deleteFrom(CGITRANSCRIPT).execute();
        context.deleteFrom(CGISTRAND).execute();
        context.deleteFrom(CGIINFO).execute();
        context.deleteFrom(CGIREGION).execute();
        context.deleteFrom(JAX).execute();
        context.deleteFrom(JAXMOLECULARPROFILE).execute();
        context.deleteFrom(JAXTHERAPY).execute();
        context.deleteFrom(JAXINDICATIONS).execute();
        context.deleteFrom(JAXREFERENCES).execute();
        context.deleteFrom(JAXTRIALS).execute();
        context.deleteFrom(JAXTRIALSINDICATIONS).execute();
        context.deleteFrom(JAXTRIALSVARIANTREQUIREMENTDETAILS).execute();
        context.deleteFrom(JAXTRIALSMOLECULARPROFILE).execute();
        context.deleteFrom(JAXTRIALSTHERAPIES).execute();
        context.deleteFrom(PMKB).execute();
        context.deleteFrom(PMKBTISSUE).execute();
        context.deleteFrom(PMKBTUMOR).execute();
        context.deleteFrom(PMKBVARIANT).execute();
        context.deleteFrom(PMKBGENE).execute();
        context.deleteFrom(ONCOKB).execute();
        context.deleteFrom(ONCOKBBIOLOGICAL).execute();
        context.deleteFrom(ONCOKBVARIANTBIOLOGICAL).execute();
        context.deleteFrom(ONCOKBCONSEQUENCESBIOLOGICAL).execute();
        context.deleteFrom(ONCOKBGENEBIOLOGICAL).execute();
        context.deleteFrom(ONCOKBGENEALIASESBIOLOGICAL).execute();
        context.deleteFrom(ONCOKBCLINICAL).execute();
        context.deleteFrom(ONCOKBDRUGABSTRACTSCLINICAL).execute();
        context.deleteFrom(ONCOKBVARIANTCLINICAL).execute();
        context.deleteFrom(ONCOKBCONSEQUENCESCLINICAL).execute();
        context.deleteFrom(ONCOKBGENECLINICAL).execute();
        context.deleteFrom(ONCOKBGENEALIASESCLINICAL).execute();
        context.deleteFrom(MOLECULARMATCHTRIALS).execute();
        context.deleteFrom(MOLECULARMATCHTRIALSALTERATIONS).execute();
        context.deleteFrom(MOLECULARMATCHTRIALSINTERVATIONS).execute();
        context.deleteFrom(MOLECULARMATCHTRIALSOTHERNAME).execute();
        context.deleteFrom(MOLECULARMATCHTRIALSOTHERGROUPLABEL).execute();
        context.deleteFrom(MOLECULARMATCHTRIALSOVERALLCONTACT).execute();
        context.deleteFrom(MOLECULARMATCHTRIALSTAGS).execute();
        context.deleteFrom(MOLECULARMATCHTRIALSLOCATIONS).execute();
        context.deleteFrom(MOLECULARMATCHTRIALSCONTACT).execute();
        context.deleteFrom(MOLECULARMATCHTRIALSGEO).execute();
        context.deleteFrom(MOLECULARMATCHTRIALSLOCATION).execute();
        context.deleteFrom(MOLECULARMATCHTRIALSCOORDINATES).execute();
        context.deleteFrom(CIVIC).execute();
        context.deleteFrom(CIVICASSERTIONS).execute();
        context.deleteFrom(CIVICHGVSEXPRESSIONS).execute();
        context.deleteFrom(CIVICCLINVARENTRIES).execute();
        context.deleteFrom(CIVICVARIANTALIASES).execute();
        context.deleteFrom(CIVICVARIANTTYPES).execute();
        context.deleteFrom(CIVICDESCRIPTION).execute();
        context.deleteFrom(CIVICCOORDINATES).execute();
        context.deleteFrom(CIVICVARIANTSGROUPS).execute();
        context.deleteFrom(CIVICVARIANTSGROUPSCOORDINATES).execute();
        context.deleteFrom(CIVICVARIANTSGROUPSTYPES).execute();
        context.deleteFrom(CIVICEVIDENCEITEMS).execute();
        context.deleteFrom(CIVICDRUGS).execute();
        context.deleteFrom(CIVICDISEASE).execute();
        context.deleteFrom(CIVICEVIDENCEITEMSSOURCE).execute();
        context.deleteFrom(CIVICEVIDENCEITEMSPUBLICATION).execute();
        context.deleteFrom(CIVICEVIDENCEITEMSCLINICALTRIAL).execute();
        context.deleteFrom(CIVICSOURCE).execute();
        context.deleteFrom(CIVICERROR).execute();
        context.deleteFrom(CIVICPUBLICATION).execute();
        context.deleteFrom(CIVICCLINICALTRIAL).execute();
        context.deleteFrom(CIVICLIFECYCLEACTIONS).execute();
        context.deleteFrom(CIVICLASTCOMMENTEDON).execute();
        context.deleteFrom(CIVICLASTCOMMENTEDONUSER).execute();
        context.deleteFrom(CIVICLASTCOMMENTEDONAVATARS).execute();
        context.deleteFrom(CIVICLASTCOMMENTEDONORGANIZATION).execute();
        context.deleteFrom(CIVICLASTCOMMENTEDONPROFILEIMAGE).execute();
        context.deleteFrom(CIVICLASTMODIFIED).execute();
        context.deleteFrom(CIVICLASTMODIFIEDUSER).execute();
        context.deleteFrom(CIVICLASTMODIFIEDAVATARS).execute();
        context.deleteFrom(CIVICLASTMODIFIEDORGANIZATION).execute();
        context.deleteFrom(CIVICLASTMODIFIEDPROFILEIMAGE).execute();
        context.deleteFrom(CIVICLASTREVIEWED).execute();
        context.deleteFrom(CIVICLASTREVIEWEDUSER).execute();
        context.deleteFrom(CIVICLASTREVIEWEDAVATARS).execute();
        context.deleteFrom(CIVICLASTREVIEWEDORGANIZATION).execute();
        context.deleteFrom(CIVICLASTREVIEWEDPROFILEIMAGE).execute();
    }
}
