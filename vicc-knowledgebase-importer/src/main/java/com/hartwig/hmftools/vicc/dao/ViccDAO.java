package com.hartwig.hmftools.vicc.dao;

import static com.hartwig.hmftools.vicc.database.Tables.APPROVEDCOUNTRIES;
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
import static com.hartwig.hmftools.vicc.database.Tables.LINKS;
import static com.hartwig.hmftools.vicc.database.Tables.PHENOTYPE;
import static com.hartwig.hmftools.vicc.database.Tables.PHENOTYPETYPE;
import static com.hartwig.hmftools.vicc.database.Tables.PROVENANCE;
import static com.hartwig.hmftools.vicc.database.Tables.PUBLICATIONURL;
import static com.hartwig.hmftools.vicc.database.Tables.SEQUENCEONTOLOGY;
import static com.hartwig.hmftools.vicc.database.Tables.SYNONYMS;
import static com.hartwig.hmftools.vicc.database.Tables.TAG;
import static com.hartwig.hmftools.vicc.database.Tables.TAXONOMY;
import static com.hartwig.hmftools.vicc.database.Tables.VICCENTRY;

import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.SQLException;
import java.util.List;

import com.hartwig.hmftools.vicc.ViccJsonToSQLImporter;
import com.hartwig.hmftools.vicc.datamodel.Association;
import com.hartwig.hmftools.vicc.datamodel.EnvironmentalContext;
import com.hartwig.hmftools.vicc.datamodel.Evidence;
import com.hartwig.hmftools.vicc.datamodel.EvidenceInfo;
import com.hartwig.hmftools.vicc.datamodel.EvidenceType;
import com.hartwig.hmftools.vicc.datamodel.Feature;
import com.hartwig.hmftools.vicc.datamodel.GeneIdentifier;
import com.hartwig.hmftools.vicc.datamodel.Phenotype;
import com.hartwig.hmftools.vicc.datamodel.PhenotypeType;
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

        writeTags(id, viccEntry.tags());
        writeDevTags(id, viccEntry.devTags());
        writeGeneIdentifiers(id, viccEntry.geneIdentifiers());
        writeGenes(id, viccEntry.genes());
        writeFeatureNames(id, viccEntry.featureNames());
        writeFeatures(id, viccEntry.features());
        writeAssociation(id, viccEntry.association());
    }

    private void writeAssociation(int viccEntryId, @NotNull Association association) {
        int id = context.insertInto(ASSOCIATION,
                ASSOCIATION.VARIANTNAME,
                ASSOCIATION.EVIDENCELEVEL,
                ASSOCIATION.EVIDENCELABEL,
                ASSOCIATION.RESPONSETYPE,
                ASSOCIATION.DRUGLABELS,
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
        writePublicaionsUrls(id, association.publicationUrls());
        writePhenotype(id, association.phenotype());
        writeEnvironmentalContexts(id, association.environmentalContexts());
    }

    private void writeEnvironmentalContexts(int associationEntryId, @NotNull List<EnvironmentalContext> environmentalContexts) {
        for (EnvironmentalContext environmentalContext : environmentalContexts) {
            int id = context.insertInto(ENVIRONMENTALCONTEXT,
                    ENVIRONMENTALCONTEXT.TERM,
                    ENVIRONMENTALCONTEXT.DESCRIPTION,
                    ENVIRONMENTALCONTEXT.SOURCE,
                    ENVIRONMENTALCONTEXT.USANSTEM,
                    ENVIRONMENTALCONTEXT.IDENVIRONMENTALCONTEXTS,
                    ENVIRONMENTALCONTEXT.ASSOCIATIONENTRYID)
                    .values(environmentalContext.term(),
                            environmentalContext.description(),
                            environmentalContext.source(),
                            environmentalContext.usanStem(),
                            environmentalContext.id(),
                            associationEntryId)
                    .returning(ENVIRONMENTALCONTEXT.ID)
                    .fetchOne()
                    .getValue(ENVIRONMENTALCONTEXT.ID);
            writeApprovedCountries(id, environmentalContext.approvedCountries());
            writeTaxonomy(id, environmentalContext.taxonomy());
        }
    }

    private void writeTaxonomy(int environmentalContextsEntryId, @NotNull Taxonomy taxonomy) {
        context.insertInto(TAXONOMY,
                TAXONOMY.KINGDOM,
                TAXONOMY.DIRECTPARENT,
                TAXONOMY.CLASS,
                TAXONOMY.SUBCLASS,
                TAXONOMY.SUPERCLASS,
                TAXONOMY.ENVIRONMENTALCONTEXTSENTRYID)
                .values(taxonomy.kingdom(),
                        taxonomy.directParent(),
                        taxonomy.classs(),
                        taxonomy.subClass(),
                        taxonomy.superClass(),
                        environmentalContextsEntryId)
                .execute();
    }

    private void writeApprovedCountries(int environmentalContextsEntryId, @NotNull List<String> approvedCountries) {
        for (String approvesCountry : approvedCountries) {
            context.insertInto(APPROVEDCOUNTRIES, APPROVEDCOUNTRIES.APPROVEDCOUNTRIES_, APPROVEDCOUNTRIES.ENVIRONMENTALCONTEXTSENTRYID)
                    .values(approvesCountry, environmentalContextsEntryId)
                    .execute();
        }
    }

    private void writePhenotype(int associationEntryId, @NotNull Phenotype phenotype) {
        int id = context.insertInto(PHENOTYPE, PHENOTYPE.DESCRIPTION, PHENOTYPE.FAMILY, PHENOTYPE.IDPHENOTYPE, PHENOTYPE.ASSOCIATIONENTRYID)
                .values(phenotype.description(), phenotype.family(), phenotype.id(), associationEntryId)
                .returning(PHENOTYPE.ID)
                .fetchOne()
                .getValue(PHENOTYPE.ID);
        writePhenotypeType(id, phenotype.type());
    }

    private void writePhenotypeType(int phenoTypeEntryId, @NotNull PhenotypeType phenotypeType) {
        context.insertInto(PHENOTYPETYPE,
                PHENOTYPETYPE.SOURCE,
                PHENOTYPETYPE.TERM,
                PHENOTYPETYPE.IDPHENOTYPETYPE,
                PHENOTYPETYPE.PHENOTYPEENTRYID)
                .values(phenotypeType.source(), phenotypeType.term(), phenotypeType.id(), phenoTypeEntryId)
                .execute();
    }

    private void writePublicaionsUrls(int associationEntryId, @NotNull List<String> publicationsUrls) {
        for (String publicationUrl : publicationsUrls) {
            context.insertInto(PUBLICATIONURL, PUBLICATIONURL.PUBLICATIONURLS, PUBLICATIONURL.ASSOCIATIONENTRYID)
                    .values(publicationUrl, associationEntryId)
                    .execute();
        }
    }

    private void writeEvidence(int associationEntryId, @NotNull List<Evidence> evidences) {
        for (Evidence evidence : evidences) {
            int id = context.insertInto(EVIDENCE, EVIDENCE.DESCRIPTION, EVIDENCE.ASSOCIATIONENTRYID)
                    .values(evidence.description(), associationEntryId)
                    .returning(EVIDENCE.ID)
                    .fetchOne()
                    .getValue(EVIDENCE.ID);
            writeEvidenceInfo(id, evidence.info());
            writeEvidenceType(id, evidence.evidenceType());
        }
    }

    private void writeEvidenceInfo(int evidenceEntryId, @NotNull EvidenceInfo evidenceInfo) {
        for (String publication : evidenceInfo.publications()) {
            context.insertInto(EVIDENCEINFO, EVIDENCEINFO.PUBLICATIONS, EVIDENCEINFO.EVIDENCEENTRYID)
                    .values(publication, evidenceEntryId)
                    .execute();
        }
    }

    private void writeEvidenceType(int evidenceEntryId, @NotNull EvidenceType evidenceType) {
        context.insertInto(EVIDENCETYPE, EVIDENCETYPE.SOURCENAME, EVIDENCETYPE.IDEVIDENCETYPE, EVIDENCETYPE.EVIDENCEENTRYID)
                .values(evidenceType.sourceName(), evidenceType.id(), evidenceEntryId)
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

    private void writeProvenance(int featureEntryId, @NotNull List<String> provenances) {
        for (String provenance : provenances) {
            context.insertInto(PROVENANCE, PROVENANCE.PROVENANCE_, PROVENANCE.FEATUREENTRYID).values(provenance, featureEntryId).execute();
        }
    }

    private void writeSynonyms(int featureEntryId, @NotNull List<String> synonyms) {
        for (String synonym : synonyms) {
            context.insertInto(SYNONYMS, SYNONYMS.SYNONYMS_, SYNONYMS.FEATUREENTRYID).values(synonym, featureEntryId).execute();
        }
    }

    private void writeLinks(int featureEntryId, @NotNull List<String> links) {
        for (String link : links) {
            context.insertInto(LINKS, LINKS.LINKS_, LINKS.FEATUREENTRYID).values(link, featureEntryId).execute();
        }
    }

    private void writeSequenceOntology(int featureId, @Nullable SequenceOntology sequenceOntologies) {
        int id = context.insertInto(SEQUENCEONTOLOGY,
                SEQUENCEONTOLOGY.SOID,
                SEQUENCEONTOLOGY.PARENTSOID,
                SEQUENCEONTOLOGY.NAME,
                SEQUENCEONTOLOGY.PARENTNAME,
                SEQUENCEONTOLOGY.FEATUREENTRYID)
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

    private void writeHierarchy(int sequenceOntologyId, @Nullable List<String> hierchies) {
        for (String hierchy : hierchies) {
            context.insertInto(HIERARCHY, HIERARCHY.HIERARCHY_, HIERARCHY.SEQUENCEONTOLOGYENTRYID)
                    .values(hierchy, sequenceOntologyId)
                    .execute();
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
                context.insertInto(FEATURENAME, FEATURENAME.FEATURENAME_, FEATURENAME.VICCENTRYID)
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
        context.deleteFrom(SYNONYMS).execute();
        context.deleteFrom(LINKS).execute();
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
        context.deleteFrom(APPROVEDCOUNTRIES).execute();
        context.deleteFrom(TAXONOMY).execute();
        context.deleteFrom(VICCENTRY).execute();
    }
}
