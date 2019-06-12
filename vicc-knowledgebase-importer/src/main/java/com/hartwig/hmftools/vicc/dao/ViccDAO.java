package com.hartwig.hmftools.vicc.dao;

import static com.hartwig.hmftools.vicc.database.Tables.ASSOCIATION;
import static com.hartwig.hmftools.vicc.database.Tables.DEVTAG;
import static com.hartwig.hmftools.vicc.database.Tables.ENVIRONMENTALCONTEXT;
import static com.hartwig.hmftools.vicc.database.Tables.EVIDENCE;
import static com.hartwig.hmftools.vicc.database.Tables.FEATURE;
import static com.hartwig.hmftools.vicc.database.Tables.FEATURENAME;
import static com.hartwig.hmftools.vicc.database.Tables.GENE;
import static com.hartwig.hmftools.vicc.database.Tables.GENEIDENTIFIER;
import static com.hartwig.hmftools.vicc.database.Tables.PHENOTYPE;
import static com.hartwig.hmftools.vicc.database.Tables.SEQUENCEONTOLOGY;
import static com.hartwig.hmftools.vicc.database.Tables.TAG;
import static com.hartwig.hmftools.vicc.database.Tables.VICCENTRY;

import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.SQLException;
import java.util.List;

import com.hartwig.hmftools.vicc.ViccJsonToSQLImporter;
import com.hartwig.hmftools.vicc.datamodel.Association;
import com.hartwig.hmftools.vicc.datamodel.EnvironmentalContext;
import com.hartwig.hmftools.vicc.datamodel.Evidence;
import com.hartwig.hmftools.vicc.datamodel.Feature;
import com.hartwig.hmftools.vicc.datamodel.GeneIdentifier;
import com.hartwig.hmftools.vicc.datamodel.SequenceOntology;
import com.hartwig.hmftools.vicc.datamodel.ViccEntry;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.util.Strings;
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
        writeFeature(id, viccEntry.features());
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
                ASSOCIATION.PUBLICATIONURLS,
                ASSOCIATION.DESCRIPTION,
                ASSOCIATION.ONCOGENIC,
                ASSOCIATION.VICCENTRYID)
                .values(association.variantName(),
                        association.evidenceLevel(),
                        association.evidenceLabel(),
                        association.responseType(),
                        association.drugLabels(),
                        association.sourceLink(),
                        "",
                        association.description(),
                        association.oncogenic(),
                        viccEntryId)
                .returning(ASSOCIATION.IDASSOCIATION)
                .fetchOne()
                .getValue(ASSOCIATION.IDASSOCIATION);
        writeEvidence(id, association);
        writePhenotype(id, association);
        writeEnvironmentalContext(id, association);
    }

    private void writeEvidence(int associationId, @Nullable Association association) {
        for (Evidence evidence : association.evidence()) {
            context.insertInto(EVIDENCE,
                    EVIDENCE.EVIDENCEINFOPUBLICATIONS,
                    EVIDENCE.EVIDENCETYPESOURCENAME,
                    EVIDENCE.EVIDENCETYPEID,
                    EVIDENCE.DESCRIPTION,
                    EVIDENCE.ASSOCIATIONENTRYID)
                    .values("", evidence.evidenceType().sourceName(), evidence.evidenceType().id(), evidence.description(), associationId)
                    .execute();
        }

    }

    private void writePhenotype(int associationId, @Nullable Association association) {
        context.insertInto(PHENOTYPE,
                PHENOTYPE.PHENOTYPETYPESOURCE,
                PHENOTYPE.PHENOTYPETYPETERM,
                PHENOTYPE.PHENOTYPETYPEID,
                PHENOTYPE.DESCRIPTION,
                PHENOTYPE.FAMILY,
                PHENOTYPE.ID,
                PHENOTYPE.ASSOCIATIONENTRYID)
                .values(association.phenotype().type().source(),
                        association.phenotype().type().term(),
                        association.phenotype().type().id(),
                        association.phenotype().description(),
                        association.phenotype().family(),
                        association.phenotype().id(),
                        associationId)
                .execute();
    }

    private void writeEnvironmentalContext(int associationId, @Nullable Association association) {
        for (EnvironmentalContext environmentalContext : association.environmentalContexts()) {
            context.insertInto(ENVIRONMENTALCONTEXT,
                    ENVIRONMENTALCONTEXT.TERM,
                    ENVIRONMENTALCONTEXT.DESCRIPTION,
                    ENVIRONMENTALCONTEXT.TAXONOMYKINGDOM,
                    ENVIRONMENTALCONTEXT.TAXONOMYDIRECTPARENT,
                    ENVIRONMENTALCONTEXT.TAXONOMYCLASS,
                    ENVIRONMENTALCONTEXT.TAXONOMYSUBCLASS,
                    ENVIRONMENTALCONTEXT.TAXONOMYSUPERCLASS,
                    ENVIRONMENTALCONTEXT.SOURCE,
                    ENVIRONMENTALCONTEXT.USANSTEM,
                    ENVIRONMENTALCONTEXT.APPROVEDCOUNTRIES,
                    ENVIRONMENTALCONTEXT.ID,
                    ENVIRONMENTALCONTEXT.ASSOCIATIONENTRYID)
                    .values(environmentalContext.term(),
                            environmentalContext.description(),
                            environmentalContext.taxonomy().kingdom(),
                            environmentalContext.taxonomy().directParent(),
                            environmentalContext.taxonomy().classs(),
                            environmentalContext.taxonomy().subClass(),
                            environmentalContext.taxonomy().superClass(),
                            environmentalContext.source(),
                            environmentalContext.usanStem(),
                            "",
                            environmentalContext.id(),
                            associationId)
                    .execute();
        }
    }

    private void writeSequenceOntology(int featureId, @Nullable SequenceOntology sequenceOntologies) {
        String hierarchy = Strings.EMPTY;
        context.insertInto(SEQUENCEONTOLOGY,
                SEQUENCEONTOLOGY.HIERARCHY,
                SEQUENCEONTOLOGY.SOID,
                SEQUENCEONTOLOGY.PARENTSOID,
                SEQUENCEONTOLOGY.NAME,
                SEQUENCEONTOLOGY.PARENTNAME,
                SEQUENCEONTOLOGY.FEATUREENTRYID).values(hierarchy, "", "", "", "", featureId).execute();
    }

    private void writeFeature(int viccEntryId, @NotNull List<Feature> features) {
        int idFeature = 1;
        for (Feature feature : features) {
            String linksString = Strings.EMPTY;
            String provenanceString = Strings.EMPTY;
            String synonymsString = Strings.EMPTY;

            context.insertInto(FEATURE,
                    FEATURE.NAME,
                    FEATURE.BIOMARKERTYPE,
                    FEATURE.REFERENCENAME,
                    FEATURE.CHROMOSOME,
                    FEATURE.START,
                    FEATURE.END,
                    FEATURE.REF,
                    FEATURE.ALT,
                    FEATURE.PROVENANCE,
                    FEATURE.PROVENANCERULE,
                    FEATURE.GENESYMBOL,
                    FEATURE.SYNONYMS,
                    FEATURE.ENTREZID,
                    FEATURE.LINKS,
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
                            provenanceString,
                            feature.provenanceRule(),
                            feature.geneSymbol(),
                            synonymsString,
                            feature.entrezId(),
                            linksString,
                            feature.description(),
                            viccEntryId)
                    .execute();
            writeSequenceOntology(idFeature, feature.sequenceOntology());
            idFeature++;
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
        context.deleteFrom(SEQUENCEONTOLOGY).execute();
        context.deleteFrom(ASSOCIATION).execute();
        context.deleteFrom(EVIDENCE).execute();
        context.deleteFrom(PHENOTYPE).execute();
        context.deleteFrom(ENVIRONMENTALCONTEXT).execute();
        context.deleteFrom(VICCENTRY).execute();
    }
}
