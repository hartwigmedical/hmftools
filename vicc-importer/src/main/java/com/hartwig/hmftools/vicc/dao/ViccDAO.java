package com.hartwig.hmftools.vicc.dao;

import static com.hartwig.hmftools.vicc.database.Tables.APPROVEDCOUNTRY;
import static com.hartwig.hmftools.vicc.database.Tables.ASSOCIATION;
import static com.hartwig.hmftools.vicc.database.Tables.ASSOCIATIONVARIANT;
import static com.hartwig.hmftools.vicc.database.Tables.DEVTAG;
import static com.hartwig.hmftools.vicc.database.Tables.ENVIRONMENTALCONTEXT;
import static com.hartwig.hmftools.vicc.database.Tables.EVIDENCE;
import static com.hartwig.hmftools.vicc.database.Tables.EVIDENCEINFO;
import static com.hartwig.hmftools.vicc.database.Tables.EVIDENCETYPE;
import static com.hartwig.hmftools.vicc.database.Tables.FEATURE;
import static com.hartwig.hmftools.vicc.database.Tables.FEATUREATTRIBUTE;
import static com.hartwig.hmftools.vicc.database.Tables.FEATUREINFO;
import static com.hartwig.hmftools.vicc.database.Tables.FEATURENAME;
import static com.hartwig.hmftools.vicc.database.Tables.GENE;
import static com.hartwig.hmftools.vicc.database.Tables.GENEIDENTIFIER;
import static com.hartwig.hmftools.vicc.database.Tables.HIERARCHY;
import static com.hartwig.hmftools.vicc.database.Tables.LINK;
import static com.hartwig.hmftools.vicc.database.Tables.PHENOTYPE;
import static com.hartwig.hmftools.vicc.database.Tables.PHENOTYPETYPE;
import static com.hartwig.hmftools.vicc.database.Tables.PROVENANCE;
import static com.hartwig.hmftools.vicc.database.Tables.PUBLICATIONURL;
import static com.hartwig.hmftools.vicc.database.Tables.SEQUENCEONTOLOGY;
import static com.hartwig.hmftools.vicc.database.Tables.SYNONYM;
import static com.hartwig.hmftools.vicc.database.Tables.TAG;
import static com.hartwig.hmftools.vicc.database.Tables.TAXONOMY;
import static com.hartwig.hmftools.vicc.database.Tables.VICCENTRY;

import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.SQLException;
import java.util.List;

import com.hartwig.hmftools.vicc.datamodel.Association;
import com.hartwig.hmftools.vicc.datamodel.EnvironmentalContext;
import com.hartwig.hmftools.vicc.datamodel.Evidence;
import com.hartwig.hmftools.vicc.datamodel.EvidenceInfo;
import com.hartwig.hmftools.vicc.datamodel.EvidenceType;
import com.hartwig.hmftools.vicc.datamodel.Feature;
import com.hartwig.hmftools.vicc.datamodel.FeatureAttribute;
import com.hartwig.hmftools.vicc.datamodel.FeatureInfo;
import com.hartwig.hmftools.vicc.datamodel.GeneIdentifier;
import com.hartwig.hmftools.vicc.datamodel.KbSpecificObject;
import com.hartwig.hmftools.vicc.datamodel.Phenotype;
import com.hartwig.hmftools.vicc.datamodel.PhenotypeType;
import com.hartwig.hmftools.vicc.datamodel.SequenceOntology;
import com.hartwig.hmftools.vicc.datamodel.Taxonomy;
import com.hartwig.hmftools.vicc.datamodel.ViccEntry;
import com.hartwig.hmftools.vicc.datamodel.brca.Brca;
import com.hartwig.hmftools.vicc.datamodel.cgi.Cgi;
import com.hartwig.hmftools.vicc.datamodel.civic.Civic;
import com.hartwig.hmftools.vicc.datamodel.jax.Jax;
import com.hartwig.hmftools.vicc.datamodel.jaxtrials.JaxTrials;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.MolecularMatch;
import com.hartwig.hmftools.vicc.datamodel.molecularmatchtrials.MolecularMatchTrials;
import com.hartwig.hmftools.vicc.datamodel.oncokb.OncoKb;
import com.hartwig.hmftools.vicc.datamodel.pmkb.Pmkb;
import com.hartwig.hmftools.vicc.datamodel.sage.Sage;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;
import org.jooq.DSLContext;
import org.jooq.SQLDialect;
import org.jooq.conf.MappedSchema;
import org.jooq.conf.RenderMapping;
import org.jooq.conf.Settings;
import org.jooq.impl.DSL;

public class ViccDAO {

    private static final Logger LOGGER = LogManager.getLogger(ViccDAO.class);

    private static final String DEV_CATALOG = "vicc_test";

    @NotNull
    private final DSLContext context;

    @NotNull
    public static ViccDAO connectToViccDAO(@NotNull String userName, @NotNull String password, @NotNull String url) throws SQLException {
        System.setProperty("org.jooq.no-logo", "true");
        System.setProperty("org.jooq.no-tips", "true");

        Connection conn = DriverManager.getConnection(url, userName, password);
        String catalog = conn.getCatalog();
        LOGGER.info("Connecting to database '{}'", catalog);

        return new ViccDAO(DSL.using(conn, SQLDialect.MYSQL, settings(catalog)));
    }

    @Nullable
    private static Settings settings(@NotNull String catalog) {
        if (catalog.equals(DEV_CATALOG)) {
            return null;
        }

        return new Settings().withRenderMapping(new RenderMapping().withSchemata(new MappedSchema().withInput(DEV_CATALOG)
                .withOutput(catalog)));
    }

    private ViccDAO(@NotNull final DSLContext context) {
        this.context = context;
    }

    public void writeViccEntry(@NotNull ViccEntry viccEntry) {
        int id = context.insertInto(VICCENTRY, VICCENTRY.SOURCE)
                .values(viccEntry.source().display())
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
        writeKbSpecificObject(id, viccEntry.kbSpecificObject());
    }

    public void deleteAll() {
        // Note: The order should be "from branch to root" to avoid constraint violation.
        BRCADAOFunctions.deleteAll(context);
        CgiDAOFunctions.deleteAll(context);
        CivicDAOFunctions.deleteAll(context);
        JaxDAOFunctions.deleteAll(context);
        JaxTrialsDAOFunctions.deleteAll(context);
        MolecularMatchDAOFunctions.deleteAll(context);
        MolecularMatchTrialsDAOFunctions.deleteAll(context);
        OncokbDAOFunctions.deleteAll(context);
        PmkbDAOFunctions.deleteAll(context);
        SageDAOFunctions.deleteAll(context);

        // Below tables are part of Association
        context.deleteFrom(ASSOCIATIONVARIANT).execute();
        context.deleteFrom(EVIDENCETYPE).execute();
        context.deleteFrom(EVIDENCEINFO).execute();
        context.deleteFrom(EVIDENCE).execute();
        context.deleteFrom(PUBLICATIONURL).execute();
        context.deleteFrom(PHENOTYPETYPE).execute();
        context.deleteFrom(PHENOTYPE).execute();
        context.deleteFrom(APPROVEDCOUNTRY).execute();
        context.deleteFrom(TAXONOMY).execute();
        context.deleteFrom(ENVIRONMENTALCONTEXT).execute();

        // Below tables are part of Feature
        context.deleteFrom(FEATUREINFO).execute();
        context.deleteFrom(FEATUREATTRIBUTE).execute();
        context.deleteFrom(PROVENANCE).execute();
        context.deleteFrom(SYNONYM).execute();
        context.deleteFrom(LINK).execute();
        context.deleteFrom(HIERARCHY).execute();
        context.deleteFrom(SEQUENCEONTOLOGY).execute();

        // Below tables are part of VICC Entry
        context.deleteFrom(ASSOCIATION).execute();
        context.deleteFrom(FEATURENAME).execute();
        context.deleteFrom(FEATURE).execute();
        context.deleteFrom(GENE).execute();
        context.deleteFrom(GENEIDENTIFIER).execute();
        context.deleteFrom(DEVTAG).execute();
        context.deleteFrom(TAG).execute();

        context.deleteFrom(VICCENTRY).execute();
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

    private void writeFeatureNames(int viccEntryId, @NotNull List<String> featureNames) {
        for (String featureName : featureNames) {
            context.insertInto(FEATURENAME, FEATURENAME.NAMEOFFEATURE, FEATURENAME.VICCENTRYID).values(featureName, viccEntryId).execute();
        }
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
            writeFeatureInfo(id, feature.info());
            writeFeatureAttribute(id, feature.attribute());
            writeProvenance(id, feature.provenance());
            writeSynonyms(id, feature.synonyms());
            writeLinks(id, feature.links());
            writeSequenceOntology(id, feature.sequenceOntology());
        }
    }

    private void writeFeatureInfo(int featureId, @Nullable FeatureInfo featureInfo) {
        if (featureInfo != null) {
            context.insertInto(FEATUREINFO, FEATUREINFO.GERMLINEORSOMATIC, FEATUREINFO.FEATUREID)
                    .values(featureInfo.germlineOrSomatic(), featureId)
                    .execute();
        }
    }

    private void writeFeatureAttribute(int featureId, @Nullable FeatureAttribute featureAttribute) {
        if (featureAttribute != null) {
            context.insertInto(FEATUREATTRIBUTE,
                    FEATUREATTRIBUTE.AMINOACIDCHANGE,
                    FEATUREATTRIBUTE.GERMLINE,
                    FEATUREATTRIBUTE.PARTNERGENE,
                    FEATUREATTRIBUTE.DESCRIPTION,
                    FEATUREATTRIBUTE.EXONS,
                    FEATUREATTRIBUTE.NOTES,
                    FEATUREATTRIBUTE.COSMIC,
                    FEATUREATTRIBUTE.EFFECT,
                    FEATUREATTRIBUTE.CNVTYPE,
                    FEATUREATTRIBUTE.FEATUREATTRIBUTEID,
                    FEATUREATTRIBUTE.CYTOBAND,
                    FEATUREATTRIBUTE.VARIANTTYPE,
                    FEATUREATTRIBUTE.DNACHANGE,
                    FEATUREATTRIBUTE.CODONS,
                    FEATUREATTRIBUTE.CHROMOSOMEBASEDCNV,
                    FEATUREATTRIBUTE.TRANSCRIPT,
                    FEATUREATTRIBUTE.DESCRIPTIONTYPE,
                    FEATUREATTRIBUTE.CHROMOSOME,
                    FEATUREINFO.FEATUREID)
                    .values(featureAttribute.aminoAcidChange(),
                            featureAttribute.germline(),
                            featureAttribute.partnerGene(),
                            featureAttribute.description(),
                            featureAttribute.exons(),
                            featureAttribute.notes(),
                            featureAttribute.cosmic(),
                            featureAttribute.effect(),
                            featureAttribute.cnvType(),
                            featureAttribute.id(),
                            featureAttribute.cytoband(),
                            featureAttribute.variantType(),
                            featureAttribute.dnaChange(),
                            featureAttribute.codons(),
                            featureAttribute.chromosomeBasedCnv(),
                            featureAttribute.transcript(),
                            featureAttribute.descriptionType(),
                            featureAttribute.chromosome(),
                            featureId)
                    .execute();
        }
    }

    private void writeProvenance(int featureId, @NotNull List<String> provenances) {
        for (String provenance : provenances) {
            context.insertInto(PROVENANCE, PROVENANCE.PROVENANCENAME, PROVENANCE.FEATUREID).values(provenance, featureId).execute();
        }
    }

    private void writeSynonyms(int featureId, @NotNull List<String> synonyms) {
        for (String synonym : synonyms) {
            context.insertInto(SYNONYM, SYNONYM.SYNONYMNAME, SYNONYM.FEATUREID).values(synonym, featureId).execute();
        }
    }

    private void writeLinks(int featureId, @NotNull List<String> links) {
        for (String link : links) {
            context.insertInto(LINK, LINK.LINKNAME, LINK.FEATUREID).values(link, featureId).execute();
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

    private void writeHierarchy(int sequenceOntologyId, @NotNull List<String> hierarchies) {
        for (String hierarchy : hierarchies) {
            context.insertInto(HIERARCHY, HIERARCHY.HIERARCHYNAME, HIERARCHY.SEQUENCEONTOLOGYID)
                    .values(hierarchy, sequenceOntologyId)
                    .execute();
        }
    }

    private void writeAssociation(int viccEntryId, @NotNull Association association) {
        int id = context.insertInto(ASSOCIATION,
                ASSOCIATION.EVIDENCELEVEL,
                ASSOCIATION.EVIDENCELABEL,
                ASSOCIATION.RESPONSETYPE,
                ASSOCIATION.DRUGLABELS,
                ASSOCIATION.SOURCELINK,
                ASSOCIATION.DESCRIPTION,
                ASSOCIATION.ONCOGENIC,
                ASSOCIATION.VICCENTRYID)
                .values(association.evidenceLevel(),
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
        writeVariantNames(id, association.variantNames());
        writeEvidence(id, association.evidence());
        writePublicationsUrls(id, association.publicationUrls());
        writePhenotype(id, association.phenotype());
        writeEnvironmentalContexts(id, association.environmentalContexts());
    }

    private void writeVariantNames(int associationId, @NotNull List<String> variantNames) {
        for (String variant : variantNames) {
            context.insertInto(ASSOCIATIONVARIANT, ASSOCIATIONVARIANT.VARIANTNAME, ASSOCIATIONVARIANT.ASSOCIATIONID)
                    .values(variant, associationId)
                    .execute();
        }
    }

    private void writeEvidence(int associationId, @NotNull Evidence evidence) {
        int id = context.insertInto(EVIDENCE, EVIDENCE.DESCRIPTION, EVIDENCE.ASSOCIATIONID)
                .values(evidence.description(), associationId)
                .returning(EVIDENCE.ID)
                .fetchOne()
                .getValue(EVIDENCE.ID);
        writeEvidenceInfo(id, evidence.info());
        writeEvidenceType(id, evidence.evidenceType());
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

    private void writePublicationsUrls(int associationId, @NotNull List<String> publicationsUrls) {
        for (String publicationUrl : publicationsUrls) {
            context.insertInto(PUBLICATIONURL, PUBLICATIONURL.URLOFPUBLICATION, PUBLICATIONURL.ASSOCIATIONID)
                    .values(publicationUrl, associationId)
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

    private void writeEnvironmentalContexts(int associationId, @Nullable List<EnvironmentalContext> environmentalContexts) {
        if (environmentalContexts != null) {
            for (EnvironmentalContext environmentalContext : environmentalContexts) {
                int id = context.insertInto(ENVIRONMENTALCONTEXT,
                        ENVIRONMENTALCONTEXT.TERM,
                        ENVIRONMENTALCONTEXT.DESCRIPTION,
                        ENVIRONMENTALCONTEXT.SOURCE,
                        ENVIRONMENTALCONTEXT.USANSTEM,
                        ENVIRONMENTALCONTEXT.TOXICITY,
                        ENVIRONMENTALCONTEXT.IDENVIRONMENTALCONTEXT,
                        ENVIRONMENTALCONTEXT.ASSOCIATIONID)
                        .values(environmentalContext.term(),
                                environmentalContext.description(),
                                environmentalContext.source(),
                                environmentalContext.usanStem(),
                                environmentalContext.toxicity(),
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

    private void writeApprovedCountries(int environmentalContextsId, @NotNull List<String> approvedCountries) {
        for (String approvesCountry : approvedCountries) {
            context.insertInto(APPROVEDCOUNTRY, APPROVEDCOUNTRY.APPROVEDCOUNTRYNAME, APPROVEDCOUNTRY.ENVIRONMENTALCONTEXTID)
                    .values(approvesCountry, environmentalContextsId)
                    .execute();
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
                    TAXONOMY.ENVIRONMENTALCONTEXTID)
                    .values(taxonomy.kingdom(),
                            taxonomy.directParent(),
                            taxonomy.classs(),
                            taxonomy.subClass(),
                            taxonomy.superClass(),
                            environmentalContextsId)
                    .execute();
        }
    }

    private void writeKbSpecificObject(int viccEntryId, @NotNull KbSpecificObject object) {
        if (object instanceof Sage) {
            SageDAOFunctions.write(context, viccEntryId, (Sage) object);
        } else if (object instanceof Brca) {
            BRCADAOFunctions.write(context, viccEntryId, (Brca) object);
        } else if (object instanceof Cgi) {
            CgiDAOFunctions.write(context, viccEntryId, (Cgi) object);
        } else if (object instanceof Jax) {
            JaxDAOFunctions.write(context, viccEntryId, (Jax) object);
        } else if (object instanceof JaxTrials) {
            JaxTrialsDAOFunctions.write(context, viccEntryId, (JaxTrials) object);
        } else if (object instanceof Pmkb) {
            PmkbDAOFunctions.write(context, viccEntryId, (Pmkb) object);
        } else if (object instanceof OncoKb) {
            OncokbDAOFunctions.write(context, viccEntryId, (OncoKb) object);
        } else if (object instanceof MolecularMatchTrials) {
            MolecularMatchTrialsDAOFunctions.write(context, viccEntryId, (MolecularMatchTrials) object);
        } else if (object instanceof Civic) {
            CivicDAOFunctions.write(context, viccEntryId, (Civic) object);
        } else if (object instanceof MolecularMatch) {
            MolecularMatchDAOFunctions.write(context, viccEntryId, (MolecularMatch) object);
        } else {
            LOGGER.warn("Could not resolve kb specific object of class {}", object.getClass().getCanonicalName());
        }
    }
}
