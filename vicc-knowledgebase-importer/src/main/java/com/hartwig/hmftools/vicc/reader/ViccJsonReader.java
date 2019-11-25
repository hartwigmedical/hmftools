package com.hartwig.hmftools.vicc.reader;

import static com.hartwig.hmftools.vicc.reader.JsonFunctions.jsonArrayToStringList;
import static com.hartwig.hmftools.vicc.reader.JsonFunctions.nullableString;
import static com.hartwig.hmftools.vicc.reader.JsonFunctions.optionalNullableString;
import static com.hartwig.hmftools.vicc.reader.JsonFunctions.optionalString;
import static com.hartwig.hmftools.vicc.reader.JsonFunctions.optionalStringList;

import java.io.FileReader;
import java.io.IOException;
import java.util.List;

import javax.annotation.Nullable;

import com.google.common.collect.Lists;
import com.google.gson.JsonArray;
import com.google.gson.JsonElement;
import com.google.gson.JsonObject;
import com.google.gson.JsonParser;
import com.google.gson.stream.JsonReader;
import com.google.gson.stream.JsonToken;
import com.hartwig.hmftools.vicc.datamodel.Association;
import com.hartwig.hmftools.vicc.datamodel.EnvironmentalContext;
import com.hartwig.hmftools.vicc.datamodel.Evidence;
import com.hartwig.hmftools.vicc.datamodel.EvidenceInfo;
import com.hartwig.hmftools.vicc.datamodel.EvidenceType;
import com.hartwig.hmftools.vicc.datamodel.Feature;
import com.hartwig.hmftools.vicc.datamodel.FeatureAttribute;
import com.hartwig.hmftools.vicc.datamodel.FeatureInfo;
import com.hartwig.hmftools.vicc.datamodel.GeneIdentifier;
import com.hartwig.hmftools.vicc.datamodel.ImmutableAssociation;
import com.hartwig.hmftools.vicc.datamodel.ImmutableEnvironmentalContext;
import com.hartwig.hmftools.vicc.datamodel.ImmutableEvidence;
import com.hartwig.hmftools.vicc.datamodel.ImmutableEvidenceInfo;
import com.hartwig.hmftools.vicc.datamodel.ImmutableEvidenceType;
import com.hartwig.hmftools.vicc.datamodel.ImmutableFeature;
import com.hartwig.hmftools.vicc.datamodel.ImmutableFeatureAttribute;
import com.hartwig.hmftools.vicc.datamodel.ImmutableFeatureInfo;
import com.hartwig.hmftools.vicc.datamodel.ImmutableGeneIdentifier;
import com.hartwig.hmftools.vicc.datamodel.ImmutablePhenotype;
import com.hartwig.hmftools.vicc.datamodel.ImmutablePhenotypeType;
import com.hartwig.hmftools.vicc.datamodel.ImmutableSequenceOntology;
import com.hartwig.hmftools.vicc.datamodel.ImmutableTaxonomy;
import com.hartwig.hmftools.vicc.datamodel.ImmutableViccEntry;
import com.hartwig.hmftools.vicc.datamodel.Phenotype;
import com.hartwig.hmftools.vicc.datamodel.PhenotypeType;
import com.hartwig.hmftools.vicc.datamodel.SequenceOntology;
import com.hartwig.hmftools.vicc.datamodel.Taxonomy;
import com.hartwig.hmftools.vicc.datamodel.ViccEntry;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public final class ViccJsonReader {

    private static final Logger LOGGER = LogManager.getLogger(ViccJsonReader.class);

    private ViccJsonReader() {
    }

    @NotNull
    public static List<ViccEntry> readViccKnowledgebaseJsonFile(@NotNull String jsonPath) throws IOException {
        List<ViccEntry> entries = Lists.newArrayList();

        JsonParser parser = new JsonParser();
        JsonReader reader = new JsonReader(new FileReader(jsonPath));
        reader.setLenient(true);

        while (reader.peek() != JsonToken.END_DOCUMENT) {
            JsonObject viccEntryObject = parser.parse(reader).getAsJsonObject();
            ViccDatamodelCheckerFactory.viccEntryChecker().check(viccEntryObject);

            ImmutableViccEntry.Builder viccEntryBuilder = ImmutableViccEntry.builder();
            viccEntryBuilder.source(viccEntryObject.getAsJsonPrimitive("source").getAsString());
            viccEntryBuilder.genes(jsonArrayToStringList(viccEntryObject.getAsJsonArray("genes")));
            viccEntryBuilder.geneIdentifiers(createGeneIdentifiers(viccEntryObject.getAsJsonArray("gene_identifiers")));

            // SAGE records have no "feature names" while all other knowledgebases do have it.
            if (viccEntryObject.has("feature_names")) {
                JsonElement featureNames = viccEntryObject.get("feature_names");
                if (featureNames.isJsonArray()) {
                    viccEntryBuilder.featureNames(jsonArrayToStringList(featureNames.getAsJsonArray()));
                } else if (featureNames.isJsonPrimitive()) {
                    viccEntryBuilder.featureNames(Lists.newArrayList(featureNames.getAsJsonPrimitive().getAsString()));
                }
            }

            viccEntryBuilder.features(createFeatures(viccEntryObject.getAsJsonArray("features")));
            viccEntryBuilder.association(createAssociation(viccEntryObject.getAsJsonObject("association")));
            viccEntryBuilder.tags(jsonArrayToStringList(viccEntryObject.getAsJsonArray("tags")));
            viccEntryBuilder.devTags(jsonArrayToStringList(viccEntryObject.getAsJsonArray("dev_tags")));

            if (viccEntryObject.has("cgi")) {
                viccEntryBuilder.KbSpecificObject(CgiObjectFactory.create(viccEntryObject.getAsJsonObject("cgi")));
            } else if (viccEntryObject.has("brca")) {
                viccEntryBuilder.KbSpecificObject(BRCAObjectFactory.create(viccEntryObject.getAsJsonObject("brca")));
            } else if (viccEntryObject.has("sage")) {
                viccEntryBuilder.KbSpecificObject(SageObjectFactory.create(viccEntryObject.getAsJsonObject("sage")));
            } else if (viccEntryObject.has("pmkb")) {
                viccEntryBuilder.KbSpecificObject(PmkbObjectFactory.create(viccEntryObject.getAsJsonObject("pmkb")));
            } else if (viccEntryObject.has("oncokb")) {
                viccEntryBuilder.KbSpecificObject(OncokbObjectFactory.create(viccEntryObject.getAsJsonObject("oncokb")));
            } else if (viccEntryObject.has("jax")) {
                viccEntryBuilder.KbSpecificObject(JaxObjectFactory.create(viccEntryObject.getAsJsonObject("jax")));
            } else if (viccEntryObject.has("jax_trials")) {
                viccEntryBuilder.KbSpecificObject(JaxTrialsObjectFactory.create(viccEntryObject.getAsJsonObject("jax_trials")));
            } else if (viccEntryObject.has("molecularmatch")) {
                viccEntryBuilder.KbSpecificObject(MolecularMatchObjectFactory.create(viccEntryObject.getAsJsonObject("molecularmatch")));
            } else if (viccEntryObject.has("molecularmatch_trials")) {
                viccEntryBuilder.KbSpecificObject(MolecularMatchTrialsObjectFactory.create(viccEntryObject.getAsJsonObject(
                        "molecularmatch_trials")));
            } else if (viccEntryObject.has("civic")) {
                viccEntryBuilder.KbSpecificObject(CivicObjectFactory.create(viccEntryObject.getAsJsonObject("civic")));
            } else {
                LOGGER.warn("Could not resolve kb specific object for {}", viccEntryObject);
            }

            entries.add(viccEntryBuilder.build());
        }

        reader.close();

        return entries;
    }

    @NotNull
    private static List<GeneIdentifier> createGeneIdentifiers(@NotNull JsonArray geneIdentifierArray) {
        List<GeneIdentifier> geneIdentifierList = Lists.newArrayList();

        for (JsonElement geneIdentifierElement : geneIdentifierArray) {
            JsonObject geneIdentifierObject = geneIdentifierElement.getAsJsonObject();
            ViccDatamodelCheckerFactory.geneIdentifierChecker().check(geneIdentifierObject);

            geneIdentifierList.add(ImmutableGeneIdentifier.builder()
                    .symbol(geneIdentifierObject.getAsJsonPrimitive("symbol").getAsString())
                    .entrezId(geneIdentifierObject.getAsJsonPrimitive("entrez_id").getAsString())
                    .ensemblGeneId(!geneIdentifierObject.get("ensembl_gene_id").isJsonNull() ? geneIdentifierObject.getAsJsonPrimitive(
                            "ensembl_gene_id").getAsString() : null)
                    .build());
        }

        return geneIdentifierList;
    }

    @NotNull
    private static List<Feature> createFeatures(@NotNull JsonArray featureArray) {
        List<Feature> featureList = Lists.newArrayList();
        ViccDatamodelChecker featureChecker = ViccDatamodelCheckerFactory.featureChecker();

        for (JsonElement featureElement : featureArray) {
            JsonObject featureObject = featureElement.getAsJsonObject();
            featureChecker.check(featureObject);

            featureList.add(ImmutableFeature.builder()
                    .name(optionalString(featureObject, "name"))
                    .biomarkerType(optionalString(featureObject, "biomarker_type"))
                    .referenceName(optionalString(featureObject, "referenceName"))
                    .chromosome(optionalString(featureObject, "chromosome"))
                    .start(optionalNullableString(featureObject, "start"))
                    .end(optionalNullableString(featureObject, "end"))
                    .ref(optionalNullableString(featureObject, "ref"))
                    .alt(optionalNullableString(featureObject, "alt"))
                    .provenance(optionalStringList(featureObject, "provenance"))
                    .provenanceRule(optionalString(featureObject, "provenance_rule"))
                    .geneSymbol(optionalNullableString(featureObject, "geneSymbol"))
                    .synonyms(optionalStringList(featureObject, "synonyms"))
                    .entrezId(optionalString(featureObject, "entrez_id"))
                    .sequenceOntology(featureObject.has("sequence_ontology") ? createSequenceOntology(featureObject.getAsJsonObject(
                            "sequence_ontology")) : null)
                    .links(optionalStringList(featureObject, "links"))
                    .description(optionalString(featureObject, "description"))
                    .info(featureObject.has("info") ? createFeatureInfo(featureObject.getAsJsonObject("info")) : null)
                    .attribute(featureObject.has("attributes") ? createFeatureAttribute(featureObject.getAsJsonObject("attributes")) : null)
                    .build());
        }

        return featureList;
    }

    @NotNull
    private static FeatureInfo createFeatureInfo(@NotNull JsonObject featureInfoObject) {
        ViccDatamodelCheckerFactory.featureInfoChecker().check(featureInfoObject);

        return ImmutableFeatureInfo.builder().germlineOrSomatic(featureInfoObject.get("germline_or_somatic").getAsString()).build();
    }

    @NotNull
    private static FeatureAttribute createFeatureAttribute(@NotNull JsonObject featureAttributeObject) {
        ViccDatamodelCheckerFactory.featureAttributeChecker().check(featureAttributeObject);

        return ImmutableFeatureAttribute.builder()
                .aminoAcidChange(extractStringValueFromAttribute(featureAttributeObject.getAsJsonObject("amino_acid_change")))
                .germline(extractStringValueFromAttribute(featureAttributeObject.getAsJsonObject("germline")))
                .partnerGene(extractStringValueFromAttribute(featureAttributeObject.getAsJsonObject("partner_gene")))
                .description(extractStringValueFromAttribute(featureAttributeObject.getAsJsonObject("description")))
                .exons(extractStringValueFromAttribute(featureAttributeObject.getAsJsonObject("exons")))
                .notes(extractStringValueFromAttribute(featureAttributeObject.getAsJsonObject("notes")))
                .cosmic(extractStringValueFromAttribute(featureAttributeObject.getAsJsonObject("cosmic")))
                .effect(extractStringValueFromAttribute(featureAttributeObject.getAsJsonObject("effect")))
                .cnvType(extractStringValueFromAttribute(featureAttributeObject.getAsJsonObject("cnv_type")))
                .id(extractStringValueFromAttribute(featureAttributeObject.getAsJsonObject("id")))
                .cytoband(extractStringValueFromAttribute(featureAttributeObject.getAsJsonObject("cytoband")))
                .variantType(extractStringValueFromAttribute(featureAttributeObject.getAsJsonObject("variant_type")))
                .dnaChange(extractStringValueFromAttribute(featureAttributeObject.getAsJsonObject("dna_change")))
                .codons(extractStringValueFromAttribute(featureAttributeObject.getAsJsonObject("codons")))
                .chromosomeBasedCnv(extractStringValueFromAttribute(featureAttributeObject.getAsJsonObject("chromosome_based_cnv")))
                .transcript(extractStringValueFromAttribute(featureAttributeObject.getAsJsonObject("transcript")))
                .descriptionType(extractStringValueFromAttribute(featureAttributeObject.getAsJsonObject("description_type")))
                .chromosome(extractStringValueFromAttribute(featureAttributeObject.getAsJsonObject("chromosome")))
                .build();
    }

    @Nullable
    private static String extractStringValueFromAttribute(@NotNull JsonObject attributeObject) {
        ViccDatamodelCheckerFactory.featureAttributeObjectChecker().check(attributeObject);

        return nullableString(attributeObject, "string_value");
    }

    @NotNull
    private static SequenceOntology createSequenceOntology(@NotNull JsonObject sequenceOntologyObject) {
        ViccDatamodelCheckerFactory.sequenceOntologyChecker().check(sequenceOntologyObject);

        return ImmutableSequenceOntology.builder()
                .hierarchy(sequenceOntologyObject.has("hierarchy")
                        ? jsonArrayToStringList(sequenceOntologyObject.getAsJsonArray("hierarchy"))
                        : null)
                .soid(sequenceOntologyObject.getAsJsonPrimitive("soid").getAsString())
                .parentSoid(sequenceOntologyObject.getAsJsonPrimitive("parent_soid").getAsString())
                .name(sequenceOntologyObject.getAsJsonPrimitive("name").getAsString())
                .parentName(sequenceOntologyObject.getAsJsonPrimitive("parent_name").getAsString())
                .build();
    }

    @NotNull
    private static Association createAssociation(@NotNull JsonObject associationObject) {
        ViccDatamodelCheckerFactory.associationChecker().check(associationObject);

        return ImmutableAssociation.builder()
                .variantName(associationObject.has("variant_name") && associationObject.get("variant_name").isJsonArray()
                        ? Strings.EMPTY
                        : associationObject.has("variant_name") && associationObject.get("variant_name").isJsonPrimitive()
                                ? associationObject.getAsJsonPrimitive("variant_name").getAsString()
                                : null)
                .evidence(createEvidence(associationObject.getAsJsonArray("evidence")))
                .evidenceLevel(associationObject.has("evidence_level") ? associationObject.getAsJsonPrimitive("evidence_level")
                        .getAsString() : null)
                .evidenceLabel(
                        associationObject.has("evidence_label") && !associationObject.get("evidence_label").isJsonNull() ? associationObject
                                .getAsJsonPrimitive("evidence_label")
                                .getAsString() : null)
                .responseType(associationObject.has("response_type") && !associationObject.get("response_type").isJsonNull()
                        ? associationObject.getAsJsonPrimitive("response_type").getAsString()
                        : null)
                .drugLabels(associationObject.has("drug_labels") ? associationObject.getAsJsonPrimitive("drug_labels").getAsString() : null)
                .sourceLink(associationObject.has("source_link") ? associationObject.getAsJsonPrimitive("source_link").getAsString() : null)
                .publicationUrls(associationObject.has("publication_url") && associationObject.get("publication_url").isJsonPrimitive()
                        ? Lists.newArrayList(associationObject.getAsJsonPrimitive("publication_url").getAsString())
                        : associationObject.has("publication_url") && associationObject.get("publication_url").isJsonArray()
                                ? Lists.newArrayList(associationObject.getAsJsonArray("publication_url").getAsString())
                                : null)
                .phenotype(associationObject.has("phenotype") ? createPhenotype(associationObject.getAsJsonObject("phenotype")) : null)
                .description(associationObject.getAsJsonPrimitive("description").getAsString())
                .environmentalContexts(associationObject.has("environmentalContexts")
                        ? createEnvironmentalContexts(associationObject.getAsJsonArray("environmentalContexts"))
                        : null)
                .oncogenic(associationObject.has("oncogenic") ? associationObject.getAsJsonPrimitive("oncogenic").getAsString() : null)
                .build();
    }

    @NotNull
    private static List<Evidence> createEvidence(@NotNull JsonArray evidenceArray) {
        List<Evidence> evidenceList = Lists.newArrayList();
        ViccDatamodelChecker evidenceChecker = ViccDatamodelCheckerFactory.evidenceChecker();

        for (JsonElement evidenceElement : evidenceArray) {
            JsonObject evidenceObject = evidenceElement.getAsJsonObject();
            evidenceChecker.check(evidenceObject);

            evidenceList.add(ImmutableEvidence.builder()
                    .info(!evidenceObject.get("info").isJsonNull() ? createEvidenceInfo(evidenceObject.getAsJsonObject("info")) : null)
                    .evidenceType(createEvidenceType(evidenceObject.getAsJsonObject("evidenceType")))
                    .description(!evidenceObject.get("description").isJsonNull() ? evidenceObject.getAsJsonPrimitive("description")
                            .getAsString() : null)
                    .build());
        }
        return evidenceList;
    }

    @NotNull
    private static EvidenceInfo createEvidenceInfo(@NotNull JsonObject evidenceInfoObject) {
        ViccDatamodelCheckerFactory.evidenceInfoChecker().check(evidenceInfoObject);

        return ImmutableEvidenceInfo.builder()
                .publications(jsonArrayToStringList(evidenceInfoObject.getAsJsonArray("publications")))
                .build();
    }

    @NotNull
    private static EvidenceType createEvidenceType(@NotNull JsonObject evidenceTypeObject) {
        ViccDatamodelCheckerFactory.evidenceTypeChecker().check(evidenceTypeObject);

        return ImmutableEvidenceType.builder()
                .sourceName(evidenceTypeObject.getAsJsonPrimitive("sourceName").getAsString())
                .id(evidenceTypeObject.has("id") ? evidenceTypeObject.getAsJsonPrimitive("id").getAsString() : null)
                .build();
    }

    @NotNull
    private static List<EnvironmentalContext> createEnvironmentalContexts(@NotNull JsonArray environmentalContextArray) {
        List<EnvironmentalContext> environmentalContextList = Lists.newArrayList();
        ViccDatamodelChecker environmentalContextChecker = ViccDatamodelCheckerFactory.environmentalContextChecker();

        for (JsonElement environmentalContextElement : environmentalContextArray) {
            JsonObject environmentalContextObject = environmentalContextElement.getAsJsonObject();
            environmentalContextChecker.check(environmentalContextObject);

            environmentalContextList.add(ImmutableEnvironmentalContext.builder()
                    .term(environmentalContextObject.has("term")
                            ? environmentalContextObject.getAsJsonPrimitive("term").getAsString()
                            : null)
                    .description(environmentalContextObject.getAsJsonPrimitive("description").getAsString())
                    .taxonomy(environmentalContextObject.has("taxonomy") ? createTaxonomy(environmentalContextObject.getAsJsonObject(
                            "taxonomy")) : null)
                    .source(environmentalContextObject.has("source")
                            ? environmentalContextObject.getAsJsonPrimitive("source").getAsString()
                            : null)
                    .usanStem(environmentalContextObject.has("usan_stem") ? environmentalContextObject.getAsJsonPrimitive("usan_stem")
                            .getAsString() : null)
                    .approvedCountries(environmentalContextObject.has("approved_countries") ? jsonArrayToStringList(
                            environmentalContextObject.getAsJsonArray("approved_countries")) : Lists.newArrayList())
                    .toxicity(environmentalContextObject.has("toxicity") ? environmentalContextObject.getAsJsonPrimitive("toxicity")
                            .getAsString() : null)
                    .id(environmentalContextObject.has("id") && !environmentalContextObject.get("id").isJsonNull()
                            ? environmentalContextObject.getAsJsonPrimitive("id").getAsString()
                            : null)
                    .build());
        }
        return environmentalContextList;
    }

    @NotNull
    private static Taxonomy createTaxonomy(@NotNull JsonObject taxonomyObject) {
        ViccDatamodelCheckerFactory.taxonomyChecker().check(taxonomyObject);

        return ImmutableTaxonomy.builder()
                .kingdom(taxonomyObject.getAsJsonPrimitive("kingdom").getAsString())
                .directParent(taxonomyObject.getAsJsonPrimitive("direct-parent").getAsString())
                .classs(taxonomyObject.getAsJsonPrimitive("class").getAsString())
                .subClass(taxonomyObject.has("subclass") ? taxonomyObject.getAsJsonPrimitive("subclass").getAsString() : null)
                .superClass(taxonomyObject.getAsJsonPrimitive("superclass").getAsString())
                .build();
    }

    @NotNull
    private static Phenotype createPhenotype(@NotNull JsonObject phenotypeObject) {
        ViccDatamodelCheckerFactory.phenotypeChecker().check(phenotypeObject);

        return ImmutablePhenotype.builder()
                .type(phenotypeObject.has("type") ? createPhenotypeType(phenotypeObject.getAsJsonObject("type")) : null)
                .description(phenotypeObject.getAsJsonPrimitive("description").getAsString())
                .family(phenotypeObject.getAsJsonPrimitive("family").getAsString())
                .id(phenotypeObject.has("id") ? phenotypeObject.getAsJsonPrimitive("id").getAsString() : null)
                .build();
    }

    @NotNull
    private static PhenotypeType createPhenotypeType(@NotNull JsonObject phenotypeTypeObject) {
        ViccDatamodelCheckerFactory.phenotypeTypeChecker().check(phenotypeTypeObject);

        return ImmutablePhenotypeType.builder()
                .source(!phenotypeTypeObject.get("source").isJsonNull()
                        ? phenotypeTypeObject.getAsJsonPrimitive("source").getAsString()
                        : null)
                .term(phenotypeTypeObject.getAsJsonPrimitive("term").getAsString())
                .id(phenotypeTypeObject.getAsJsonPrimitive("id").getAsString())
                .build();
    }
}
