package com.hartwig.hmftools.vicc;

import java.io.FileReader;
import java.io.IOException;
import java.util.List;
import java.util.Set;

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
import com.hartwig.hmftools.vicc.datamodel.GeneIdentifier;
import com.hartwig.hmftools.vicc.datamodel.ImmutableAssociation;
import com.hartwig.hmftools.vicc.datamodel.ImmutableEnvironmentalContext;
import com.hartwig.hmftools.vicc.datamodel.ImmutableEvidence;
import com.hartwig.hmftools.vicc.datamodel.ImmutableEvidenceInfo;
import com.hartwig.hmftools.vicc.datamodel.ImmutableEvidenceType;
import com.hartwig.hmftools.vicc.datamodel.ImmutableGeneIdentifier;
import com.hartwig.hmftools.vicc.datamodel.ImmutablePhenotype;
import com.hartwig.hmftools.vicc.datamodel.ImmutablePhenotypeType;
import com.hartwig.hmftools.vicc.datamodel.ImmutableTaxonomy;
import com.hartwig.hmftools.vicc.datamodel.ImmutableViccEntry;
import com.hartwig.hmftools.vicc.datamodel.Phenotype;
import com.hartwig.hmftools.vicc.datamodel.PhenotypeType;
import com.hartwig.hmftools.vicc.datamodel.Taxonomy;
import com.hartwig.hmftools.vicc.datamodel.ViccEntry;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public final class ViccFactory {
    private static final Logger LOGGER = LogManager.getLogger(ViccFactory.class);

    // SAGE records hold 8 field (no "feature names") while all other knowledgebases hold 9 records.
    private static final List<Integer> EXPECTED_VICC_ENTRY_SIZES = Lists.newArrayList(8, 9);

    private static final List<Integer> EXPECTED_ASSOCIATION_ELEMENT_SIZES = Lists.newArrayList(9, 10);

    private ViccFactory() {
    }

    public static List<ViccEntry> readViccKnowledgebaseJsonFile(@NotNull String jsonPath) throws IOException {
        JsonParser parser = new JsonParser();
        JsonReader reader = new JsonReader(new FileReader(jsonPath));
        reader.setLenient(true);

        List<ViccEntry> entries = Lists.newArrayList();
        LOGGER.info("Reading VICC knowledgebase from " + jsonPath);
        while (reader.peek() != JsonToken.END_DOCUMENT) {
            JsonObject viccEntryElement = parser.parse(reader).getAsJsonObject();

            if (!EXPECTED_VICC_ENTRY_SIZES.contains(viccEntryElement.size())) {
                LOGGER.warn("Found " + viccEntryElement.size() + " elements in a vicc entry rather than the expected "
                        + EXPECTED_VICC_ENTRY_SIZES);
                LOGGER.warn(viccEntryElement);
            }

            ImmutableViccEntry.Builder viccEntryBuilder = ImmutableViccEntry.builder();
            viccEntryBuilder.source(viccEntryElement.getAsJsonPrimitive("source").getAsString());
            viccEntryBuilder.genes(jsonArrayToStringList(viccEntryElement.getAsJsonArray("genes")));

            viccEntryBuilder.geneIdentifiers(createGeneIdentifiers(viccEntryElement));

            if (viccEntryElement.has("feature_names")) {
                JsonElement featureNames = viccEntryElement.get("feature_names");
                if (featureNames.isJsonArray()) {
                    viccEntryBuilder.featureNames(jsonArrayToStringList(featureNames.getAsJsonArray()));
                } else if (featureNames.isJsonPrimitive()) {
                    viccEntryBuilder.featureNames(Lists.newArrayList(featureNames.getAsJsonPrimitive().getAsString()));
                }
            }

            viccEntryBuilder.features(Lists.newArrayList());

            JsonElement elementAssociation = viccEntryElement.get("association");
            Set<String> keysAssociation = elementAssociation.getAsJsonObject().keySet();

            if (!EXPECTED_ASSOCIATION_ELEMENT_SIZES.contains(keysAssociation.size())) {
                LOGGER.warn("Found " + keysAssociation.size() + " elements in a vicc entry rather than the expected "
                        + EXPECTED_ASSOCIATION_ELEMENT_SIZES);
                LOGGER.warn(keysAssociation);
            } else {
                viccEntryBuilder.association(createAssociation(viccEntryElement));
            }

            viccEntryBuilder.tags(jsonArrayToStringList(viccEntryElement.getAsJsonArray("tags")));
            viccEntryBuilder.devTags(jsonArrayToStringList(viccEntryElement.getAsJsonArray("dev_tags")));

            entries.add(viccEntryBuilder.build());
        }
        reader.close();

        return entries;
    }

    @NotNull
    private static List<String> jsonArrayToStringList(@NotNull JsonArray array) {
        List<String> values = Lists.newArrayList();
        for (JsonElement element : array) {
            values.add(element.getAsString());
        }
        return values;
    }

    @NotNull
    private static List<GeneIdentifier> createGeneIdentifiers(JsonObject viccEntryElement) {
        JsonArray geneIdentifiers = viccEntryElement.getAsJsonArray("gene_identifiers");
        List<GeneIdentifier> listGeneIdentifiers = Lists.newArrayList();

        for (JsonElement elementGeneIdentifier : geneIdentifiers) {
            listGeneIdentifiers.add(toGeneIdentifier(elementGeneIdentifier.getAsJsonObject()));

        }
        return listGeneIdentifiers;
    }

    @NotNull
    private static GeneIdentifier toGeneIdentifier(@NotNull JsonObject objectGeneIdentifier) {
        return ImmutableGeneIdentifier.builder()
                .symbol(objectGeneIdentifier.getAsJsonPrimitive("symbol").toString())
                .entrezId(objectGeneIdentifier.getAsJsonPrimitive("entrez_id").toString())
                .ensemblGeneId(objectGeneIdentifier.getAsJsonPrimitive("ensembl_gene_id").toString())
                .build();
    }

    @NotNull
    private static Association createAssociation(JsonObject viccEntryElement) {
        JsonObject objectAssociation = viccEntryElement.getAsJsonObject("association");
        return ImmutableAssociation.builder()
                .variantName(objectAssociation.has("variant_name") ? objectAssociation.getAsJsonObject()
                        .getAsJsonPrimitive("variant_name")
                        .toString() : null)

                .evidence(Lists.newArrayList(createEvidence(objectAssociation)))
                .evidenceLevel(objectAssociation.getAsJsonPrimitive("evidence_level").getAsString())
                .evidenceLabel(objectAssociation.getAsJsonPrimitive("evidence_label").getAsString())
                .responseType(objectAssociation.getAsJsonPrimitive("response_type").getAsString())
                .drugLabels(objectAssociation.getAsJsonPrimitive("drug_labels").getAsString())
                .sourceLink(objectAssociation.has("source_link") ? objectAssociation.getAsJsonPrimitive("source_link").getAsString() : null)
                .publicationUrls(Lists.newArrayList(objectAssociation.getAsJsonObject()
                        .getAsJsonPrimitive("publication_url")
                        .getAsString()))
                .phenotype(createPhenotype(objectAssociation))
                .description(objectAssociation.getAsJsonPrimitive("description").getAsString())
                .environmentalContexts(createEnvironmentalContexts(objectAssociation))
                .oncogenic(objectAssociation.has("oncogenic") ? objectAssociation.getAsJsonPrimitive("oncogenic").getAsString() : null)
                .build();
    }

    @NotNull
    private static List<EnvironmentalContext> createEnvironmentalContexts(@NotNull JsonObject objectAssociation) {
        List<EnvironmentalContext> environmentalContexts = Lists.newArrayList();

        for (JsonElement elementEnvironmentContext : objectAssociation.getAsJsonArray("environmentalContexts")) {
            JsonObject objectEnvironmentContext = elementEnvironmentContext.getAsJsonObject();

            List<String> approvedCountries = Lists.newArrayList();
            if (objectEnvironmentContext.has("approved_countries")) {
                for (JsonElement approvedCountriesElement : objectEnvironmentContext.getAsJsonArray("approved_countries")) {
                    approvedCountries.add(approvedCountriesElement.getAsString());
                }
            }

            environmentalContexts.add(ImmutableEnvironmentalContext.builder()
                    .term(objectEnvironmentContext.has("term") ? objectEnvironmentContext.get("term").getAsString() : null)
                    .description(objectEnvironmentContext.getAsJsonPrimitive("description").getAsString())
                    .taxonomy(objectEnvironmentContext.has("taxonomy")
                            ? createTaxonomy(objectEnvironmentContext.getAsJsonObject("taxonomy"))
                            : null)
                    .source(objectEnvironmentContext.has("source")
                            ? objectEnvironmentContext.getAsJsonPrimitive("source").getAsString()
                            : null)
                    .usanStem(objectEnvironmentContext.has("usan_stem") ? objectEnvironmentContext.getAsJsonPrimitive("usan_stem")
                            .getAsString() : null)
                    .approvedCountries(approvedCountries)
                    .id(objectEnvironmentContext.has("id") ? objectEnvironmentContext.getAsJsonPrimitive("id").getAsString() : null)
                    .build());
        }
        return environmentalContexts;
    }

    private static Taxonomy createTaxonomy(JsonObject elementEnvironmentContext) {
        return ImmutableTaxonomy.builder()
                .kingdom(elementEnvironmentContext.getAsJsonPrimitive("kingdom").getAsString())
                .directParent(elementEnvironmentContext.getAsJsonPrimitive("direct-parent").getAsString())
                .classs(elementEnvironmentContext.getAsJsonPrimitive("class").getAsString())
                .subClass(elementEnvironmentContext.has("subclass")
                        ? elementEnvironmentContext.getAsJsonPrimitive("subclass").getAsString()
                        : null)
                .superClass(elementEnvironmentContext.getAsJsonPrimitive("superclass").getAsString())
                .build();
    }

    private static Evidence createEvidence(JsonElement elementAssociation) {
        ImmutableEvidence builderEvidence = ImmutableEvidence.builder()
                .info(ImmutableEvidenceInfo.builder().publications(Lists.newArrayList()).build())
                .evidenceType(ImmutableEvidenceType.builder().id(Strings.EMPTY).sourceName(Strings.EMPTY).build())
                .description(Strings.EMPTY)
                .build();

        for (JsonElement elementEvidence : elementAssociation.getAsJsonObject().getAsJsonPrimitive("evidence").getAsJsonArray()) {
            builderEvidence = ImmutableEvidence.builder()
                    .info(createEvidenceInfo(elementEvidence.getAsJsonObject().getAsJsonPrimitive("info")))
                    .evidenceType(createEvidenceType(elementEvidence.getAsJsonObject().getAsJsonPrimitive("evidenceType")))
                    .description(elementEvidence.getAsJsonObject().getAsJsonPrimitive("description").getAsString())
                    .build();
        }
        return builderEvidence;
    }

    private static EvidenceType createEvidenceType(JsonElement elementAssociation) {
        return ImmutableEvidenceType.builder()
                .sourceName(elementAssociation.getAsJsonObject().getAsJsonPrimitive("sourceName").toString())
                .id(elementAssociation.getAsJsonObject().has("id") ? elementAssociation.getAsJsonObject().getAsJsonPrimitive("id").toString() : null)
                .build();
    }

    private static EvidenceInfo createEvidenceInfo(JsonElement elementAssociation) {
        return ImmutableEvidenceInfo.builder()
                .publications(Lists.newArrayList(elementAssociation.getAsJsonObject().getAsJsonPrimitive("publications").toString()))
                .build();
    }

    private static PhenotypeType createPhenotypeType(JsonElement elementAssociation) {
        return ImmutablePhenotypeType.builder()
                .source(elementAssociation.getAsJsonObject().getAsJsonPrimitive("source").getAsString())
                .term(elementAssociation.getAsJsonObject().getAsJsonPrimitive("term").getAsString())
                .id(elementAssociation.getAsJsonObject().getAsJsonPrimitive("id").getAsString())
                .build();
    }

    private static Phenotype createPhenotype(JsonElement elementAssociation) {
        return ImmutablePhenotype.builder()
                .type(createPhenotypeType(elementAssociation.getAsJsonObject().getAsJsonObject("phenotype").getAsJsonPrimitive("type")))
                .description(elementAssociation.getAsJsonObject().getAsJsonObject("phenotype").getAsJsonPrimitive("description").toString())
                .family(elementAssociation.getAsJsonObject().getAsJsonObject("phenotype").getAsJsonPrimitive("family").toString())
                .build();
    }
}
