package com.hartwig.hmftools.common.doid;

import static com.hartwig.hmftools.common.utils.json.JsonFunctions.optionalJsonArray;
import static com.hartwig.hmftools.common.utils.json.JsonFunctions.optionalJsonObject;
import static com.hartwig.hmftools.common.utils.json.JsonFunctions.optionalString;
import static com.hartwig.hmftools.common.utils.json.JsonFunctions.optionalStringList;
import static com.hartwig.hmftools.common.utils.json.JsonFunctions.string;
import static com.hartwig.hmftools.common.utils.json.JsonFunctions.stringList;

import java.io.FileReader;
import java.io.IOException;
import java.util.List;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.google.gson.JsonArray;
import com.google.gson.JsonElement;
import com.google.gson.JsonObject;
import com.google.gson.JsonParser;
import com.google.gson.stream.JsonReader;
import com.google.gson.stream.JsonToken;
import com.hartwig.hmftools.common.utils.json.JsonDatamodelChecker;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class DiseaseOntology {

    private DiseaseOntology() {
    }

    private static final Logger LOGGER = LogManager.getLogger(DiseaseOntology.class);

    @NotNull
    public static DoidEntry readDoidJsonFile(@NotNull String doidJsonFile) throws IOException {
        JsonParser parser = new JsonParser();
        JsonReader reader = new JsonReader(new FileReader(doidJsonFile));
        reader.setLenient(true);
        List<DoidNode> doidNodes = Lists.newArrayList();
        ImmutableDoidEntry.Builder doidDoidEntryBuilder = ImmutableDoidEntry.builder();

        JsonDatamodelChecker doidEntryChecker = DoidDatamodelCheckerFactory.doidEntryChecker();
        while (reader.peek() != JsonToken.END_DOCUMENT) {
            JsonObject doidObject = parser.parse(reader).getAsJsonObject();
            doidEntryChecker.check(doidObject);
            JsonDatamodelChecker doidEntryGraphChecker = DoidDatamodelCheckerFactory.doidEntryGraphChecker();

            JsonArray graphArray = doidObject.getAsJsonArray("graphs");
            if (graphArray.size() == 1) {
                for (JsonElement graph : graphArray) {
                    DoidEdges result = new DoidEdges();

                    doidEntryGraphChecker.check(graph.getAsJsonObject());
                    JsonArray nodeArray = graph.getAsJsonObject().getAsJsonArray("nodes");

                    JsonDatamodelChecker doidEntryNodesChecker = DoidDatamodelCheckerFactory.doidEntryNodesChecker();

                    JsonArray edgesArray = graph.getAsJsonObject().getAsJsonArray("edges");
                    for (JsonElement edgeElement : edgesArray) {
                        JsonObject edge = edgeElement.getAsJsonObject();
                        String predicate = string(edge, "pred");
                        if (predicate.equals("is_a")) {
                            String child = extractDoid(string(edge, "sub"));
                            String parent = extractDoid(string(edge, "obj"));
                            result.isA(child, parent);
                        }
                    }

                    //createMetaNodes(graph.getAsJsonObject().getAsJsonObject("meta"));
                    List<DoidEquivalentNodesSets> doidEquivalentNodesSets =
                            extractDoidEquivalentNodesSets(graph.getAsJsonObject().getAsJsonArray("equivalentNodesSets"));
                    List<DoidLogicalDefinitionAxioms> doidLogicalDefinitionAxioms =
                            extractDoidLogicalDefinitionAxioms(graph.getAsJsonObject().getAsJsonArray("logicalDefinitionAxioms"));
                    List<DoidDomainRangeAxioms> doidDomainRangeAxioms =
                            extractDoidDomainRangeAxioms(graph.getAsJsonObject().getAsJsonArray("domainRangeAxioms"));
                    List<DoidPropertyChainAxioms> doidPropertyChainAxioms =
                            extractDoidPropertyChainAxioms(graph.getAsJsonObject().getAsJsonArray("propertyChainAxioms"));

                    for (JsonElement nodeElement : nodeArray) {
                        JsonObject node = nodeElement.getAsJsonObject();
                        doidEntryNodesChecker.check(node);
                        String url = string(node, "id");
                        doidNodes.add(ImmutableDoidNode.builder()
                                .doid(extractDoid(url))
                                .url(url)
                                .doidMetadata(extractDoidMetadata(optionalJsonObject(node, "meta")))
                                .type(optionalString(node, "type"))
                                .doidTerm(optionalString(node, "lbl"))
                                .build());
                    }
                    doidDoidEntryBuilder.doidNodes(doidNodes);
                    doidDoidEntryBuilder.edges(result);
                    doidDoidEntryBuilder.id("");
                    doidDoidEntryBuilder.equivalentNodesSets(doidEquivalentNodesSets);
                    doidDoidEntryBuilder.meta(Lists.newArrayList());
                    doidDoidEntryBuilder.logicalDefinitionAxioms(doidLogicalDefinitionAxioms);
                    doidDoidEntryBuilder.domainRangeAxioms(doidDomainRangeAxioms);
                    doidDoidEntryBuilder.propertyChainAxioms(doidPropertyChainAxioms);
                }
            } else {
                LOGGER.info(" Size of graph elements are not correct {}", graphArray.size());
            }
        }
        return doidDoidEntryBuilder.build();
    }

    @Nullable
    private static List<String> createMetaNodes(@Nullable JsonObject metaNodesObject) {
        return Lists.newArrayList();

    }

    @Nullable
    private static List<DoidPropertyChainAxioms> extractDoidPropertyChainAxioms(@Nullable JsonArray propertyChainAxiomsArray) {
        if (propertyChainAxiomsArray == null) {
            return null;
        }

        List<DoidPropertyChainAxioms> doidPropertyChainAxioms = Lists.newArrayList();
        JsonDatamodelChecker doidPropertyChainAxiomsChecker = DoidDatamodelCheckerFactory.doidPropertyChainAxiomsChecker();

        for (JsonElement propertyChainAxioms : propertyChainAxiomsArray) {
            JsonObject propertyChainAxiomsObject = propertyChainAxioms.getAsJsonObject();
            doidPropertyChainAxiomsChecker.check(propertyChainAxiomsObject);

            doidPropertyChainAxioms.add(ImmutableDoidPropertyChainAxioms.builder().build());
        }

        return doidPropertyChainAxioms;
    }

    @Nullable
    private static List<DoidDomainRangeAxioms> extractDoidDomainRangeAxioms(@Nullable JsonArray domainRangeAxiomsArray) {
        if (domainRangeAxiomsArray == null) {
            return null;
        }

        List<DoidDomainRangeAxioms> doidDomainRangeAxioms = Lists.newArrayList();
        JsonDatamodelChecker doidDomainRangeAxiomsChecker = DoidDatamodelCheckerFactory.doidDomainRangeAxiomsChecker();

        for (JsonElement domainRangeAxioms : domainRangeAxiomsArray) {
            JsonObject domainRangeAxiomsObject = domainRangeAxioms.getAsJsonObject();
            doidDomainRangeAxiomsChecker.check(domainRangeAxiomsObject);

            doidDomainRangeAxioms.add(ImmutableDoidDomainRangeAxioms.builder().build());
        }

        return doidDomainRangeAxioms;
    }

    @Nullable
    private static List<DoidLogicalDefinitionAxioms> extractDoidLogicalDefinitionAxioms(@Nullable JsonArray logicalDefinitionAxiomsArray) {
        if (logicalDefinitionAxiomsArray == null) {
            return null;
        }

        List<DoidLogicalDefinitionAxioms> doidLogicalDefinitionAxioms = Lists.newArrayList();
        JsonDatamodelChecker doidLogicalDefinitionAxiomsChecker = DoidDatamodelCheckerFactory.doidLogicalDefinitionAxiomChecker();

        for (JsonElement logicalDefinitionAxioms : logicalDefinitionAxiomsArray) {
            JsonObject logicalDefinitionAxiomsObject = logicalDefinitionAxioms.getAsJsonObject();
            doidLogicalDefinitionAxiomsChecker.check(logicalDefinitionAxiomsObject);

            doidLogicalDefinitionAxioms.add(ImmutableDoidLogicalDefinitionAxioms.builder().build());
        }

        return doidLogicalDefinitionAxioms;
    }

    @Nullable
    private static List<DoidEquivalentNodesSets> extractDoidEquivalentNodesSets(@Nullable JsonArray equivalentNodesSetsArray) {
        if (equivalentNodesSetsArray == null) {
            return null;
        }

        List<DoidEquivalentNodesSets> doidEquivalentNodesSets = Lists.newArrayList();
        JsonDatamodelChecker doidEquivalentNodesSetsChecker = DoidDatamodelCheckerFactory.doidEquivalentNodesSetsChecker();

        for (JsonElement equivalentNodesSets : equivalentNodesSetsArray) {
            JsonObject equivalentNodesSetsObject = equivalentNodesSets.getAsJsonObject();
            doidEquivalentNodesSetsChecker.check(equivalentNodesSetsObject);

            doidEquivalentNodesSets.add(ImmutableDoidEquivalentNodesSets.builder().build());
        }

        return doidEquivalentNodesSets;
    }

    @VisibleForTesting
    @NotNull
    public static String extractDoid(@NotNull String url) {
        return url.replace("http://purl.obolibrary.org/obo/DOID_", "");
    }

    @Nullable
    private static DoidMetadata extractDoidMetadata(@Nullable JsonObject metadataObject) {
        if (metadataObject == null) {
            return null;
        }

        JsonDatamodelChecker metadataChecker = DoidDatamodelCheckerFactory.doidMetadataChecker();
        metadataChecker.check(metadataObject);

        JsonArray xrefArray = metadataObject.getAsJsonArray("xrefs");
        List<DoidXref> xrefValList = Lists.newArrayList();
        if (xrefArray != null) {
            JsonDatamodelChecker xrefChecker = DoidDatamodelCheckerFactory.doidMetadataXrefChecker();

            for (JsonElement xref : xrefArray) {
                xrefChecker.check(xref.getAsJsonObject());
                xrefValList.add(ImmutableDoidXref.builder().val(string(xref.getAsJsonObject(), "val")).build());
            }
        }

        ImmutableDoidMetadata.Builder doidMetadataBuilder = ImmutableDoidMetadata.builder();
        doidMetadataBuilder.synonyms(extractDoidSynonyms(optionalJsonArray(metadataObject, "synonyms")));
        doidMetadataBuilder.basicPropertyValues(extractBasicPropertyValues(optionalJsonArray(metadataObject, "basicPropertyValues")));
        doidMetadataBuilder.doidDefinition(extractDoidDefinition(optionalJsonObject(metadataObject, "definition")));
        doidMetadataBuilder.subsets(optionalStringList(metadataObject, "subsets"));
        doidMetadataBuilder.xrefs(xrefValList);
        return doidMetadataBuilder.build();
    }

    @Nullable
    private static List<DoidSynonym> extractDoidSynonyms(@Nullable JsonArray synonymArray) {
        if (synonymArray == null) {
            return null;
        }

        JsonDatamodelChecker doidSynonymsChecker = DoidDatamodelCheckerFactory.doidSynonymsChecker();

        List<DoidSynonym> doidSynonymList = Lists.newArrayList();
        for (JsonElement synonymElement : synonymArray) {
            JsonObject synonym = synonymElement.getAsJsonObject();
            doidSynonymsChecker.check(synonym);

            doidSynonymList.add(ImmutableDoidSynonym.builder()
                    .pred(string(synonym, "pred"))
                    .val(string(synonym, "val"))
                    .xrefs(stringList(synonym, "xrefs"))
                    .build());
        }

        return doidSynonymList;
    }

    @Nullable
    private static List<DoidBasicPropertyValue> extractBasicPropertyValues(@Nullable JsonArray basicPropertyValueArray) {
        if (basicPropertyValueArray == null) {
            return null;
        }

        List<DoidBasicPropertyValue> doidBasicPropertyValueList = Lists.newArrayList();
        JsonDatamodelChecker basicPropertyValuesChecker = DoidDatamodelCheckerFactory.doidBasicPropertyValuesChecker();

        for (JsonElement basicPropertyElement : basicPropertyValueArray) {
            JsonObject basicProperty = basicPropertyElement.getAsJsonObject();
            basicPropertyValuesChecker.check(basicProperty);

            doidBasicPropertyValueList.add(ImmutableDoidBasicPropertyValue.builder()
                    .pred(string(basicProperty, "pred"))
                    .val(string(basicProperty, "val"))
                    .build());
        }

        return doidBasicPropertyValueList;
    }

    @Nullable
    private static DoidDefinition extractDoidDefinition(@Nullable JsonObject definitionObject) {
        if (definitionObject == null) {
            return null;
        }

        JsonDatamodelChecker doidDefinitionChecker = DoidDatamodelCheckerFactory.doidDefinitionChecker();
        doidDefinitionChecker.check(definitionObject);

        ImmutableDoidDefinition.Builder doidDefinitionBuilder = ImmutableDoidDefinition.builder();
        doidDefinitionBuilder.definitionVal(string(definitionObject, "val"));
        doidDefinitionBuilder.definitionXrefs(stringList(definitionObject, "xrefs"));
        return doidDefinitionBuilder.build();
    }
}
