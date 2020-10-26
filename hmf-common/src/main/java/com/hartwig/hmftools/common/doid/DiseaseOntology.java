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

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class DiseaseOntology {

    private DiseaseOntology() {
    }

    @NotNull
    public static List<DoidEntry> readDoidJsonFile(@NotNull String doidJsonFile) throws IOException {
        JsonParser parser = new JsonParser();
        JsonReader reader = new JsonReader(new FileReader(doidJsonFile));
        reader.setLenient(true);
        List<DoidEntry> doids = Lists.newArrayList();
        JsonDatamodelChecker doidEntryChecker = DoidDatamodelCheckerFactory.doidEntryChecker();

        while (reader.peek() != JsonToken.END_DOCUMENT) {
            JsonObject doidObject = parser.parse(reader).getAsJsonObject();
            doidEntryChecker.check(doidObject);
            JsonDatamodelChecker doidEntryGraphChecker = DoidDatamodelCheckerFactory.doidEntryGraphChecker();

            JsonArray graphArray = doidObject.getAsJsonArray("graphs");

            for (JsonElement graph : graphArray) {
                doidEntryGraphChecker.check(graph.getAsJsonObject());
                JsonArray nodeArray = graph.getAsJsonObject().getAsJsonArray("nodes");

                JsonDatamodelChecker doidEntryNodesChecker = DoidDatamodelCheckerFactory.doidEntryNodesChecker();

                for (JsonElement nodeElement : nodeArray) {
                    JsonObject node = nodeElement.getAsJsonObject();
                    doidEntryNodesChecker.check(node);
                    String url = string(node, "id");
                    doids.add(ImmutableDoidEntry.builder()
                            .doid(extractDoid(url))
                            .url(url)
                            .doidMetadata(extractDoidMetadata(optionalJsonObject(node, "meta")))
                            .type(optionalString(node, "type"))
                            .doidTerm(optionalString(node, "lbl"))
                            .idNodes("")
                            .edges(Lists.newArrayList())
                            .metaNodes(Lists.newArrayList())
                            .equivalentNodesSets(Lists.newArrayList())
                            .logicalDefinitionAxioms(Lists.newArrayList())
                            .domainRangeAxioms(Lists.newArrayList())
                            .propertyChainAxioms(Lists.newArrayList())
                            .build());
                }
            }
        }
        return doids;
    }

    @NotNull
    public static DoidEdges readDoidJsonFileEdges(@NotNull String doidJsonFile) throws IOException {
        DoidEdges result = new DoidEdges();

        JsonParser parser = new JsonParser();
        JsonReader reader = new JsonReader(new FileReader(doidJsonFile));
        reader.setLenient(true);
        while (reader.peek() != JsonToken.END_DOCUMENT) {
            JsonObject doidObject = parser.parse(reader).getAsJsonObject();

            JsonArray graphArray = doidObject.getAsJsonArray("graphs");
            for (JsonElement graph : graphArray) {
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
            }
        }

        return result;
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
