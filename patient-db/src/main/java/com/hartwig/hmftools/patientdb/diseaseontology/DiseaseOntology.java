package com.hartwig.hmftools.patientdb.diseaseontology;

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

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;


public class DiseaseOntology {
    private static final Logger LOGGER = LogManager.getLogger(DiseaseOntology.class);

    @Nullable
    private static DoidDefinition extractDoidDefinition(@NotNull JsonObject definitionJsonObject) {
        ImmutableDoidDefinition.Builder doidDefinition = ImmutableDoidDefinition.builder();
        doidDefinition.definitionVal(JsonFunctions.nullableString(definitionJsonObject, "val"));
        doidDefinition.definitionXref(JsonFunctions.stringList(definitionJsonObject, "xrefs"));
        return doidDefinition.build();

    }

    @NotNull
    private static List<DoidSynonyms> extractDoidSynonyms(@NotNull JsonArray arraySynonyms) {
        List<DoidSynonyms> doidSynonymsList = Lists.newArrayList();
        for (JsonElement synonyms : arraySynonyms) {
            doidSynonymsList.add(ImmutableDoidSynonyms.builder()
                    .pred(JsonFunctions.toString(synonyms.getAsJsonObject(), "pred"))
                    .val(JsonFunctions.toString(synonyms.getAsJsonObject(), "val"))
                    .xrefs(JsonFunctions.stringList(synonyms.getAsJsonObject(), "xrefs"))
                    .build());
        }

        return doidSynonymsList;
    }

    @NotNull
    private static List<DoidBasicPropertyValues> extractBasicProportiesValues(@NotNull JsonArray arrayBasicPropertyValues) {
        List<DoidBasicPropertyValues> doidBasicPropertyValuesList = Lists.newArrayList();
        for (JsonElement basicProperty : arrayBasicPropertyValues) {
            doidBasicPropertyValuesList.add(ImmutableDoidBasicPropertyValues.builder()
                    .pred(JsonFunctions.toString(basicProperty.getAsJsonObject(), "pred"))
                    .val(JsonFunctions.toString(basicProperty.getAsJsonObject(), "val"))
                    .build());
        }

        return doidBasicPropertyValuesList;

    }

    @NotNull
    private static DoidMetadata extractDoidMetadata(@NotNull JsonObject metaDataJsonObject) {
        ImmutableDoidMetadata.Builder doidMetadata = ImmutableDoidMetadata.builder();
        List<DoidXref> xrefValList = Lists.newArrayList();

        JsonArray arrayXref = metaDataJsonObject.getAsJsonArray("xrefs");
        if (arrayXref != null) {
            for (JsonElement xref : arrayXref) {
                xrefValList.add(ImmutableDoidXref.builder().val(JsonFunctions.toString(xref.getAsJsonObject(), "val")).build());
            }
        }
        doidMetadata.synonyms(metaDataJsonObject.getAsJsonArray("synonyms") == null
                ? null
                : extractDoidSynonyms(metaDataJsonObject.getAsJsonArray("synonyms")));
        doidMetadata.basicPropertyValues(metaDataJsonObject.getAsJsonArray("basicPropertyValues") == null
                ? null
                : extractBasicProportiesValues(metaDataJsonObject.getAsJsonArray("basicPropertyValues")));
        doidMetadata.doidDefinition(metaDataJsonObject.getAsJsonObject("definition") == null
                ? null
                : extractDoidDefinition(metaDataJsonObject.getAsJsonObject("definition")));
        doidMetadata.subset(JsonFunctions.optionalStringList(metaDataJsonObject, "subsets"));
        doidMetadata.xref(xrefValList);
        return doidMetadata.build();

    }

    @NotNull
    private static String extractDoidId(@NotNull String id) {
        return id.replace("http://purl.obolibrary.org/obo/DOID_", "");
    }

    @VisibleForTesting
    public static List<Doid> readDoidJsonFile(@NotNull String doidJsonFile) throws IOException {
        JsonParser parser = new JsonParser();
        JsonReader reader = new JsonReader(new FileReader(doidJsonFile));
        reader.setLenient(true);
        List<Doid> doids = Lists.newArrayList();
        while (reader.peek() != JsonToken.END_DOCUMENT) {
            JsonObject doidObject = parser.parse(reader).getAsJsonObject();

            JsonArray arrayGraphs = doidObject.getAsJsonObject().getAsJsonArray("graphs").getAsJsonArray();
            for (JsonElement graph : arrayGraphs) {
                JsonArray arrayNodes = graph.getAsJsonObject().getAsJsonArray("nodes");
                for (JsonElement nodes : arrayNodes) {
                    doids.add(ImmutableDoid.builder()
                            .id(JsonFunctions.toString(nodes.getAsJsonObject(), "id"))
                            .doid(extractDoidId(JsonFunctions.toString(nodes.getAsJsonObject(), "id")))
                            .doidMetadata(nodes.getAsJsonObject().getAsJsonObject("meta") == null
                                    ? null
                                    : extractDoidMetadata(nodes.getAsJsonObject().getAsJsonObject("meta")))
                            .type(nodes.getAsJsonObject().getAsJsonPrimitive("type") == null
                                    ? null
                                    : JsonFunctions.toString(nodes.getAsJsonObject(), "type"))
                            .doidTerm(nodes.getAsJsonObject().getAsJsonPrimitive("lbl") == null
                                    ? null
                                    : JsonFunctions.nullableString(nodes.getAsJsonObject(), "lbl"))
                            .build());
                }
            }
        }
        return doids;
    }

}
