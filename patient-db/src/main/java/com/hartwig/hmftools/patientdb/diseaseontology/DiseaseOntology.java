package com.hartwig.hmftools.patientdb.diseaseontology;

import static com.hartwig.hmftools.common.utils.json.JsonFunctions.nullableString;
import static com.hartwig.hmftools.common.utils.json.JsonFunctions.optionalJsonArray;
import static com.hartwig.hmftools.common.utils.json.JsonFunctions.optionalJsonObject;
import static com.hartwig.hmftools.common.utils.json.JsonFunctions.optionalString;
import static com.hartwig.hmftools.common.utils.json.JsonFunctions.optionalStringList;
import static com.hartwig.hmftools.common.utils.json.JsonFunctions.string;
import static com.hartwig.hmftools.common.utils.json.JsonFunctions.stringList;

import java.io.FileReader;
import java.io.IOException;
import java.util.List;

import com.google.common.collect.Lists;
import com.google.gson.JsonArray;
import com.google.gson.JsonElement;
import com.google.gson.JsonObject;
import com.google.gson.JsonParser;
import com.google.gson.stream.JsonReader;
import com.google.gson.stream.JsonToken;
import com.hartwig.hmftools.common.utils.json.DatamodelChecker;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class DiseaseOntology {

    private DiseaseOntology() {
    }

    @NotNull
    public static List<Doid> readDoidJsonFile(@NotNull String doidJsonFile) throws IOException {
        JsonParser parser = new JsonParser();
        JsonReader reader = new JsonReader(new FileReader(doidJsonFile));
        reader.setLenient(true);
        List<Doid> doids = Lists.newArrayList();
        while (reader.peek() != JsonToken.END_DOCUMENT) {
            JsonObject doidObject = parser.parse(reader).getAsJsonObject();

            JsonArray graphsArray = doidObject.getAsJsonObject().getAsJsonArray("graphs").getAsJsonArray();
            for (JsonElement graph : graphsArray) {
                JsonArray nodeArray = graph.getAsJsonObject().getAsJsonArray("nodes");
                for (JsonElement node : nodeArray) {
                    JsonObject nodeObject = node.getAsJsonObject();
                    String id = string(nodeObject, "id");
                    doids.add(ImmutableDoid.builder()
                            .id(id)
                            .doid(extractDoidId(id))
                            .doidMetadata(extractDoidMetadata(optionalJsonObject(nodeObject, "meta")))
                            .type(optionalString(nodeObject, "type"))
                            .doidTerm(optionalString(nodeObject, "lbl"))
                            .build());
                }
            }
        }
        return doids;
    }

    @NotNull
    private static String extractDoidId(@NotNull String id) {
        return id.replace("http://purl.obolibrary.org/obo/DOID_", "");
    }

    @Nullable
    private static DoidMetadata extractDoidMetadata(@Nullable JsonObject metaDataJsonObject) {
        if (metaDataJsonObject == null) {
            return null;
        }

        DatamodelChecker doidMetadataChecker = DoidDatamodelCheckerFactory.doidMetadataChecker();
        doidMetadataChecker.check(metaDataJsonObject);

        JsonArray arrayXref = metaDataJsonObject.getAsJsonArray("xrefs");
        List<DoidXref> xrefValList = Lists.newArrayList();
        if (arrayXref != null) {
            for (JsonElement xref : arrayXref) {
                xrefValList.add(ImmutableDoidXref.builder().val(string(xref.getAsJsonObject(), "val")).build());
            }
        }

        ImmutableDoidMetadata.Builder doidMetadataBuilder = ImmutableDoidMetadata.builder();
        doidMetadataBuilder.synonyms(extractDoidSynonyms(optionalJsonArray(metaDataJsonObject, "synonyms")));
        doidMetadataBuilder.basicPropertyValues(extractBasicPropertyValues(optionalJsonArray(metaDataJsonObject, "basicPropertyValues")));
        doidMetadataBuilder.doidDefinition(extractDoidDefinition(optionalJsonObject(metaDataJsonObject, "definition")));
        doidMetadataBuilder.subset(optionalStringList(metaDataJsonObject, "subsets"));
        doidMetadataBuilder.xref(xrefValList);
        return doidMetadataBuilder.build();
    }

    @Nullable
    private static List<DoidSynonym> extractDoidSynonyms(@Nullable JsonArray synonymArray) {
        if (synonymArray == null) {
            return null;
        }

        List<DoidSynonym> doidSynonymList = Lists.newArrayList();
        for (JsonElement synonym : synonymArray) {
            JsonObject synonymObject = synonym.getAsJsonObject();
            doidSynonymList.add(ImmutableDoidSynonym.builder()
                    .pred(string(synonymObject, "pred"))
                    .val(string(synonymObject, "val"))
                    .xrefs(stringList(synonymObject, "xrefs"))
                    .build());
        }

        return doidSynonymList;
    }

    @NotNull
    private static List<DoidBasicPropertyValue> extractBasicPropertyValues(@NotNull JsonArray basicPropertyValueArray) {
        List<DoidBasicPropertyValue> doidBasicPropertyValueList = Lists.newArrayList();
        DatamodelChecker doidBasicPropertyValuesChecker = DoidDatamodelCheckerFactory.doidBasicPropertyValuesChecker();

        for (JsonElement basicProperty : basicPropertyValueArray) {
            JsonObject basicPropertyObject = basicProperty.getAsJsonObject();
            doidBasicPropertyValuesChecker.check(basicPropertyObject);

            doidBasicPropertyValueList.add(ImmutableDoidBasicPropertyValue.builder()
                    .pred(string(basicPropertyObject, "pred"))
                    .val(string(basicPropertyObject, "val"))
                    .build());
        }

        return doidBasicPropertyValueList;
    }

    @Nullable
    private static DoidDefinition extractDoidDefinition(@Nullable JsonObject definitionObject) {
        if (definitionObject == null) {
            return null;
        }

        ImmutableDoidDefinition.Builder doidDefinitionBuilder = ImmutableDoidDefinition.builder();
        doidDefinitionBuilder.definitionVal(nullableString(definitionObject, "val"));
        doidDefinitionBuilder.definitionXref(stringList(definitionObject, "xrefs"));
        return doidDefinitionBuilder.build();
    }
}
