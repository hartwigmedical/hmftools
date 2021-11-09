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

    private static final Logger LOGGER = LogManager.getLogger(DiseaseOntology.class);

    @VisibleForTesting
    static final String ID_TO_READ = "http://purl.obolibrary.org/obo/doid.owl";

    private DiseaseOntology() {
    }

    @NotNull
    public static DoidEntry readDoidOwlEntryFromDoidJson(@NotNull String doidJsonFile) throws IOException {
        JsonParser parser = new JsonParser();
        JsonReader reader = new JsonReader(new FileReader(doidJsonFile));
        reader.setLenient(true);

        JsonObject doidObject = parser.parse(reader).getAsJsonObject();
        JsonDatamodelChecker doidObjectChecker = DoidDatamodelCheckerFactory.doidObjectChecker();
        doidObjectChecker.check(doidObject);

        ImmutableDoidEntry.Builder doidEntryBuilder = ImmutableDoidEntry.builder();
        JsonArray graphArray = doidObject.getAsJsonArray("graphs");
        for (JsonElement graphElement : graphArray) {
            JsonObject graphObject = graphElement.getAsJsonObject();
            JsonDatamodelChecker doidGraphsChecker = DoidDatamodelCheckerFactory.doidGraphsChecker();
            doidGraphsChecker.check(graphObject);

            String id = string(graphObject, "id");
            if (id.equals(ID_TO_READ)) {
                LOGGER.debug(" Reading DOID entry with ID '{}'", id);

                // Add data to doid entry
                doidEntryBuilder.id(id);
                doidEntryBuilder.nodes(extractNodes(graphObject.getAsJsonArray("nodes")));
                doidEntryBuilder.edges(extractEdges(graphObject.getAsJsonArray("edges")));
                doidEntryBuilder.meta(extractGraphMetaNode(graphObject.getAsJsonObject("meta")));
                doidEntryBuilder.logicalDefinitionAxioms(extractDoidLogicalDefinitionAxioms(graphObject.getAsJsonArray(
                        "logicalDefinitionAxioms")));

                // Below always seem to be empty string lists
                doidEntryBuilder.equivalentNodesSets(optionalStringList(graphObject, "equivalentNodesSets"));
                doidEntryBuilder.domainRangeAxioms(optionalStringList(graphObject, "domainRangeAxioms"));
                doidEntryBuilder.propertyChainAxioms(optionalStringList(graphObject, ("propertyChainAxioms")));
            }
        }

        if (reader.peek() != JsonToken.END_DOCUMENT) {
            LOGGER.warn("More data found in {} after reading main JSON object!", doidJsonFile);
        }

        return doidEntryBuilder.build();
    }

    @NotNull
    static String extractDoid(@NotNull String url) {
        return url.replace("http://purl.obolibrary.org/obo/DOID_", "");
    }

    @NotNull
    private static List<DoidNode> extractNodes(@NotNull JsonArray nodeArray) {
        List<DoidNode> nodes = Lists.newArrayList();

        for (JsonElement nodeElement : nodeArray) {
            JsonObject nodeObject = nodeElement.getAsJsonObject();
            JsonDatamodelChecker doidNodesChecker = DoidDatamodelCheckerFactory.doidNodeChecker();
            doidNodesChecker.check(nodeObject);

            String url = string(nodeObject, "id");
            nodes.add(ImmutableDoidNode.builder()
                    .doid(extractDoid(url))
                    .url(url)
                    .doidMetadata(extractDoidMetadata(optionalJsonObject(nodeObject, "meta")))
                    .type(optionalString(nodeObject, "type"))
                    .doidTerm(optionalString(nodeObject, "lbl"))
                    .build());
        }

        return nodes;
    }

    @NotNull
    private static List<DoidEdge> extractEdges(@NotNull JsonArray edgeArray) {
        List<DoidEdge> edges = Lists.newArrayList();

        for (JsonElement edgeElement : edgeArray) {
            JsonObject edgeObject = edgeElement.getAsJsonObject();
            JsonDatamodelChecker doidEdgeChecker = DoidDatamodelCheckerFactory.doidEdgeChecker();
            doidEdgeChecker.check(edgeObject);

            String predicate = string(edgeObject, "pred");
            String object = string(edgeObject, "obj");
            String subject = string(edgeObject, "sub");

            edges.add(ImmutableDoidEdge.builder().predicate(predicate).subject(subject).object(object).build());
        }

        return edges;
    }

    @NotNull
    private static DoidGraphMetaData extractGraphMetaNode(@Nullable JsonObject metaObject) {
        JsonDatamodelChecker doidGraphMetaDataChecker = DoidDatamodelCheckerFactory.doidGraphMetaDataChecker();
        doidGraphMetaDataChecker.check(metaObject);

        JsonArray xrefArray = metaObject.getAsJsonArray("xrefs");
        List<DoidXref> xrefValList = Lists.newArrayList();
        if (xrefArray != null) {
            JsonDatamodelChecker xrefChecker = DoidDatamodelCheckerFactory.doidMetadataXrefChecker();

            for (JsonElement xrefElement : xrefArray) {
                JsonObject xrefObject = xrefElement.getAsJsonObject();

                xrefChecker.check(xrefObject);
                xrefValList.add(ImmutableDoidXref.builder().val(string(xrefObject, "val")).build());
            }
        }

        ImmutableDoidGraphMetaData.Builder DoidGraphMetaDataBuilder = ImmutableDoidGraphMetaData.builder();
        DoidGraphMetaDataBuilder.basicPropertyValues(extractBasicPropertyValues(optionalJsonArray(metaObject, "basicPropertyValues")));
        DoidGraphMetaDataBuilder.subsets(optionalStringList(metaObject, "subsets"));
        DoidGraphMetaDataBuilder.xrefs(xrefValList);
        DoidGraphMetaDataBuilder.version(optionalString(metaObject, "version"));
        return DoidGraphMetaDataBuilder.build();
    }

    @Nullable
    private static List<DoidBasicPropertyValue> extractBasicPropertyValues(@Nullable JsonArray basicPropertyValueArray) {
        if (basicPropertyValueArray == null) {
            return null;
        }

        List<DoidBasicPropertyValue> basicPropertyValueList = Lists.newArrayList();

        for (JsonElement basicPropertyElement : basicPropertyValueArray) {
            JsonObject basicPropertyObject = basicPropertyElement.getAsJsonObject();
            JsonDatamodelChecker basicPropertyValuesChecker = DoidDatamodelCheckerFactory.doidBasicPropertyValueChecker();
            basicPropertyValuesChecker.check(basicPropertyObject);

            basicPropertyValueList.add(ImmutableDoidBasicPropertyValue.builder()
                    .pred(string(basicPropertyObject, "pred"))
                    .val(string(basicPropertyObject, "val"))
                    .build());
        }

        return basicPropertyValueList;
    }

    @Nullable
    private static List<DoidLogicalDefinitionAxioms> extractDoidLogicalDefinitionAxioms(@Nullable JsonArray logicalDefinitionAxiomArray) {
        if (logicalDefinitionAxiomArray == null) {
            return null;
        }

        List<DoidLogicalDefinitionAxioms> logicalDefinitionAxioms = Lists.newArrayList();

        for (JsonElement logicalDefinitionAxiomElement : logicalDefinitionAxiomArray) {
            JsonObject logicalDefinitionAxiomObject = logicalDefinitionAxiomElement.getAsJsonObject();

            JsonDatamodelChecker doidLogicalDefinitionAxiomsChecker = DoidDatamodelCheckerFactory.doidLogicalDefinitionAxiomChecker();
            doidLogicalDefinitionAxiomsChecker.check(logicalDefinitionAxiomObject);

            List<DoidRestriction> restrictionList = Lists.newArrayList();
            for (JsonElement restrictionElement : logicalDefinitionAxiomObject.getAsJsonArray("restrictions")) {
                JsonDatamodelChecker doidRestrictionChecker = DoidDatamodelCheckerFactory.doidRestrictionChecker();

                if (restrictionElement.isJsonObject()) {
                    JsonObject restrictionObject = restrictionElement.getAsJsonObject();
                    doidRestrictionChecker.check(restrictionObject.getAsJsonObject());

                    restrictionList.add(ImmutableDoidRestriction.builder()
                            .propertyId(string(restrictionObject, "propertyId"))
                            .fillerId(string(restrictionObject, "fillerId"))
                            .build());
                }
            }

            List<String> genusIdList = Lists.newArrayList();
            for (JsonElement genusIdElement : logicalDefinitionAxiomObject.getAsJsonArray("genusIds")) {
                genusIdList.add(genusIdElement.getAsString());
            }

            logicalDefinitionAxioms.add(ImmutableDoidLogicalDefinitionAxioms.builder()
                    .definedClassId(string(logicalDefinitionAxiomObject, "definedClassId"))
                    .genusIds(genusIdList)
                    .restrictions(restrictionList)
                    .build());
        }

        return logicalDefinitionAxioms;
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
            for (JsonElement xrefElement : xrefArray) {
                JsonObject xrefObject = xrefElement.getAsJsonObject();
                xrefChecker.check(xrefObject);

                xrefValList.add(ImmutableDoidXref.builder().val(string(xrefObject, "val")).build());
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

        List<DoidSynonym> synonymList = Lists.newArrayList();
        for (JsonElement synonymElement : synonymArray) {
            JsonObject synonymObject = synonymElement.getAsJsonObject();
            JsonDatamodelChecker doidSynonymChecker = DoidDatamodelCheckerFactory.doidSynonymChecker();
            doidSynonymChecker.check(synonymObject);

            synonymList.add(ImmutableDoidSynonym.builder()
                    .pred(string(synonymObject, "pred"))
                    .val(string(synonymObject, "val"))
                    .xrefs(stringList(synonymObject, "xrefs"))
                    .build());
        }

        return synonymList;
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
