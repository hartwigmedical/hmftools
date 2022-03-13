package com.hartwig.hmftools.common.doid;

import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.utils.json.JsonDatamodelChecker;

import org.jetbrains.annotations.NotNull;

final class DoidDatamodelCheckerFactory {

    private DoidDatamodelCheckerFactory() {
    }

    @NotNull
    static JsonDatamodelChecker doidObjectChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("graphs", true);

        return new JsonDatamodelChecker("DoidObject", map);
    }

    @NotNull
    static JsonDatamodelChecker doidGraphsChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("nodes", true);
        map.put("edges", true);
        map.put("id", true);
        map.put("meta", true);
        map.put("equivalentNodesSets", true);
        map.put("logicalDefinitionAxioms", true);
        map.put("domainRangeAxioms", true);
        map.put("propertyChainAxioms", true);

        return new JsonDatamodelChecker("DoidGraphs", map);
    }

    @NotNull
    static JsonDatamodelChecker doidNodeChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("type", false);
        map.put("lbl", false);
        map.put("id", true);
        map.put("meta", false);

        return new JsonDatamodelChecker("DoidNode", map);
    }

    @NotNull
    static JsonDatamodelChecker doidEdgeChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("sub", true);
        map.put("pred", true);
        map.put("obj", true);

        return new JsonDatamodelChecker("DoidEdge", map);
    }

    @NotNull
    static JsonDatamodelChecker doidGraphMetaDataChecker() {
        Map<String, Boolean> map = Maps.newHashMap();

        map.put("xrefs", true);
        map.put("basicPropertyValues", true);
        map.put("version", false);
        map.put("subsets", true);

        return new JsonDatamodelChecker("DoidGraphMetaData", map);
    }

    @NotNull
    static JsonDatamodelChecker doidLogicalDefinitionAxiomChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("definedClassId", true);
        map.put("genusIds", true);
        map.put("restrictions", true);

        return new JsonDatamodelChecker("DoidLogicalDefinitionAxiom", map);
    }

    @NotNull
    static JsonDatamodelChecker doidRestrictionChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("propertyId", true);
        map.put("fillerId", true);

        return new JsonDatamodelChecker("DoidRestriction", map);
    }

    @NotNull
    static JsonDatamodelChecker doidSynonymChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("pred", true);
        map.put("val", true);
        map.put("xrefs", true);
        return new JsonDatamodelChecker("DoidSynonym", map);
    }

    @NotNull
    static JsonDatamodelChecker doidDefinitionChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("xrefs", true);
        map.put("val", true);
        return new JsonDatamodelChecker("DoidDefinition", map);
    }

    @NotNull
    static JsonDatamodelChecker doidBasicPropertyValueChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("pred", true);
        map.put("val", true);
        return new JsonDatamodelChecker("DoidBasicPropertyValue", map);
    }

    @NotNull
    static JsonDatamodelChecker doidMetadataXrefChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("val", true);

        return new JsonDatamodelChecker("DoidMetadataXref", map);
    }

    @NotNull
    static JsonDatamodelChecker doidMetadataChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("xrefs", false);
        map.put("synonyms", false);
        map.put("basicPropertyValues", false);
        map.put("definition", false);
        map.put("subsets", false);
        map.put("deprecated", false);
        map.put("comments", false);

        return new JsonDatamodelChecker("DoidMetadata", map);
    }
}
