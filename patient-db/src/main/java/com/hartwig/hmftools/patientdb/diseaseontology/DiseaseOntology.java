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

public class DiseaseOntology {

    private static final Logger LOGGER = LogManager.getLogger(DiseaseOntology.class);

    @VisibleForTesting
    public static List<Doid> readDoidJsonFile(@NotNull String doidJsonFile)throws IOException {
        List<Doid> doids = Lists.newArrayList();
        JsonParser parser = new JsonParser();
        JsonReader reader = new JsonReader(new FileReader(doidJsonFile));
        reader.setLenient(true);

        while (reader.peek() != JsonToken.END_DOCUMENT) {
            JsonObject doidObject = parser.parse(reader).getAsJsonObject();
            LOGGER.info("doidObject {} ",  doidObject);

            JsonArray arrayGraphs = doidObject.getAsJsonObject().getAsJsonArray("graphs").getAsJsonArray();
            for (JsonElement graph : arrayGraphs) {
                JsonArray arrayNodes = graph.getAsJsonObject().getAsJsonArray("nodes");
                for (JsonElement nodes: arrayNodes) {
                    LOGGER.info(nodes.getAsJsonObject().getAsJsonPrimitive("id"));
                    LOGGER.info(nodes.getAsJsonObject().getAsJsonObject("meta"));
                    LOGGER.info(nodes.getAsJsonObject().getAsJsonObject("meta").getAsJsonObject("definition").getAsJsonPrimitive("val"));
                    LOGGER.info(nodes.getAsJsonObject().getAsJsonObject("meta").getAsJsonObject("definition").getAsJsonArray("xrefs"));
                    LOGGER.info(nodes.getAsJsonObject().getAsJsonObject("meta").getAsJsonArray("subsets"));
                    JsonArray arrayXref = nodes.getAsJsonObject().getAsJsonObject("meta").getAsJsonArray("xrefs");
                    if (arrayXref != null) {
                        for (JsonElement xref: arrayXref) {
                            LOGGER.info(xref.getAsJsonObject().getAsJsonPrimitive("val"));
                        }
                    }
                    LOGGER.info(nodes.getAsJsonObject().getAsJsonObject("meta").getAsJsonArray("synonyms"));
                    JsonArray arraySynonyms = nodes.getAsJsonObject().getAsJsonObject("meta").getAsJsonArray("synonyms");
                    for (JsonElement synonyms: arraySynonyms) {
                        LOGGER.info(synonyms.getAsJsonObject().getAsJsonPrimitive("pred"));
                        LOGGER.info(synonyms.getAsJsonObject().getAsJsonPrimitive("val"));
                        LOGGER.info(synonyms.getAsJsonObject().getAsJsonArray("xrefs"));
                    }
                    LOGGER.info(nodes.getAsJsonObject().getAsJsonObject("meta").getAsJsonArray("basicPropertyValues"));
                    JsonArray arrayBasicPropertyValues = nodes.getAsJsonObject().getAsJsonObject("meta").getAsJsonArray("basicPropertyValues");
                    for (JsonElement basicPropertyValues: arrayBasicPropertyValues) {
                        LOGGER.info(basicPropertyValues.getAsJsonObject().getAsJsonPrimitive("pred"));
                        LOGGER.info(basicPropertyValues.getAsJsonObject().getAsJsonPrimitive("val"));
                    }


                    LOGGER.info(nodes.getAsJsonObject().getAsJsonPrimitive("type"));
                    LOGGER.info(nodes.getAsJsonObject().getAsJsonPrimitive("lbl"));

                }
            }

            //doids.add(doidObject);
        }
        return doids;
    }

}
