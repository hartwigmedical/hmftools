package com.hartwig.hmftools.patientdb.diseaseontology;

import java.io.FileReader;
import java.io.IOException;
import java.util.List;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
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
    public static void readDoidJsonFile(@NotNull String doidJsonFile)throws IOException {
        List<Doid> doids = Lists.newArrayList();
        JsonParser parser = new JsonParser();
        JsonReader reader = new JsonReader(new FileReader(doidJsonFile));
        reader.setLenient(true);

        while (reader.peek() != JsonToken.END_DOCUMENT) {
            JsonObject doidObject = parser.parse(reader).getAsJsonObject();
            LOGGER.info(doidObject);
            //doids.add(doidObject);
        }
    }

}
