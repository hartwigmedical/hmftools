package com.hartwig.hmftools.common.vicc;

import static org.junit.Assert.assertNotNull;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;

import com.google.gson.JsonObject;
import com.google.gson.JsonParser;
import com.google.gson.stream.JsonReader;
import com.google.gson.stream.JsonToken;

import org.junit.Ignore;
import org.junit.Test;

public class ViccFactoryTest {

    @Test
    @Ignore
    public void convertAll() throws IOException {
        final String baseDir = System.getProperty("user.home") + File.separator + "hmf" + File.separator + "tmp";
        final String inputFile = baseDir + File.separator + "all.json";
        final String outputCsvFileName = baseDir + File.separator + "all.csv";

        ViccFactory.extractAllFileSpecificFields(inputFile, outputCsvFileName);
    }

    @Test
    @Ignore
    public void testMe() throws IOException {
        final String baseDir = System.getProperty("user.home") + File.separator + "hmf" + File.separator + "tmp";
        final String inputFile = baseDir + File.separator + "all.json";

        JsonParser parser = new JsonParser();
        JsonReader reader = new JsonReader(new FileReader(inputFile));
        reader.setLenient(true);

        while (reader.peek() != JsonToken.END_DOCUMENT) {
            JsonObject object = parser.parse(reader).getAsJsonObject();
            if (object.getAsJsonObject("civic") == null) {
                int x = 1;
            }
            assertNotNull(object);
        }
        reader.close();
    }
}