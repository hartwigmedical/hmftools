package com.hartwig.hmftools.common.vicc;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;

import com.google.gson.JsonArray;
import com.google.gson.JsonObject;
import com.google.gson.JsonParser;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class ViccFactory {
    private static final Logger LOGGER = LogManager.getLogger(ViccFactory.class);

    private static final String ALL_JSON_FILE = "all.json";
    private static final String BRCA_JSON_FILE = "brca.json";
    private static final String CGI_JSON_FILE = "cgi.json";
    private static final String CIVIC_JSON_FILE = "civic.json";
    private static final String JAX_JSON_FILE = "jax.json";
    private static final String JAXTRIALS_JSON_FILE = "jax_trials.json";
    private static final String MOLECULARMATCH_JSON_FILE = "molecularmatch.json";
    private static final String MOLECULARMATCHTRIALS_JSON_FILE = "molecularmatch_trials.json";
    private static final String ONCOKB_JSON_FILE = "oncokb.json";
    private static final String PMKB_JSON_FILE = "pmkb.json";
    private static final String SAGE_JSON_FILE = "sage.json";

    private ViccFactory() {
    }

    public static void VICCDirectory(@NotNull final String viccDirectory) throws IOException {
        final String brcaJsonPath = viccDirectory + File.separator + BRCA_JSON_FILE;

        readingBRCAFile(brcaJsonPath);

    }

    public static void readingBRCAFile(@NotNull String brcaJsonPath) throws IOException {

      //  final Gson gson = VICCGsonAdaptor.buildSampleGson();
        final Object object = new JsonParser().parse(new FileReader(brcaJsonPath));
        JsonArray jsonObject = (JsonArray) object;
        for (Object dataJson : jsonObject ){
            JsonObject brca = (JsonObject) dataJson;
            brca.getAsJsonObject("brca").get("Variant_frequency_LOVD");
        }

    }

}
