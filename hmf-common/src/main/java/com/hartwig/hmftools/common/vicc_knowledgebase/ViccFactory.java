package com.hartwig.hmftools.common.vicc_knowledgebase;

import java.io.File;
import java.io.IOException;
import java.time.LocalDate;
import java.util.Map;

import com.google.common.collect.Maps;

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

    @NotNull
    public static VICC VICCDirectory(@NotNull final String viccDirectory) throws IOException {
        final String brcaJsonPath = viccDirectory + File.separator + BRCA_JSON_FILE;

        Map<String, FileBRCA> BRCA = readingBRCAFile(brcaJsonPath);


        return new VICC(BRCA);
    }

    @NotNull
    public static Map<String, FileBRCA> readingBRCAFile(@NotNull String brcaJsonPath) {
        return Maps.newHashMap();
    }
}
