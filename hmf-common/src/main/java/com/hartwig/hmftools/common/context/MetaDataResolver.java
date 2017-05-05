package com.hartwig.hmftools.common.context;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;

import com.google.gson.Gson;
import com.google.gson.GsonBuilder;
import com.google.gson.JsonElement;
import com.google.gson.JsonObject;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

final class MetaDataResolver {

    private static final Logger LOGGER = LogManager.getLogger(MetaDataResolver.class);
    private static final String METADATA_FILE = "metadata";

    private static final String REF_SAMPLE_FIELD = "ref_sample";
    private static final String TUMOR_SAMPLE_FIELD = "tumor_sample";
    private static final String SET_NAME_FIELD = "set_name";

    private static final String NO_TUMOR_SAMPLE = "NA";

    private static final Gson GSON = new GsonBuilder().create();

    private MetaDataResolver() {
    }

    @Nullable
    static RunContext fromMetaDataFile(@NotNull final String runDirectory) {
        final String metaDataFilePath = runDirectory + File.separator + METADATA_FILE;
        final JsonObject json;
        try {
            json = GSON.fromJson(new FileReader(metaDataFilePath), JsonObject.class);
        } catch (FileNotFoundException exception) {
            return null;
        }

        final String refSample = fieldValue(json, REF_SAMPLE_FIELD);
        final String setName = fieldValue(json, SET_NAME_FIELD);

        if (refSample == null) {
            LOGGER.warn("Could not find " + REF_SAMPLE_FIELD + " in metadata file!");
            return null;
        } else if (setName == null) {
            LOGGER.warn("Could not find " + SET_NAME_FIELD + " in metadata file!");
            return null;
        } else {
            final String tumorSample = fieldValue(json, TUMOR_SAMPLE_FIELD);
            final boolean isSomaticRun = tumorSample != null && !tumorSample.equals(NO_TUMOR_SAMPLE);

            return new RunContextImpl(runDirectory, setName, refSample, isSomaticRun ? tumorSample : Strings.EMPTY,
                    isSomaticRun);
        }
    }

    @Nullable
    private static String fieldValue(@NotNull final JsonObject object, @NotNull final String fieldName) {
        final JsonElement element = object.get(fieldName);
        return element != null ? element.getAsString() : null;
    }
}
