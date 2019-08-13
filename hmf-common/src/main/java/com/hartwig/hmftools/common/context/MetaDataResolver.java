package com.hartwig.hmftools.common.context;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;

import com.google.gson.Gson;
import com.google.gson.GsonBuilder;
import com.google.gson.JsonElement;
import com.google.gson.JsonNull;
import com.google.gson.JsonObject;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

final class MetaDataResolver {

    private static final Logger LOGGER = LogManager.getLogger(MetaDataResolver.class);
    private static final String METADATA_FILE_P4 = "metadata";
    private static final String METADATA_FILE_P5 = "metadata.json";

    private static final String REF_SAMPLE_FIELD_P4 = "ref_sample";
    private static final String TUMOR_SAMPLE_FIELD_P4 = "tumor_sample";

    private static final String REF_SAMPLE_FIELD_P5 = "reference";
    private static final String TUMOR_SAMPLE_FIELD_P5 = "tumor";

    private static final String SET_NAME_FIELD_P4 = "set_name";
    private static final String SET_NAME_FIELD_P5 = "runName";

    private static final String NO_TUMOR_SAMPLE = "NA";

    private static final Gson GSON = new GsonBuilder().create();

    private MetaDataResolver() {
    }

    @Nullable
    static RunContext fromMetaDataFile(@NotNull final String runDirectory) {
        final String metaDataFilePathP4 = runDirectory + File.separator + METADATA_FILE_P4;
        final String metaDataFilePathP5 = runDirectory + File.separator + METADATA_FILE_P5;
        final JsonObject json;
        File fileP4 = new File(metaDataFilePathP4);
        File fileP5 = new File(metaDataFilePathP5);
        String refSample4 = Strings.EMPTY;
        String refSample5 = Strings.EMPTY;
        String setName4 = Strings.EMPTY;
        String setName5 = Strings.EMPTY;
        String tumorSample4 = Strings.EMPTY;
        String tumorSample5 = Strings.EMPTY;
        String metadata = Strings.EMPTY;

        if (fileP4.exists()) {
            try {
                json = GSON.fromJson(new FileReader(metaDataFilePathP4), JsonObject.class);
                refSample4 = fieldValue(json, REF_SAMPLE_FIELD_P4);
                setName4 = fieldValue(json, SET_NAME_FIELD_P4);
                tumorSample4 = fieldValue(json, TUMOR_SAMPLE_FIELD_P4);
            } catch (FileNotFoundException exception) {
                LOGGER.warn("Could not find meta data file: " + metaDataFilePathP4);
                return null;
            }
        } else if (fileP5.exists()) {
            try {
                json = GSON.fromJson(new FileReader(metaDataFilePathP5), JsonObject.class);
                refSample5 = sampleIdP5(json, REF_SAMPLE_FIELD_P5);
                setName5 = fieldValue(json, SET_NAME_FIELD_P5);
                tumorSample5 = sampleIdP5(json, TUMOR_SAMPLE_FIELD_P5);
            } catch (FileNotFoundException exception) {
                LOGGER.warn("Could not find meta data file: " + metaDataFilePathP5);
                return null;
            }
        } else {
            metadata = null;
            LOGGER.info("ERROR no metadata file");
        }

        if (setName4 == null && fileP4.exists()) {
            LOGGER.warn("Could not find " + SET_NAME_FIELD_P4 + " in metadata file!");
            return null;
        } else if (setName5 == null && fileP5.exists()) {
            LOGGER.warn("Could not find " + SET_NAME_FIELD_P5 + " in metadata file!");
            return null;
        }

        if (refSample4 == null && fileP4.exists()) {
            LOGGER.warn("Could not find " + REF_SAMPLE_FIELD_P4 + " in metadata file!");
            return null;
        } else if (refSample5 == null && fileP5.exists()) {
            LOGGER.warn("Could not find " + REF_SAMPLE_FIELD_P5 + " in metadata file!");
            return null;
        }

        if (metadata == null) {
            LOGGER.info("No metadata present");
            return null;
        } else {
            String tumorSample = Strings.EMPTY;
            String setName = Strings.EMPTY;
            String refSample = Strings.EMPTY;

            if (tumorSample4 != null && fileP4.exists()) {
                tumorSample = tumorSample4;
            }
            if (tumorSample5 != null && fileP5.exists()) {
                tumorSample = tumorSample5;
            }

            if (refSample4 != null && fileP4.exists()) {
                refSample = refSample4;
            }
            if (refSample5 != null && fileP5.exists()) {
                refSample = refSample5;
            }

            if (setName4 != null && fileP4.exists()) {
                setName = setName4;
            }
            if( setName5 != null && fileP5.exists()){
                setName = setName5;
            }
            final boolean isSomaticRun = !tumorSample.equals(NO_TUMOR_SAMPLE);
            return new RunContextImpl(runDirectory, setName, refSample, isSomaticRun ? tumorSample : Strings.EMPTY, isSomaticRun);
        }
    }

    @Nullable
    private static String fieldValue(@NotNull final JsonObject object, @NotNull final String fieldName) {
        final JsonElement element = object.get(fieldName);
        return element != null && !(element instanceof JsonNull) ? element.getAsString() : null;
    }

    @Nullable
    private static String sampleIdP5(@NotNull final JsonObject object, @NotNull final String fieldName) {
        final JsonElement element = object.get(fieldName);
        JsonElement sampleId = element.getAsJsonObject().get("sampleName");
        return sampleId != null && !(sampleId instanceof JsonNull) ? sampleId.getAsString() : null;
    }
}
