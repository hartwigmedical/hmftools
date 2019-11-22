package com.hartwig.hmftools.patientdb.context;

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
    private static final String BARCODE_TUMOR_SAMPLE_FIELD_P5 = "sampleId";

    private static final String SET_NAME_FIELD_P4 = "set_name";
    private static final String SET_NAME_FIELD_P5 = "runName";

    private static final String NO_TUMOR_SAMPLE = "NA";
    private static final String BARCODE_START = "FR";

    private static final Gson GSON = new GsonBuilder().create();

    private MetaDataResolver() {
    }

    @Nullable
    static RunContext fromMetaDataFile(@NotNull final String runDirectory) {
        File metaDataFileP4 = new File(runDirectory + File.separator + METADATA_FILE_P4);
        File metaDataFileP5 = new File(runDirectory + File.separator + METADATA_FILE_P5);

        if (metaDataFileP4.exists()) {
            try {
                return fromPv4MetaData(runDirectory, metaDataFileP4);
            } catch (FileNotFoundException exception) {
                LOGGER.warn("Could not find meta data file {} for run dir {}.", METADATA_FILE_P4, runDirectory);
                return null;
            }
        } else if (metaDataFileP5.exists()) {
            try {
                return fromPv5MetaData(runDirectory, metaDataFileP5);
            } catch (FileNotFoundException exception) {
                LOGGER.warn("Could not find meta data file {} for run dir {}.", METADATA_FILE_P5, runDirectory);
                return null;
            }
        } else {
            LOGGER.warn("ERROR no metadata file found for run dir {}.", runDirectory);
            return null;
        }
    }

    @Nullable
    private static RunContext fromPv4MetaData(@NotNull String runDirectory, @NotNull File pv4MetadataFile) throws FileNotFoundException {
        JsonObject json = GSON.fromJson(new FileReader(pv4MetadataFile), JsonObject.class);

        String refSample = fieldValue(json, REF_SAMPLE_FIELD_P4);
        String tumorSample = fieldValue(json, TUMOR_SAMPLE_FIELD_P4);
        String setName = fieldValue(json, SET_NAME_FIELD_P4);
        String tumorBarcodeSample = Strings.EMPTY;

        if (refSample == null) {
            LOGGER.warn("Could not find " + REF_SAMPLE_FIELD_P4 + " in metadata file!");
            return null;
        } else if (setName == null) {
            LOGGER.warn("Could not find " + SET_NAME_FIELD_P4 + " in metadata file!");
            return null;
        } else if (tumorSample == null) {
            LOGGER.warn("Could not find " + TUMOR_SAMPLE_FIELD_P4 + " in metadata file!");
            return null;
        }

        boolean containsFR = false;
        for (String setNamePart : setName.split("_")) {
            if (setNamePart.startsWith(BARCODE_START)) {
                containsFR = true;
                tumorBarcodeSample = setNamePart;
            }
        }
        if (!containsFR) {
            LOGGER.warn("No tumor barcode is known for set set '{}'", setName);
        }

        return new RunContextImpl(runDirectory, setName, refSample, tumorSample, tumorBarcodeSample);
    }

    @Nullable
    private static RunContext fromPv5MetaData(@NotNull String runDirectory, @NotNull File pv5MetadataFile) throws FileNotFoundException {
        JsonObject json = GSON.fromJson(new FileReader(pv5MetadataFile), JsonObject.class);

        String refSample = sampleIdP5(json, REF_SAMPLE_FIELD_P5);
        String tumorSample = sampleIdP5(json, TUMOR_SAMPLE_FIELD_P5);
        String tumorBarcodeSample = sampleBarcodeP5(json, TUMOR_SAMPLE_FIELD_P5);
        String setName = fieldValue(json, SET_NAME_FIELD_P5);

        if (refSample == null) {
            LOGGER.warn("Could not find " + REF_SAMPLE_FIELD_P5 + " in metadata file!");
            return null;
        } else if (tumorSample == null) {
            LOGGER.warn("Could not find " + TUMOR_SAMPLE_FIELD_P5 + " in metadata file!");
            return null;
        } else if (setName == null) {
            LOGGER.warn("Could not find " + SET_NAME_FIELD_P5 + " in metadata file!");
            return null;
        } else if (tumorBarcodeSample == null) {
            LOGGER.warn("Could not find " + BARCODE_TUMOR_SAMPLE_FIELD_P5 + " in metadata file!");
            return null;
        }
        return new RunContextImpl(runDirectory, setName, refSample, tumorSample, tumorBarcodeSample);
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

    @Nullable
    private static String sampleBarcodeP5(@NotNull final JsonObject object, @NotNull final String fieldName) {
        final JsonElement element = object.get(fieldName);
        JsonElement sampleBarcodeId = element.getAsJsonObject().get("sampleId");
        return sampleBarcodeId != null && !(sampleBarcodeId instanceof JsonNull) ? sampleBarcodeId.getAsString() : null;
    }
}
