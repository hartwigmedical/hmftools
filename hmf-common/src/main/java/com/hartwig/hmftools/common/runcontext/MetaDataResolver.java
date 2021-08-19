package com.hartwig.hmftools.common.runcontext;

import static com.hartwig.hmftools.common.utils.json.JsonFunctions.nullableString;
import static com.hartwig.hmftools.common.utils.json.JsonFunctions.optionalNullableString;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;

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

public final class MetaDataResolver {

    private static final Logger LOGGER = LogManager.getLogger(MetaDataResolver.class);

    private static final String METADATA_FILE_P4 = "metadata";
    private static final String METADATA_FILE_P5 = "metadata.json";

    private static final String REF_SAMPLE_ID_FIELD_P4 = "ref_sample";
    private static final String TUMOR_SAMPLE_ID_FIELD_P4 = "tumor_sample";
    private static final String SET_NAME_FIELD_P4 = "set_name";

    private static final String REF_SAMPLE_OBJECT_P5 = "reference";
    private static final String TUMOR_SAMPLE_OBJECT_P5 = "tumor";
    private static final String SET_NAME_FIELD_P5 = "runName";
    private static final String SET_NAME_FIELD_POST_P5_15 = "set";

    private static final String BARCODE_START = "FR";
    private static final String BARCODE_START_OLD = "HMF";

    private static final Gson GSON = new GsonBuilder().create();

    private MetaDataResolver() {
    }

    @Nullable
    public static RunContext fromMetaDataFile(@NotNull String runDirectory, @NotNull String pipelineVersionFile) throws IOException {
        File metaDataFileP4 = new File(runDirectory + File.separator + METADATA_FILE_P4);
        File metaDataFileP5 = new File(runDirectory + File.separator + METADATA_FILE_P5);

        if (metaDataFileP4.exists()) {
            try {
                return fromPv4MetaData(runDirectory, metaDataFileP4);
            } catch (FileNotFoundException exception) {
                LOGGER.warn("Could not find meta data file '{}' for run dir '{}'.", METADATA_FILE_P4, runDirectory);
                return null;
            }
        } else if (metaDataFileP5.exists()) {
            try {
                return fromPv5MetaData(runDirectory, metaDataFileP5, new File(runDirectory + File.separator + pipelineVersionFile));
            } catch (FileNotFoundException exception) {
                LOGGER.warn("Could not find meta data file '{}' for run dir '{}'.", METADATA_FILE_P5, runDirectory);
                return null;
            }
        } else {
            LOGGER.warn("No metadata file found for run dir '{}'.", runDirectory);
            return null;
        }
    }

    @Nullable
    public static String majorDotMinorVersion(@NotNull File pipelineVersionFile) throws IOException {
        String version = readPipelineVersion(pipelineVersionFile);
        if (version == null) {
            return null;
        } else {
            String[] versionSplit = version.split("\\.");
            return versionSplit[0] + "." + versionSplit[1];
        }
    }

    @Nullable
    private static RunContext fromPv4MetaData(@NotNull String runDirectory, @NotNull File pv4MetadataFile) throws FileNotFoundException {
        JsonObject json = GSON.fromJson(new FileReader(pv4MetadataFile), JsonObject.class);

        String refSample = optionalNullableString(json, REF_SAMPLE_ID_FIELD_P4);
        String tumorSample = optionalNullableString(json, TUMOR_SAMPLE_ID_FIELD_P4);
        String setName = optionalNullableString(json, SET_NAME_FIELD_P4);

        if (refSample == null) {
            LOGGER.warn("Could not find '{}' in metadata file!", REF_SAMPLE_ID_FIELD_P4);
            return null;
        } else if (tumorSample == null) {
            LOGGER.warn("Could not find '{}' in metadata file!", TUMOR_SAMPLE_ID_FIELD_P4);
            return null;
        } else if (setName == null) {
            LOGGER.warn("Could not find '{}' in metadata file!", SET_NAME_FIELD_P4);
            return null;
        }

        String tumorBarcodeSample = Strings.EMPTY;
        // Always take the final (second) barcode of setName (assume this is the tumor barcode)
        boolean containsBarcode = false;
        for (String setNamePart : setName.split("_")) {
            if (setNamePart.startsWith(BARCODE_START) || setNamePart.startsWith(BARCODE_START_OLD)) {
                containsBarcode = true;
                tumorBarcodeSample = setNamePart;
            }
        }

        if (!containsBarcode) {
            LOGGER.warn("No tumor barcode could be derived from set name for '{}'", setName);
        }

        return new RunContextImpl(runDirectory, setName, refSample, tumorSample, tumorBarcodeSample);
    }

    @Nullable
    private static RunContext fromPv5MetaData(@NotNull String runDirectory, @NotNull File pv5MetadataFile,
            @NotNull File pipelineVersionFile) throws IOException {
        JsonObject json = GSON.fromJson(new FileReader(pv5MetadataFile), JsonObject.class);

        String pipelineVersion = readPipelineVersion(pipelineVersionFile);

        String refSample = sampleIdP5(json, REF_SAMPLE_OBJECT_P5);
        String tumorSample = sampleIdP5(json, TUMOR_SAMPLE_OBJECT_P5);
        String tumorBarcodeSample = sampleBarcodeP5(json, TUMOR_SAMPLE_OBJECT_P5, pipelineVersion);

        String setName;
        if (isPostP5dot15(pipelineVersion)) {
            setName = nullableString(json, SET_NAME_FIELD_POST_P5_15);
        } else {
            // this is pre 5.15 pipelines
            setName = nullableString(json, SET_NAME_FIELD_P5);
        }

        if (refSample == null) {
            LOGGER.warn("Could not find ref sample id in metadata object '{}'!", REF_SAMPLE_OBJECT_P5);
            return null;
        } else if (tumorSample == null) {
            LOGGER.warn("Could not find tumor sample id in metadata object '{}'!", TUMOR_SAMPLE_OBJECT_P5);
            return null;
        } else if (tumorBarcodeSample == null) {
            LOGGER.warn("Could not find tumor sample barcode in metadata object '{}'!", TUMOR_SAMPLE_OBJECT_P5);
            return null;
        } else if (setName == null) {
            LOGGER.warn("Could not find '{}' in metadata file!", SET_NAME_FIELD_P5);
            return null;
        }

        return new RunContextImpl(runDirectory, setName, refSample, tumorSample, tumorBarcodeSample);
    }

    @Nullable
    private static String sampleIdP5(@NotNull JsonObject metadata, @NotNull String objectName) {
        JsonObject object = metadata.getAsJsonObject(objectName);
        if (object == null) {
            return null;
        }
        JsonElement sampleId = object.get("sampleName");
        return sampleId != null && !(sampleId instanceof JsonNull) ? sampleId.getAsString() : null;
    }

    @Nullable
    private static String sampleBarcodeP5(@NotNull JsonObject metadata, @NotNull String objectName, @NotNull String pipelineVersion) {
        JsonObject object = metadata.getAsJsonObject(objectName);
        if (object == null) {
            return null;
        }

        if (isPostP5dot15(pipelineVersion)) {
            return optionalNullableString(object, "barcode");
        } else {
            // this is pre 5.15 pipelines
            return optionalNullableString(object, "sampleId");
        }
    }

    @Nullable
    public static String readPipelineVersion(@NotNull File pipelineVersionFile) throws IOException {
        List<String> lines = Files.readAllLines(pipelineVersionFile.toPath());
        if (lines.isEmpty()) {
            throw new IOException("Pipeline version file seems empty on " + pipelineVersionFile.getPath());
        } else {
            if (lines.size() == 1) {
                return lines.get(0);
            } else {
                LOGGER.warn("Too many lines in pipeline version file {}!", pipelineVersionFile);
                return null;
            }
        }
    }

    private static boolean isPostP5dot15(@NotNull String pipelineVersion) {
        String minorVersion = pipelineVersion.substring(2, 4);
        return minorVersion.matches("[0-9]+") && Integer.parseInt(minorVersion) >= 15;
    }
}
