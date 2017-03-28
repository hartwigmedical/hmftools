package com.hartwig.hmftools.healthchecker.context;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;

import com.google.gson.Gson;
import com.google.gson.GsonBuilder;
import com.google.gson.JsonObject;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

final class MetaDataResolver {

    private static final String METADATA_FILE = "metadata";

    private static final String REF_SAMPLE_FIELD = "ref_sample";
    private static final String TUMOR_SAMPLE_FIELD = "tumor_sample";
    private static final String SET_NAME_FIELD = "set_name";

    private static final String NO_TUMOR_SAMPLE = "NA";

    private static final Gson GSON = new GsonBuilder().create();

    private MetaDataResolver() {
    }

    @Nullable
    static RunContext fromMetaDataFile(@NotNull final String runDirectory) throws FileNotFoundException {
        final String metaDataFilePath = runDirectory + File.separator + METADATA_FILE;
        if (new File(metaDataFilePath).exists()) {
            final JsonObject json = GSON.fromJson(new FileReader(metaDataFilePath), JsonObject.class);
            final String refSample = json.get(REF_SAMPLE_FIELD).getAsString();
            final String tumorSample = json.get(TUMOR_SAMPLE_FIELD).getAsString();
            final String setName = json.get(SET_NAME_FIELD).getAsString();

            final boolean isSomaticRun = !tumorSample.equals(NO_TUMOR_SAMPLE);

            return new RunContextImpl(runDirectory, setName, refSample, isSomaticRun ? tumorSample : Strings.EMPTY,
                    isSomaticRun);
        }
        return null;
    }
}
