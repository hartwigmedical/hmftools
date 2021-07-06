package com.hartwig.hmftools.patientdb.clinical.readers;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.lang.reflect.Type;
import java.util.ArrayList;
import java.util.List;

import com.google.common.reflect.TypeToken;
import com.google.gson.Gson;
import com.google.gson.GsonBuilder;
import com.hartwig.hmftools.patientdb.clinical.context.RunContext;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class RunsJsonReader {

    private static final Gson GSON = new GsonBuilder().create();

    @NotNull
    public static List<RunContext> extractRunContexts(@NotNull final File jsonFile) {
        try {
            List<RunContext> runContexts = GSON.fromJson(new FileReader(jsonFile), listOfRunContext());
            if (runContexts != null) {
                return runContexts;
            } else {
                throw new IllegalArgumentException(String.format("File [%s] did not contain run contexts", jsonFile.getPath()));
            }
        } catch (FileNotFoundException e) {
            throw new IllegalArgumentException(e);
        }
    }

    @SuppressWarnings("UnstableApiUsage")
    @NotNull
    private static Type listOfRunContext() {
        return new TypeToken<ArrayList<JsonRunContext>>() {
        }.getType();
    }

    static class JsonRunContext implements RunContext {
        @Nullable
        private String setName;
        @Nullable
        private String refSample;
        @Nullable
        private String tumorSample;
        @Nullable
        private String tumorBarcodeSample;

        public void setSetName(@NotNull String setName) {
            this.setName = setName;
        }

        public void setRefSample(@NotNull String refSample) {
            this.refSample = refSample;
        }

        public void setTumorSample(@NotNull String tumorSample) {
            this.tumorSample = tumorSample;
        }

        public void setTumorBarcodeSample(@NotNull String tumorBarcodeSample) {
            this.tumorBarcodeSample = tumorBarcodeSample;
        }

        @NotNull
        @Override
        public String runDirectory() {
            return "json_file";
        }

        @NotNull
        @Override
        public String setName() {
            return throwWhenNull(setName, "setName");
        }

        @NotNull
        @Override
        public String refSample() {
            return throwWhenNull(refSample, "refSample");
        }

        @NotNull
        @Override
        public String tumorSample() {
            return throwWhenNull(tumorSample, "tumorSample");
        }

        @NotNull
        @Override
        public String tumorBarcodeSample() {
            return throwWhenNull(tumorBarcodeSample, "tumorBarcodeSample");
        }

        @NotNull
        static String throwWhenNull(@Nullable String nullable, @NotNull String field) {
            if (nullable == null) {
                throw new IllegalStateException(String.format("Field [%s] was null in the input JSON", field));
            }
            return nullable;
        }
    }
}
