package com.hartwig.hmftools.patientdb.clinical.readers;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.lang.reflect.Type;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import com.google.common.reflect.TypeToken;
import com.google.gson.Gson;
import com.google.gson.GsonBuilder;
import com.hartwig.hmftools.common.runcontext.RunContext;

import org.jetbrains.annotations.NotNull;

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
    private static Type listOfRunContext() {
        return new TypeToken<ArrayList<JsonRunContext>>() {
        }.getType();
    }

    static class JsonRunContext implements RunContext {
        private String setName;
        private String refSample;
        private String tumorSample;
        private String tumorBarcodeSample;

        public void setSetName(final String setName) {
            this.setName = setName;
        }

        public void setRefSample(final String refSample) {
            this.refSample = refSample;
        }

        public void setTumorSample(final String tumorSample) {
            this.tumorSample = tumorSample;
        }

        public void setTumorBarcodeSample(final String tumorBarcodeSample) {
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
            return setName;
        }

        @NotNull
        @Override
        public String refSample() {
            return refSample;
        }

        @NotNull
        @Override
        public String tumorSample() {
            return tumorSample;
        }

        @NotNull
        @Override
        public String tumorBarcodeSample() {
            return tumorBarcodeSample;
        }
    }
}
