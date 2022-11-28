package com.hartwig.hmftools.patientdb.clinical.lims;

import com.google.gson.Gson;
import com.google.gson.GsonBuilder;

import org.jetbrains.annotations.NotNull;

final class LimsGsonAdapter {

    private LimsGsonAdapter() {
    }

    @NotNull
    static Gson buildSampleGson() {
        return new GsonBuilder().registerTypeAdapterFactory(new GsonAdaptersLimsJsonSampleData()).create();
    }

    @NotNull
    static Gson buildSubmissionGson() {
        return new GsonBuilder().registerTypeAdapterFactory(new GsonAdaptersLimsJsonSubmissionData()).create();
    }
}
