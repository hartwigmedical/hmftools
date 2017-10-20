package com.hartwig.hmftools.common.lims;

import com.google.gson.Gson;
import com.google.gson.GsonBuilder;

import org.jetbrains.annotations.NotNull;

final class LimsGsonAdapter {

    private LimsGsonAdapter() {
    }

    @NotNull
    static Gson buildGson() {
        return new GsonBuilder().registerTypeAdapterFactory(new GsonAdaptersLimsJsonData()).create();
    }
}
