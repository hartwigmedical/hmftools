package com.hartwig.hmftools.common.lims;

import com.google.gson.Gson;
import com.google.gson.GsonBuilder;

public class LimsGsonAdapter {
    public static Gson buildGson() {
        return new GsonBuilder().registerTypeAdapterFactory(new GsonAdaptersLimsJsonData()).create();
    }
}
