package com.hartwig.hmftools.apiclients.diseaseontology.data;

import com.google.gson.Gson;
import com.google.gson.GsonBuilder;

public class DiseaseOntologyGson {
    public static Gson buildGson() {
        return new GsonBuilder().registerTypeAdapterFactory(new GsonAdaptersDiseaseOntologyMetadata()).create();
    }
}
