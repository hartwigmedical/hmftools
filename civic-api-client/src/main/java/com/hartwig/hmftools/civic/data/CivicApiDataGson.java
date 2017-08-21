package com.hartwig.hmftools.civic.data;

import com.google.gson.Gson;
import com.google.gson.GsonBuilder;

public class CivicApiDataGson {
    public static Gson buildGson() {
        return new GsonBuilder().registerTypeAdapterFactory(new GsonAdaptersCivicDisease())
                .registerTypeAdapterFactory(new GsonAdaptersCivicDrug())
                .registerTypeAdapterFactory(new GsonAdaptersCivicEvidenceItem())
                .registerTypeAdapterFactory(new GsonAdaptersCivicEvidenceItemMetadata())
                .registerTypeAdapterFactory(new GsonAdaptersCivicGene())
                .registerTypeAdapterFactory(new GsonAdaptersCivicVariant())
                .registerTypeAdapterFactory(new GsonAdaptersCivicVariantCoordinates())
                .registerTypeAdapterFactory(new GsonAdaptersCivicVariantKey())
                .create();
    }
}
