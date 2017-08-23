package com.hartwig.hmftools.civic.data;

import com.google.gson.Gson;
import com.google.gson.GsonBuilder;

public class CivicApiDataGson {
    public static Gson buildGson() {
        return new GsonBuilder().registerTypeAdapterFactory(new GsonAdaptersCivicApiMetadata())
                .registerTypeAdapterFactory(new GsonAdaptersCivicDisease())
                .registerTypeAdapterFactory(new GsonAdaptersCivicDrug())
                .registerTypeAdapterFactory(new GsonAdaptersCivicEvidenceItem())
                .registerTypeAdapterFactory(new GsonAdaptersCivicEvidenceItemMetadata())
                .registerTypeAdapterFactory(new GsonAdaptersCivicEvidenceSource())
                .registerTypeAdapterFactory(new GsonAdaptersCivicGene())
                .registerTypeAdapterFactory(new GsonAdaptersCivicIndexResult())
                .registerTypeAdapterFactory(new GsonAdaptersCivicVariant())
                .registerTypeAdapterFactory(new GsonAdaptersCivicVariantCoordinates())
                .registerTypeAdapterFactory(new GsonAdaptersCivicVariantMetadata())
                .registerTypeAdapterFactory(new GsonAdaptersCivicVariantType())
                .create();
    }
}
