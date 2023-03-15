package com.hartwig.hmftools.datamodel;

import com.google.gson.Gson;
import com.google.gson.GsonBuilder;
import com.google.gson.TypeAdapterFactory;
import com.hartwig.hmftools.datamodel.orange.OrangeRecord;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ServiceLoader;

public class OrangeJson {

    @Nullable
    private static OrangeJson instance;

    private final Gson gson;

    @NotNull
    public static OrangeJson getInstance() {
        if (instance == null) {
            instance = new OrangeJson();
        }
        return instance;
    }

    private OrangeJson() {
        var gsonBuilder = new GsonBuilder();
        for (TypeAdapterFactory factory : ServiceLoader.load(TypeAdapterFactory.class)) {
            gsonBuilder.registerTypeAdapterFactory(factory);
        }
        gson = gsonBuilder
                .serializeNulls()
                .serializeSpecialFloatingPointValues()
                .create();
    }

    @NotNull
    public OrangeRecord deserialize(String orangeJson) {
        return gson.fromJson(orangeJson, OrangeRecord.class);
    }

    @NotNull
    public String serialize(OrangeRecord orangeRecord) {
        return gson.toJson(orangeRecord);
    }

    @NotNull
    public OrangeRecord read(@NotNull String orangeJsonFilePath) throws FileNotFoundException {
        var reader = new BufferedReader(new FileReader(orangeJsonFilePath));
        return gson.fromJson(reader, OrangeRecord.class);
    }

    public void write(@NotNull OrangeRecord orangeRecord, @NotNull String outputFilePath) throws IOException {
        BufferedWriter writer = new BufferedWriter(new FileWriter(outputFilePath));
        gson.toJson(orangeRecord, OrangeRecord.class, writer);
    }
}
