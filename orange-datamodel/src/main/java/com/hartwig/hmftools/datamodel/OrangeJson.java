package com.hartwig.hmftools.datamodel;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.Reader;
import java.time.LocalDate;
import java.util.ServiceLoader;

import com.google.gson.Gson;
import com.google.gson.GsonBuilder;
import com.google.gson.TypeAdapterFactory;
import com.hartwig.hmftools.datamodel.orange.OrangeRecord;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

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
        GsonBuilder gsonBuilder = new GsonBuilder().setPrettyPrinting();
        for (TypeAdapterFactory factory : ServiceLoader.load(TypeAdapterFactory.class)) {
            gsonBuilder.registerTypeAdapterFactory(factory);
        }

        gson = gsonBuilder.serializeNulls().serializeSpecialFloatingPointValues()
                .registerTypeAdapter(LocalDate.class, new LocalDateAdapter())
                .create();
    }

    @NotNull
    public OrangeRecord read(@NotNull String orangeJsonFilePath) throws IOException {
        return read(new FileReader(orangeJsonFilePath));
    }

    @NotNull
    public OrangeRecord read(@NotNull Reader reader) throws IOException {
        try (BufferedReader bufferedReader = new BufferedReader(reader)) {
            return gson.fromJson(bufferedReader, OrangeRecord.class);
        }
    }

    public void write(@NotNull OrangeRecord orangeRecord, @NotNull String outputFilePath) throws IOException {
        try (BufferedWriter writer = new BufferedWriter(new FileWriter(outputFilePath))) {
            gson.toJson(orangeRecord, OrangeRecord.class, writer);
        }
    }
}
