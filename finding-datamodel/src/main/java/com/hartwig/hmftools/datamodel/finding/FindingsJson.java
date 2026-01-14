package com.hartwig.hmftools.datamodel.finding;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.Reader;
import java.nio.file.Files;
import java.nio.file.Path;
import java.time.LocalDate;
import java.util.ServiceLoader;

import com.google.gson.Gson;
import com.google.gson.GsonBuilder;
import com.google.gson.TypeAdapterFactory;
import com.hartwig.hmftools.datamodel.LocalDateAdapter;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class FindingsJson {

    @Nullable
    private static FindingsJson instance;

    private final Gson gson;

    @NotNull
    public static FindingsJson getInstance() {
        if (instance == null) {
            instance = new FindingsJson();
        }
        return instance;
    }

    private FindingsJson() {
        GsonBuilder gsonBuilder = new GsonBuilder();
        for (TypeAdapterFactory factory : ServiceLoader.load(TypeAdapterFactory.class)) {
            gsonBuilder.registerTypeAdapterFactory(factory);
        }

        gson = gsonBuilder.serializeNulls().serializeSpecialFloatingPointValues().setPrettyPrinting()
                .registerTypeAdapter(LocalDate.class, new LocalDateAdapter())
                .create();
    }

    @NotNull
    public FindingRecord read(@NotNull Path findingsJsonFilePath) throws IOException {
        return read(Files.newBufferedReader(findingsJsonFilePath));
    }

    @NotNull
    public FindingRecord read(@NotNull InputStream inputStream) throws IOException {
        try (BufferedReader bufferedReader = new BufferedReader(new InputStreamReader(inputStream))) {
            return read(bufferedReader);
        }
    }

    @NotNull
    public FindingRecord read(@NotNull Reader reader) throws IOException {
        return gson.fromJson(reader, FindingRecord.class);
    }

    public void write(@NotNull FindingRecord findingRecord, @NotNull Path outputFilePath) throws IOException {
        try (BufferedWriter writer = Files.newBufferedWriter(outputFilePath)) {
            gson.toJson(findingRecord, FindingRecord.class, writer);
        }
    }
}
