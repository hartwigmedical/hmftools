package com.hartwig.hmftools.datamodel;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.io.OutputStreamWriter;
import java.io.Reader;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.time.LocalDate;
import java.util.ServiceLoader;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

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
        GsonBuilder gsonBuilder = new GsonBuilder();
        for (TypeAdapterFactory factory : ServiceLoader.load(TypeAdapterFactory.class)) {
            gsonBuilder.registerTypeAdapterFactory(factory);
        }

        gson = gsonBuilder
                .serializeNulls()
                .serializeSpecialFloatingPointValues()
                .setPrettyPrinting()
                .registerTypeAdapter(LocalDate.class, new LocalDateAdapter())
                .create();
    }

    @NotNull
    public OrangeRecord read(@NotNull String orangeJsonFilePath) throws IOException {
        return read(Path.of(orangeJsonFilePath));
    }

    @NotNull
    public OrangeRecord read(@NotNull Path orangeJsonFilePath) throws IOException {
        // NOTE: must use Files.newInputStream to support other Path types such as gs://
        InputStream inputStream = Files.newInputStream(orangeJsonFilePath);
        if(orangeJsonFilePath.endsWith(".gz")) {
            inputStream = new GZIPInputStream(inputStream);
        }
        return read(new InputStreamReader(inputStream, StandardCharsets.UTF_8));
    }

    @NotNull
    public OrangeRecord read(@NotNull Reader reader) throws IOException {
        try (BufferedReader bufferedReader = new BufferedReader(reader)) {
            return gson.fromJson(bufferedReader, OrangeRecord.class);
        }
    }

    public void write(@NotNull OrangeRecord orangeRecord, @NotNull String outputFilePath) throws IOException {
        write(orangeRecord, Path.of(outputFilePath));
    }

    public void write(@NotNull OrangeRecord orangeRecord, @NotNull Path outputFilePath) throws IOException {
        boolean isGzipped = outputFilePath.endsWith(".gz");
        // NOTE: must use Files.newOutputStream to support other Path types such as gs://
        OutputStream outputStream = Files.newOutputStream(outputFilePath);
        if(isGzipped) {
            outputStream = new GZIPOutputStream(outputStream);
        }
        try (BufferedWriter writer = new BufferedWriter(new OutputStreamWriter(outputStream, StandardCharsets.UTF_8))) {
            gson.toJson(orangeRecord, OrangeRecord.class, writer);
        }
    }
}
