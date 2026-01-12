package com.hartwig.hmftools.datamodel;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.io.OutputStreamWriter;
import java.io.Reader;
import java.nio.charset.StandardCharsets;
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

        gson = gsonBuilder.serializeNulls().serializeSpecialFloatingPointValues()
                .registerTypeAdapter(LocalDate.class, new LocalDateAdapter())
                .create();
    }

    @NotNull
    public OrangeRecord read(@NotNull String orangeJsonFilePath) throws IOException {
        InputStream inputStream = new FileInputStream(orangeJsonFilePath);
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
        boolean isGzipped = outputFilePath.endsWith(".gz");
        Gson gsonToUse = gson;
        OutputStream outputStream = new FileOutputStream(outputFilePath);
        if(isGzipped) {
            outputStream = new GZIPOutputStream(outputStream);

            // Pretty print if gzipped
            gsonToUse = gson.newBuilder().setPrettyPrinting().create();
        }
        try (BufferedWriter writer = new BufferedWriter(new OutputStreamWriter(outputStream, StandardCharsets.UTF_8))) {
            gsonToUse.toJson(orangeRecord, OrangeRecord.class, writer);
        }
    }
}
