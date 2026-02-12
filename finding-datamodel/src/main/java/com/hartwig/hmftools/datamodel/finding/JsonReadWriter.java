package com.hartwig.hmftools.datamodel.finding;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.io.OutputStreamWriter;
import java.io.Reader;
import java.nio.file.Files;
import java.nio.file.Path;
import java.time.LocalDate;
import java.util.ServiceLoader;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

import com.google.gson.Gson;
import com.google.gson.GsonBuilder;
import com.google.gson.TypeAdapterFactory;
import com.hartwig.hmftools.datamodel.LocalDateAdapter;

import jakarta.validation.constraints.NotNull;

public class JsonReadWriter<T>
{
    private final Class<T> type;
    private final Gson gson;

    public JsonReadWriter(Class<T> type)
    {
        this.type = type;
        GsonBuilder gsonBuilder = new GsonBuilder();
        for(TypeAdapterFactory factory : ServiceLoader.load(TypeAdapterFactory.class))
        {
            gsonBuilder.registerTypeAdapterFactory(factory);
        }

        gson = gsonBuilder.serializeNulls().serializeSpecialFloatingPointValues().setPrettyPrinting()
                .registerTypeAdapter(LocalDate.class, new LocalDateAdapter())
                .create();
    }

    @NotNull
    public T read(@NotNull Path jsonFilePath) throws IOException
    {
        InputStream inputStream = Files.newInputStream(jsonFilePath);
        if(isGZIP(jsonFilePath))
        {
            inputStream = new GZIPInputStream(inputStream);
        }
        return read(inputStream);
    }

    @NotNull
    public T read(@NotNull InputStream inputStream) throws IOException
    {
        return read(new InputStreamReader(inputStream));
    }

    @NotNull
    public T read(@NotNull Reader reader) throws IOException
    {
        try(BufferedReader bufferedReader = new BufferedReader(reader))
        {
            return gson.fromJson(bufferedReader, type);
        }
    }

    public void write(@NotNull T findingRecord, @NotNull Path outputFilePath) throws IOException
    {
        OutputStream outputStream = Files.newOutputStream(outputFilePath);
        if(isGZIP(outputFilePath))
        {
            outputStream = new GZIPOutputStream(outputStream);
        }
        try(BufferedWriter writer = new BufferedWriter(new OutputStreamWriter(outputStream)))
        {
            gson.toJson(findingRecord, type, writer);
        }
    }

    private static boolean isGZIP(@NotNull Path path)
    {
        return path.toString().endsWith(".gz");
    }
}
