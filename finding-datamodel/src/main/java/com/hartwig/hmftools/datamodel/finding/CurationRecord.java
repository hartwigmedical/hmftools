package com.hartwig.hmftools.datamodel.finding;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.nio.file.Files;
import java.nio.file.Path;
import java.time.LocalDate;
import java.util.List;
import java.util.ServiceLoader;

import com.google.gson.GsonBuilder;
import com.google.gson.TypeAdapterFactory;
import com.hartwig.hmftools.datamodel.LocalDateAdapter;

import org.jetbrains.annotations.NotNull;

import io.soabase.recordbuilder.core.RecordBuilder;

@RecordBuilder
public record CurationRecord(
        @NotNull List<DriverCuration> driverCurations
)
{

    @NotNull
    public static CurationRecord read(@NotNull Path findingsJsonFilePath) throws IOException
    {
        try(BufferedReader bufferedReader = Files.newBufferedReader(findingsJsonFilePath))
        {
            return gson().fromJson(bufferedReader, CurationRecord.class);
        }
    }

    @NotNull
    public static CurationRecord read(@NotNull InputStream stream) throws IOException
    {
        try(BufferedReader bufferedReader = new BufferedReader(new InputStreamReader(stream)))
        {
            return gson().fromJson(bufferedReader, CurationRecord.class);
        }
    }

    public void write(@NotNull CurationRecord curationRecord, @NotNull Path outputFilePath) throws IOException
    {
        try(BufferedWriter writer = Files.newBufferedWriter(outputFilePath))
        {
            gson().toJson(curationRecord, CurationRecord.class, writer);
        }
    }

    private static com.google.gson.Gson gson()
    {
        GsonBuilder gsonBuilder = new GsonBuilder();
        for(TypeAdapterFactory factory : ServiceLoader.load(TypeAdapterFactory.class))
        {
            gsonBuilder.registerTypeAdapterFactory(factory);
        }

        return gsonBuilder.serializeNulls().serializeSpecialFloatingPointValues().setPrettyPrinting()
                .registerTypeAdapter(LocalDate.class, new LocalDateAdapter())
                .create();
    }
}
