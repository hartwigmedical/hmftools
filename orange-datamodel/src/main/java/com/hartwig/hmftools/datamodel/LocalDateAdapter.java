package com.hartwig.hmftools.datamodel;

import java.io.IOException;
import java.time.LocalDate;

import com.google.gson.TypeAdapter;
import com.google.gson.stream.JsonReader;
import com.google.gson.stream.JsonToken;
import com.google.gson.stream.JsonWriter;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class LocalDateAdapter extends TypeAdapter<LocalDate>
{
    @Override
    public void write(@NotNull JsonWriter jsonWriter, @Nullable LocalDate localDate) throws IOException
    {
        if(localDate == null)
        {
            jsonWriter.nullValue();
            return;
        }
        jsonWriter.value(localDate.toString());
    }

    @Override
    public @Nullable LocalDate read(@NotNull JsonReader jsonReader) throws IOException
    {
        if(jsonReader.peek() == JsonToken.NULL)
        {
            jsonReader.nextNull();
            return null;
        }
        return LocalDate.parse(jsonReader.nextString());
    }
}
