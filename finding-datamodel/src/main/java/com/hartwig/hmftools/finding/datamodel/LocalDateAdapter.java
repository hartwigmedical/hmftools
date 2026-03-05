package com.hartwig.hmftools.finding.datamodel;

import java.io.IOException;
import java.time.LocalDate;

import com.google.gson.TypeAdapter;
import com.google.gson.stream.JsonReader;
import com.google.gson.stream.JsonToken;
import com.google.gson.stream.JsonWriter;

import org.jspecify.annotations.Nullable;

import jakarta.validation.constraints.NotNull;

// TODO: This is a copy of com.hartwig.hmftools.datamodel.LocalDateAdapter but did not want to add the dependency
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
        // localDate toString always write in ISO date format:
        // https://docs.oracle.com/javase/8/docs/api/java/time/LocalDate.html#toString--
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
        // localDate parse expects ISO date format:
        // https://docs.oracle.com/javase/8/docs/api/java/time/LocalDate.html#parse-java.lang.CharSequence-
        return LocalDate.parse(jsonReader.nextString());
    }
}
