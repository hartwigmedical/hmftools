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
        jsonWriter.beginObject();
        jsonWriter.name("year").value(localDate.getYear());
        jsonWriter.name("month").value(localDate.getMonthValue());
        jsonWriter.name("day").value(localDate.getDayOfMonth());
        jsonWriter.endObject();
    }

    @Override
    public LocalDate read(@NotNull JsonReader jsonReader) throws IOException
    {
        if(jsonReader.peek() == JsonToken.NULL)
        {
            jsonReader.nextNull();
            return null;
        }

        int year = -1;
        int month = -1;
        int day = -1;

        jsonReader.beginObject();
        while(jsonReader.hasNext())
        {
            String name = jsonReader.nextName();
            if(name.equals("year"))
            {
                year = jsonReader.nextInt();
            }
            else if(name.equals("month"))
            {
                month = jsonReader.nextInt();
            }
            else if(name.equals("day"))
            {
                day = jsonReader.nextInt();
            }
            else
            {
                jsonReader.skipValue(); // Ignore unexpected fields
            }
        }
        jsonReader.endObject();

        if(year == -1 && month == -1 && day == -1)
        {
            throw new IllegalArgumentException("Invalid JSON format");
        }

        return LocalDate.of(year, month, day);
    }
}