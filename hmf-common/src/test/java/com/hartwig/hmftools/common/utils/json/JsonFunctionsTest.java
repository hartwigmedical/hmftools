package com.hartwig.hmftools.common.utils.json;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;

import com.google.gson.JsonArray;
import com.google.gson.JsonObject;

import org.junit.Test;

public class JsonFunctionsTest
{
    @Test
    public void allJsonFunctionsWork()
    {
        JsonObject object = new JsonObject();

        assertNull(JsonFunctions.optionalJsonObject(object, "any"));
        assertNull(JsonFunctions.optionalJsonArray(object, "any"));
        assertNull(JsonFunctions.optionalNullableString(object, "any"));
        assertNull(JsonFunctions.optionalString(object, "any"));
        assertNull(JsonFunctions.optionalStringList(object, "any"));

        object.add("emptyArray", new JsonArray());
        assertNotNull(JsonFunctions.optionalJsonArray(object, "emptyArray"));

        JsonArray array = new JsonArray();
        array.add("value1");
        array.add("value2");
        object.add("array", array);
        assertEquals(2, JsonFunctions.stringList(object, "array").size());

        object.addProperty("string", "value");
        assertEquals("value", JsonFunctions.string(object, "string"));

        object.addProperty("nullableString", (String) null);
        assertNull(JsonFunctions.nullableString(object, "nullableString"));

        object.add("object", new JsonObject());
        assertNotNull(JsonFunctions.optionalJsonObject(object, "object"));
    }
}