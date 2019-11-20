package com.hartwig.hmftools.vicc.reader;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.Map;

import com.google.common.collect.Maps;
import com.google.gson.JsonObject;

import org.junit.Test;

public class ViccDatamodelCheckerTest {

    @Test
    public void checkingWorksAsExpected() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("A", true);
        map.put("B", false);

        ViccDatamodelChecker checker = new ViccDatamodelChecker("test", map);

        JsonObject object = new JsonObject();
        assertFalse(checker.check(object));

        object.addProperty("A", "test A");
        assertTrue(checker.check(object));

        object.addProperty("B", "test B");
        assertTrue(checker.check(object));

        object.addProperty("C", "test C");
        assertFalse(checker.check(object));
    }
}
