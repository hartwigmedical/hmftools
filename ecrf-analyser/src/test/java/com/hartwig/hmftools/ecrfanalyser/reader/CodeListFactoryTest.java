package com.hartwig.hmftools.ecrfanalyser.reader;

import static org.junit.Assert.assertEquals;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;

import org.junit.Test;

public class CodeListFactoryTest {

    @Test
    public void canExtractValuesFromStrings() {
        List<String> codeListItems = Lists.newArrayList("1=x", "2= y");
        Map<Integer, String> values = CodeListFactory.extractValuesFromStrings(codeListItems);

        assertEquals(2, values.size());
        assertEquals("x", values.get(1));
        assertEquals("y", values.get(2));
    }
}