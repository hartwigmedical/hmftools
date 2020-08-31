package com.hartwig.hmftools.serve.hartwig;

import static org.junit.Assert.assertEquals;

import org.apache.logging.log4j.util.Strings;
import org.junit.Test;

public class HartwigProteinInterpreterTest {

    @Test
    public void canInterpretProteinAnnotation() {
        assertEquals(Strings.EMPTY, HartwigProteinInterpreter.interpretProtein(Strings.EMPTY));
        assertEquals(Strings.EMPTY, HartwigProteinInterpreter.interpretProtein(HartwigProteinInterpreter.IGNORE_PROTEIN_ANNOTATION));
        assertEquals("V600E", HartwigProteinInterpreter.interpretProtein("V600E"));
        assertEquals("V600E", HartwigProteinInterpreter.interpretProtein("p.V600E"));
        assertEquals("V600E", HartwigProteinInterpreter.interpretProtein("Val600Glu"));
        assertEquals("V600E", HartwigProteinInterpreter.interpretProtein("p.Val600Glu"));
    }
}