package com.hartwig.hmftools.common.chromosome;

import static org.junit.Assert.assertFalse;

import static junit.framework.TestCase.assertEquals;
import static junit.framework.TestCase.assertTrue;

import com.hartwig.hmftools.common.purple.gender.Gender;

import org.junit.Test;

public class HumanChromosomeTest {

    @Test
    public void testString() {
        assertEquals("1", HumanChromosome._1.toString());
        assertEquals("X", HumanChromosome._X.toString());
        assertEquals("Y", HumanChromosome._Y.toString());
    }

    @Test
    public void testSexChromosomes() {
        assertTrue(HumanChromosome._X.isAllosome());
        assertFalse(HumanChromosome._X.isAutosome());

        assertTrue(HumanChromosome._Y.isAllosome());
        assertFalse(HumanChromosome._Y.isAutosome());
    }

    @Test
    public void testIsDiploid() {
        assertTrue(HumanChromosome._X.isDiploid(Gender.FEMALE));
        assertTrue(HumanChromosome._X.isDiploid(Gender.MALE_KLINEFELTER));
        assertFalse(HumanChromosome._X.isDiploid(Gender.MALE));
        assertFalse(HumanChromosome._Y.isDiploid(Gender.MALE));
        assertFalse(HumanChromosome._Y.isDiploid(Gender.FEMALE));
    }
}
