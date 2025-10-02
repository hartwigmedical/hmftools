package com.hartwig.hmftools.cobalt.count;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotEquals;
import static org.junit.Assert.assertTrue;

import org.junit.Test;

public class DepthReadingTest
{
    @Test
    public void testConstructorAndGetters()
    {
        String chromosome = "X";
        int startPosition = 1000;
        double readDepth = 42.5;
        double readGcContent = 0.75;
        
        DepthReading reading = new DepthReading(chromosome, startPosition, readDepth, readGcContent);
        
        assertEquals(chromosome, reading.Chromosome);
        assertEquals(startPosition, reading.StartPosition);
        assertEquals(readDepth, reading.ReadDepth, 0.0001);
        assertEquals(readGcContent, reading.ReadGcContent, 0.0001);
    }
    
    @Test
    public void testEquals()
    {
        DepthReading reading1 = new DepthReading("1", 1000, 50.0, 0.5);
        DepthReading reading2 = new DepthReading("1", 1000, 50.0, 0.5);
        DepthReading reading3 = new DepthReading("2", 1000, 50.0, 0.5);
        DepthReading reading4 = new DepthReading("1", 2000, 50.0, 0.5);
        DepthReading reading5 = new DepthReading("1", 1000, 60.0, 0.5);
        DepthReading reading6 = new DepthReading("1", 1000, 50.0, 0.6);
        
        assertEquals(reading1, reading1);
        
        assertEquals(reading1, reading2);
        assertEquals(reading2, reading1);
        
        // Test with different values
        assertNotEquals(reading1, reading3);
        assertNotEquals(reading1, reading4);
        assertNotEquals(reading1, reading5);
        assertNotEquals(reading1, reading6);
        
        // Test with null and different type
        assertNotEquals(null, reading1);
        assertNotEquals("Not a DepthReading", reading1);
    }
    
    @Test
    public void testHashCode()
    {
        DepthReading reading1 = new DepthReading("1", 1000, 50.0, 0.5);
        DepthReading reading2 = new DepthReading("1", 1000, 50.0, 0.5);

        assertEquals(reading1.hashCode(), reading2.hashCode());
    }
    
    @Test
    public void testToString()
    {
        DepthReading reading = new DepthReading("1", 1000, 50.0, 0.5);
        String toString = reading.toString();
        
        assertTrue(toString.contains("Chromosome='1'"));
        assertTrue(toString.contains("StartPosition=1000"));
        assertTrue(toString.contains("ReadDepth=50.0"));
        assertTrue(toString.contains("ReadGcContent=0.5"));
    }
}