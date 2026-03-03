package com.hartwig.hmftools.amber;

import org.junit.Assert;
import org.junit.Test;

public class PositionEvidenceTest
{
    @Test
    public void equalsTest()
    {
        PositionEvidence pe0 = new PositionEvidence("1", 100, "A", "C");
        PositionEvidence pe1 = new PositionEvidence("1", 100, "A", "C");
        PositionEvidence pe2 = new PositionEvidence("1", 100, "A", "G");
        PositionEvidence pe3 = new PositionEvidence("1", 100, "T", "C");
        PositionEvidence pe4 = new PositionEvidence("1", 101, "A", "C");
        PositionEvidence pe5 = new PositionEvidence("2", 100, "A", "C");
        Assert.assertEquals(pe0, pe1);
        Assert.assertNotEquals(pe0, pe2);
        Assert.assertNotEquals(pe0, pe3);
        Assert.assertNotEquals(pe0, pe4);
        Assert.assertNotEquals(pe0, pe5);

        PositionEvidence pe6 = new PositionEvidence("1", 100, "A", "C");
        Assert.assertEquals(pe0, pe6);
        pe6.ReadDepth = 10;
        Assert.assertEquals(pe0, pe6);
        pe6.ReadDepth = 1;
        Assert.assertEquals(pe0, pe6);
        pe6.RefSupport = 90;
        Assert.assertEquals(pe0, pe6);
        pe6.AltSupport = 10;
        Assert.assertEquals(pe0, pe6);
    }

    @Test
    public void hashCodeTest()
    {
        PositionEvidence pe0 = new PositionEvidence("1", 100, "A", "C");
        PositionEvidence pe1 = new PositionEvidence("1", 100, "A", "C");
        Assert.assertEquals(pe0.hashCode(), pe1.hashCode());
        pe0.AltSupport = 10;
        pe0.RefSupport = 90;
        pe0.ReadDepth = 100;
        Assert.assertEquals(pe0.hashCode(), pe1.hashCode());
    }

    @Test
    public void convertToVafReadingTest()
    {
        PositionEvidence pe0 = new PositionEvidence("1", 1000, "A", "C");
        pe0.AltSupport = 10;
        pe0.RefSupport = 90;
        pe0.ReadDepth = 100;
        VafReading converted = pe0.convertToVafReading();
        Assert.assertEquals("1", converted.chromosome());
        Assert.assertEquals(1000, converted.position());
        Assert.assertEquals(100, converted.readDepth());
        Assert.assertEquals(90, converted.refSupport());
        Assert.assertEquals(10, converted.altSupport());
    }

    @Test
    public void vafTest()
    {
        PositionEvidence pe0 = new PositionEvidence("1", 1000, "A", "C");
        Assert.assertEquals(Double.NaN, pe0.vaf(), 0.0001);
        pe0.AltSupport = 10;
        pe0.RefSupport = 90;
        pe0.ReadDepth = 100;
        Assert.assertEquals(0.1, pe0.vaf(), 0.0001);
    }

    @Test
    public void symmetricVafTest()
    {
        PositionEvidence pe0 = new PositionEvidence("1", 1000, "A", "C");
        Assert.assertEquals(Double.NaN, pe0.symmetricVaf(), 0.0001);
        pe0.AltSupport = 90;
        pe0.RefSupport = 10;
        pe0.ReadDepth = 100;
        Assert.assertEquals(0.1, pe0.symmetricVaf(), 0.0001);
    }
}
