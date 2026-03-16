package com.hartwig.hmftools.amber;

import java.io.File;
import java.nio.file.Files;
import java.util.List;

import org.junit.Assert;
import org.junit.Test;

public class PositionEvidenceFileTest
{
    @Test
    public void writeReadTest() throws Exception
    {
        PositionEvidence pe1 = new PositionEvidence("1", 1234, "A", "G");
        pe1.AltSupport = 10;
        pe1.RefSupport = 90;
        pe1.ReadDepth = 100;
        pe1.IndelCount = 0;
        pe1.BaseQualFiltered = 1;
        pe1.MapQualFiltered = 3;

        PositionEvidence pe2 = new PositionEvidence("1", 2345, "C", "T");
        pe2.AltSupport = 90;
        pe2.RefSupport = 90;
        pe2.ReadDepth = 180;
        pe2.IndelCount = 1;
        pe2.BaseQualFiltered = 0;
        pe2.MapQualFiltered = 0;

        PositionEvidence pe3 = new PositionEvidence("X", 2345, "C", "A");
        pe3.AltSupport = 900;
        pe3.RefSupport = 90;
        pe3.ReadDepth = 990;
        pe3.IndelCount = 10;
        pe3.BaseQualFiltered = 4;
        pe3.MapQualFiltered = 5;

        List<PositionEvidence> data = List.of(pe1, pe2, pe3);

        File tempDir = Files.createTempDirectory("amber").toFile();
        tempDir.deleteOnExit();
        File destination = new File(tempDir, "pe.tsv.gz");
        PositionEvidenceFile.write(destination.getAbsolutePath(), data);

        List<PositionEvidence> read = PositionEvidenceFile.read(destination.getAbsolutePath());
        Assert.assertEquals(data, read);
    }
}
