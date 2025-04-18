package com.hartwig.hmftools.bamtools.remapper;

import static java.lang.Integer.parseInt;

import java.io.File;
import java.util.LinkedList;
import java.util.List;
import java.util.Objects;

import org.umccr.java.hellbender.utils.bwa.BwaMemAlignment;
import org.junit.Assert;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

public class RemapperTestBase
{

    static BwaMemAlignment bwa(String data, String xa, String remainingData)
    {
        String[] parts1 = data.split(",");
        String[] parts2 = remainingData.split(",");
        return new BwaMemAlignment(
                parseInt(parts1[0]),
                parseInt(parts1[1]),
                parseInt(parts1[2]),
                parseInt(parts1[3]),
                parseInt(parts1[4]),
                parseInt(parts1[5]),
                parseInt(parts1[6]),
                parseInt(parts1[7]),
                parseInt(parts1[8]),
                parseInt(parts1[9]),
                parts1[10],
                parts1[11],
                xa,
                parseInt(parts2[0]),
                parseInt(parts2[1]),
                parseInt(parts2[2])
        );
    }

    static BwaMemAlignment bwa(String data)
    {
        String[] parts = data.split(",");
        return new BwaMemAlignment(
                parseInt(parts[0]),
                parseInt(parts[1]),
                parseInt(parts[2]),
                parseInt(parts[3]),
                parseInt(parts[4]),
                parseInt(parts[5]),
                parseInt(parts[6]),
                parseInt(parts[7]),
                parseInt(parts[8]),
                parseInt(parts[9]),
                parts[10],
                parts[11],
                parts[12],
                parseInt(parts[13]),
                parseInt(parts[14]),
                parseInt(parts[15])
        );
    }

    List<SAMRecord> records = readTestFile();

    SAMFileHeader samFileHeader()
    {
        return readTestFile().get(0).getHeader();
    }

    List<SAMRecord> readTestFile()
    {
        return readSamFile(getTestFile("tiny.sam"));
    }

    List<SAMRecord> readSamFile(File samFile)
    {
        List<SAMRecord> records = new LinkedList<>();
        try(SamReader samReader = SamReaderFactory.makeDefault().open(samFile))
        {
            samReader.forEach(records::add);
        }
        catch(Exception e)
        {
            throw new RuntimeException(e);
        }
        return records;
    }

    List<SAMSequenceRecord> readDictionarySequences(File samFile)
    {
        try(SamReader samReader = SamReaderFactory.makeDefault().open(samFile))
        {
            return new LinkedList<>(samReader.getFileHeader().getSequenceDictionary().getSequences());
        }
        catch(Exception e)
        {
            throw new RuntimeException(e);
        }
    }

    File getTestFile(String name)
    {
        ClassLoader classLoader = getClass().getClassLoader();
        return new File(Objects.requireNonNull(classLoader.getResource(name)).getFile());
    }

    void check(SAMRecord expected, SAMRecord actual)
    {
        // We do the same checks as would be done by using SAMRecord.equals()
        // except that we do not check attributes.
        Assert.assertEquals(expected.getAlignmentStart(), actual.getAlignmentStart());
        Assert.assertEquals(expected.getFlags(), actual.getFlags());
        Assert.assertEquals(expected.getInferredInsertSize(), actual.getInferredInsertSize());
        Assert.assertEquals(expected.getMappingQuality(), actual.getMappingQuality());
        Assert.assertEquals(expected.getMateAlignmentStart(), actual.getMateAlignmentStart());
        Assert.assertEquals(expected.getMateReferenceIndex(), actual.getMateReferenceIndex());
        Assert.assertEquals(expected.getReferenceIndex(), actual.getReferenceIndex());
        Assert.assertEquals(expected.getReadName(), actual.getReadName());
        Assert.assertArrayEquals(expected.getBaseQualities(), actual.getBaseQualities());
        Assert.assertEquals(expected.getCigar(), actual.getCigar());
        Assert.assertEquals(expected.getMateReferenceName(), actual.getMateReferenceName());
        Assert.assertArrayEquals(expected.getReadBases(), actual.getReadBases());
        Assert.assertEquals(expected.getReferenceName(), actual.getReferenceName());
    }
}

