package com.hartwig.hmftools.finding.datamodel;

import static org.junit.Assert.assertEquals;

import java.io.IOException;
import java.nio.file.Path;

import org.junit.Test;

public class FindingsJsonTest
{
    @Test
    public void testWriteRead() throws IOException
    {
        Path filePath = Path.of("findings.json");
        FindingRecord findingRecord = TestFindingRecordFactory.createMinimalTestFindingRecord();
        FindingsJson findingsJson = new FindingsJson();
        findingsJson.write(findingRecord, filePath);
        FindingRecord readFindingRecord = findingsJson.read(filePath);
        assertEquals(findingRecord, readFindingRecord);
    }
}
