package com.hartwig.hmftools.sage.evidence;

import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.sage.common.TestUtils.READ_ID_GENERATOR;
import static com.hartwig.hmftools.sage.common.TestUtils.REF_BASES_200;
import static com.hartwig.hmftools.sage.common.TestUtils.buildCigarString;
import static com.hartwig.hmftools.sage.common.TestUtils.createSamRecord;
import static com.hartwig.hmftools.sage.common.VariantUtils.createReadContext;
import static com.hartwig.hmftools.sage.common.VariantUtils.createReadCounter;

import static org.junit.Assert.assertEquals;

import com.hartwig.hmftools.sage.common.VariantReadContext;

import org.junit.Test;

import htsjdk.samtools.SAMRecord;

public class MiscEvidenceTest
{
    @Test
    public void testReadEdgeDistancePenalty()
    {
        int position = 100;

        VariantReadContext readContext = createReadContext(position, "A", "T");

        ReadContextCounter readContextCounter = createReadCounter(0, readContext);

        String altReadBases = REF_BASES_200.substring(0, 10) + readContext.readBases() + REF_BASES_200.substring(0, 10);
        String readCigar = buildCigarString(altReadBases.length());
        int readVarIndex = 10 + readContext.VarIndex;
        int readPosStart = position - readVarIndex;

        SAMRecord altRead = createSamRecord(READ_ID_GENERATOR.nextId(), CHR_1, readPosStart, altReadBases, readCigar);

        readContextCounter.processRead(altRead, 1, null);

        assertEquals(1, readContextCounter.readCounts().Full);
        assertEquals(25, readContextCounter.readQuals().Full);

        // now a read right up against the position from the left
        altReadBases = readContext.readBases().substring(readContext.VarIndex) + REF_BASES_200.substring(0, 30);
        readCigar = buildCigarString(altReadBases.length());
        readPosStart = position;

        altRead = createSamRecord(READ_ID_GENERATOR.nextId(), CHR_1, readPosStart, altReadBases, readCigar);

        readContextCounter.processRead(altRead, 1, null);

        assertEquals(1, readContextCounter.readCounts().PartialCore);
        assertEquals(10, readContextCounter.readQuals().PartialCore);

        // now 1 base in from the edge
        altReadBases = readContext.readBases().substring(readContext.VarIndex - 1) + REF_BASES_200.substring(0, 30);
        readCigar = buildCigarString(altReadBases.length());
        readPosStart = position - 1;

        altRead = createSamRecord(READ_ID_GENERATOR.nextId(), CHR_1, readPosStart, altReadBases, readCigar);

        readContextCounter.processRead(altRead, 1, null);

        assertEquals(2, readContextCounter.readCounts().PartialCore);
        assertEquals(30, readContextCounter.readQuals().PartialCore);

        // and from the other side
        altReadBases = REF_BASES_200.substring(0, 30) + readContext.readBases().substring(0, readContext.VarIndex + 1);
        readCigar = buildCigarString(altReadBases.length());
        readPosStart = position - 30 - readContext.VarIndex;

        altRead = createSamRecord(READ_ID_GENERATOR.nextId(), CHR_1, readPosStart, altReadBases, readCigar);

        readContextCounter.processRead(altRead, 1, null);

        assertEquals(3, readContextCounter.readCounts().PartialCore);
        assertEquals(40, readContextCounter.readQuals().PartialCore);
    }
}
