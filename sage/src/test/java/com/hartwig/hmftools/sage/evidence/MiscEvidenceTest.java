package com.hartwig.hmftools.sage.evidence;

import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.MockRefGenome.generateRandomBases;
import static com.hartwig.hmftools.common.test.SamRecordTestUtils.buildDefaultBaseQuals;
import static com.hartwig.hmftools.sage.SageConstants.DEFAULT_FLANK_LENGTH;
import static com.hartwig.hmftools.sage.SageConstants.DEFAULT_MAX_READ_DEPTH;
import static com.hartwig.hmftools.sage.common.TestUtils.QUALITY_CALCULATOR;
import static com.hartwig.hmftools.sage.common.TestUtils.READ_ID_GENERATOR;
import static com.hartwig.hmftools.sage.common.TestUtils.REF_BASES_200;
import static com.hartwig.hmftools.sage.common.TestUtils.REF_SEQUENCE_200;
import static com.hartwig.hmftools.sage.common.TestUtils.TEST_CONFIG;
import static com.hartwig.hmftools.sage.common.TestUtils.TEST_SAMPLE;
import static com.hartwig.hmftools.sage.common.TestUtils.buildCigarString;
import static com.hartwig.hmftools.sage.common.TestUtils.buildSamRecord;
import static com.hartwig.hmftools.sage.common.TestUtils.createSamRecord;
import static com.hartwig.hmftools.sage.common.VariantTier.LOW_CONFIDENCE;
import static com.hartwig.hmftools.sage.common.VariantUtils.createReadContext;
import static com.hartwig.hmftools.sage.common.VariantUtils.createReadCounter;
import static com.hartwig.hmftools.sage.common.VariantUtils.createSimpleVariant;
import static com.hartwig.hmftools.sage.evidence.SplitReadSegment.formSegment;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.sage.common.RegionTaskTester;
import com.hartwig.hmftools.sage.common.SageVariant;
import com.hartwig.hmftools.sage.common.SimpleVariant;
import com.hartwig.hmftools.sage.common.VariantReadContext;
import com.hartwig.hmftools.sage.common.VariantReadContextBuilder;
import com.hartwig.hmftools.sage.pipeline.RegionTask;

import org.junit.Test;

import htsjdk.samtools.SAMRecord;
import junit.framework.TestCase;

public class MiscEvidenceTest
{
    @Test
    public void testLongInsertPartialReadMatches()
    {
        int position = 35;

        String insertedBases = "TTTTTTTTTTT";

        String refBase = REF_BASES_200.substring(35, 36);
        SimpleVariant variant = createSimpleVariant(position, refBase, refBase + insertedBases);

        String readBases = REF_BASES_200.substring(1, 35) + variant.Alt + REF_BASES_200.substring(37, 60);
        String readCigar = "35M11I24M";

        SAMRecord initialRead = createSamRecord(READ_ID_GENERATOR.nextId(), CHR_1, 1, readBases, readCigar);

        VariantReadContextBuilder builder = new VariantReadContextBuilder(DEFAULT_FLANK_LENGTH);

        VariantReadContext readContext = builder.createContext(variant, initialRead, 34, REF_SEQUENCE_200);

        assertEquals(11, readContext.VarIndex);
        assertTrue(readContext.isValid());
        assertEquals(readBases.substring(33, 48), readContext.coreStr());

        ReadContextCounter readCounter = new ReadContextCounter(
                0, readContext, LOW_CONFIDENCE,
                DEFAULT_MAX_READ_DEPTH, 1, TEST_CONFIG, QUALITY_CALCULATOR, TEST_SAMPLE, true);

        readCounter.processRead(initialRead, 1, null);

        assertEquals(1, readCounter.readCounts().Full);

        // now a read which matches the read context by 2 or more bases than then ref without being a full core match
        String altReadBases = REF_BASES_200.substring(1, 35) + variant.Alt.substring(0, 8) + "GGGG" + REF_BASES_200.substring(37, 60);
        readCigar = buildCigarString(altReadBases.length());

        SAMRecord altRead = createSamRecord(READ_ID_GENERATOR.nextId(), CHR_1, 1, altReadBases, readCigar);

        readCounter.processRead(altRead, 1, null);

        assertEquals(1, readCounter.readCounts().PartialCore);
    }

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

    @Test
    public void testLowQualSoftClip()
    {
        int position = 100;

        VariantReadContext readContext = createReadContext(position, "A", "T");

        ReadContextCounter readContextCounter = createReadCounter(0, readContext);

        String altReadBases = readContext.readBases() + REF_BASES_200.substring(0, 20);
        String readCigar = "15S20M";
        int readVarIndex = readContext.VarIndex;
        int readPosStart = position - readVarIndex + 15; // for the soft-clip

        SAMRecord altRead = createSamRecord(READ_ID_GENERATOR.nextId(), CHR_1, readPosStart, altReadBases, readCigar);

        readContextCounter.processRead(altRead, 1, null);

        assertEquals(1, readContextCounter.readCounts().Full);

        // again but with too many low-qual soft-clip bases
        for(int i = 0; i < 10; ++i)
        {
            altRead.getBaseQualities()[i] = 11;
        }

        ReadMatchType readMatchType = readContextCounter.processRead(altRead, 1, null);
        assertEquals(readMatchType, ReadMatchType.SOFT_CLIP);
        assertEquals(1, readContextCounter.readCounts().Full);
    }

    @Test
    public void testMnvBaseQuality()
    {
        int position = 100;

        VariantReadContext readContext = createReadContext(position, "ATG", "TGA");

        ReadContextCounter readContextCounter = createReadCounter(0, readContext);

        String altReadBases = REF_BASES_200.substring(0, 10) + readContext.readBases() + REF_BASES_200.substring(0, 10);
        String readCigar = buildCigarString(altReadBases.length());
        int readVarIndex = 10 + readContext.VarIndex;
        int readPosStart = position - readVarIndex;

        SAMRecord altRead = createSamRecord(READ_ID_GENERATOR.nextId(), CHR_1, readPosStart, altReadBases, readCigar);

        readContextCounter.processRead(altRead, 1, null);

        assertEquals(37, readContextCounter.qualCounters().altBaseQualityTotal());

        // min rather than average is used
        altRead.getBaseQualities()[readVarIndex] = 11;

        readContextCounter.processRead(altRead, 1, null);

        assertEquals(48, readContextCounter.qualCounters().altBaseQualityTotal());
    }

    @Test
    public void testReadContextCounterOrdering()
    {
        ChrBaseRegion region = new ChrBaseRegion(CHR_1, 1, 300);

        RegionTaskTester tester = new RegionTaskTester();
        String refBases = REF_BASES_200 + generateRandomBases(1500);
        tester.RefGenome.RefGenomeMap.put(CHR_1, refBases); // since expects region to be 1300+

        RegionTask task = tester.createRegionTask(region);

        // test a collection of variants, not phased

        // reads which establish the variants
        String readBases1 = REF_BASES_200.substring(30, 50) + "A" + REF_BASES_200.substring(51, 70);
        String readCigar = "40M";
        SAMRecord read1 = buildSamRecord(30, readCigar, readBases1, buildDefaultBaseQuals(readBases1.length()));

        String readBases2 = REF_BASES_200.substring(30, 48) + REF_BASES_200.substring(52, 70); // 4-base delete
        readCigar = "18M4D17M";
        SAMRecord read2 = buildSamRecord(30, readCigar, readBases2, buildDefaultBaseQuals(readBases2.length()));

        String readBases3 = REF_BASES_200.substring(40, 55) + "G" + REF_BASES_200.substring(56, 80);
        readCigar = "40M";
        SAMRecord read3 = buildSamRecord(40, readCigar, readBases3, buildDefaultBaseQuals(readBases3.length()));

        tester.TumorSamSlicer.ReadRecords.add(read1);
        tester.TumorSamSlicer.ReadRecords.add(read1);

        tester.TumorSamSlicer.ReadRecords.add(read2);
        tester.TumorSamSlicer.ReadRecords.add(read2);

        tester.TumorSamSlicer.ReadRecords.add(read3);
        tester.TumorSamSlicer.ReadRecords.add(read3);

        // a read beyond the delete but still considered
        String readBases4 = REF_BASES_200.substring(47, 48) + REF_BASES_200.substring(52, 80);
        readCigar = "29M"; // could have been soft-clipped or an SNV
        SAMRecord read4 = buildSamRecord(51, readCigar, readBases4, buildDefaultBaseQuals(readBases4.length()));

        tester.TumorSamSlicer.ReadRecords.add(read4);

        // Configurator.setRootLevel(Level.TRACE);

        task.run();

        TestCase.assertEquals(3, task.getVariants().size());
        SageVariant del = task.getVariants().stream().filter(x -> x.position() == 47).findFirst().orElse(null);
        SageVariant var2 = task.getVariants().stream().filter(x -> x.position() == 50).findFirst().orElse(null);
        SageVariant var3 = task.getVariants().stream().filter(x -> x.position() == 55).findFirst().orElse(null);
        TestCase.assertNotNull(del);
        TestCase.assertNotNull(var2);
        TestCase.assertNotNull(var3);

        ReadContextCounter delRcCounter = del.tumorReadCounters().get(0);

        TestCase.assertEquals(2, delRcCounter.readSupportCounts().Full);
        TestCase.assertEquals(7, delRcCounter.readSupportCounts().Total);
        TestCase.assertEquals(7, delRcCounter.depth());
    }
}
