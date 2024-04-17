package com.hartwig.hmftools.sage.old;

import static com.hartwig.hmftools.sage.SageConstants.DEFAULT_READ_CONTEXT_FLANK_SIZE;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;

import com.hartwig.hmftools.sage.common.RefSequence;
import com.hartwig.hmftools.sage.old.ReadContext;
import com.hartwig.hmftools.sage.old.ReadContextFactory;

import org.jetbrains.annotations.NotNull;
import org.junit.Before;
import org.junit.Test;

import htsjdk.samtools.SAMRecord;

public class ReadContextFactoryTest
{
    private ReadContextFactory victim;

    @Before
    public void setup()
    {
        victim = new ReadContextFactory(25);
    }

    @Test
    public void testSimpleSnvHas5BaseCore()
    {
        String refSequence = "GATCATCTAGG";
        String readSequence = "GATCACCTAGG";
        RefSequence refBases = new RefSequence(1000, refSequence.getBytes());

        SAMRecord record = buildSamRecord("11M", readSequence);
        ReadContext victim = this.victim.createSNVContext(1005, 5, record, refBases);
        assertEquals("CACCT", victim.coreString());
    }

    @Test
    public void testSimpleInsert()
    {
        String refSequence = "GATCATCTAGG";
        String readSequence = "GAGGCTCATCTAGG";
        RefSequence refBases = new RefSequence(1000, refSequence.getBytes());

        SAMRecord record = buildSamRecord("2M3I9M", readSequence);
        ReadContext victim = this.victim.createInsertContext("AGGC", 1000, 1, record.getReadBases(), refBases);
        assertEquals("GAGGCTC", victim.coreString());
    }

    @Test
    public void testInsertInRepeat()
    {
        String refSequence = "TGAAAAAAAATCT";
        String readSequence = "TGAAAAAAAAATCT";
        RefSequence refBases = new RefSequence(1000, refSequence.getBytes());

        SAMRecord record = buildSamRecord("2M1I11M", readSequence);
        ReadContext victim = this.victim.createInsertContext("GA", 1000, 1, record.getReadBases(), refBases);
        assertEquals("TGAAAAAAAAAT", victim.coreString());
    }

    @Test
    public void testInsertAtHomology()
    {
        String refSequence = "GATCATCTG";
        String readSequence = "GATCATCATCTG";
        RefSequence refBases = new RefSequence(1000, refSequence.getBytes());

        SAMRecord record = buildSamRecord("1M3I8M", readSequence);
        ReadContext victim = this.victim.createInsertContext("ATCA", 1000, 1, record.getReadBases(), refBases);
        assertEquals("GATCATCATCT", victim.coreString());
    }

    @Test
    public void testInsertAtHomologyRepeat()
    {
        String refSequence = "GATCATCATCTG";
        String readSequence = "GATCATCATCATCTG";
        RefSequence refBases = new RefSequence(1000, refSequence.getBytes());

        SAMRecord record = buildSamRecord("1M3I10M", readSequence);
        ReadContext victim = this.victim.createInsertContext("ATCA", 1000, 1, record.getReadBases(), refBases);
        assertEquals("GATCATCATCATCT", victim.coreString());
    }

    @Test
    public void testInsertAtHomologyWithAdditionalBases()
    {
        String refSequence = "ATGCGATCTTCC";
        String readSequence = "ATGCGATCAATCTTCC";
        RefSequence refBases = new RefSequence(1000, refSequence.getBytes());

        SAMRecord record = buildSamRecord("5M4I7M", readSequence);
        ReadContext victim = this.victim.createInsertContext("GATCA", 1000, 4, record.getReadBases(), refBases);
        assertEquals("GCGATCAAT", victim.coreString());
    }

    private static final ReadContextFactory READ_CONTEXT_FACTORY = new ReadContextFactory(DEFAULT_READ_CONTEXT_FLANK_SIZE);

    @Test
    public void testSimpleDelete1()
    {
        // variant: pos 1025 AT>A
        String refBases =  "ATCTCTCAATGTTGACGGACAGCCTATTTTTGCCAATATCACACTGCCAGGT";
        String readBases = "ATCTCTCAATGTTGACGGACAGCCTATTTTGCCAATATCACACTGCCAGGT";
        RefSequence refSequence = new RefSequence(1000, refBases.getBytes());

        ReadContext readContext = READ_CONTEXT_FACTORY.createDelContext("AT", 1025, 25, readBases.getBytes(), refSequence);

        assertFalse(readContext.hasIncompleteCore());
        assertEquals("CTATTTTGC", readContext.coreString());
        assertEquals("GACGGACAGC", readContext.leftFlankString());
        assertEquals("CAATATCACA", readContext.rightFlankString());
    }

    @Test
    public void testDeleteAtHomology()
    {
        String refSequence = "GATCGGATCGCTT";
        String readSequence = "GATCGCTT";
        RefSequence refBases = new RefSequence(1000, refSequence.getBytes());

        SAMRecord record = buildSamRecord("1M5D7M", readSequence);
        ReadContext victim = this.victim.createDelContext("GATCGG", 1000, 0, record.getReadBases(), refBases);
        assertEquals("GATCGC", victim.coreString());
    }

    @Test
    public void testDeleteAtHomologyRepeat()
    {
        String refSequence = "GATCACCATCTG";
        String readSequence = "GATCATCTG";
        RefSequence refBases = new RefSequence(1000, refSequence.getBytes());

        SAMRecord record = buildSamRecord("1M3D8M", readSequence);
        ReadContext victim = this.victim.createDelContext("ATCA", 1000, 1, record.getReadBases(), refBases);
        assertEquals("GATCATCT", victim.coreString());
    }

    @Test
    public void testDeleteAtRepeatInRef()
    {
        String refSequence = "GATCATCATCTG";
        String readSequence = "GATCATCTG";
        RefSequence refBases = new RefSequence(1000, refSequence.getBytes());

        SAMRecord record = buildSamRecord("1M3D8M", readSequence);
        ReadContext victim = this.victim.createDelContext("ATCA", 1000, 1, record.getReadBases(), refBases);
        assertEquals("GATCATCTG", victim.coreString());
    }

    @Test
    public void testDeleteOneBase()
    {
        String refSequence = "GATCATCTAGG";
        String readSequence = "GTCATCTAGG";
        RefSequence refBases = new RefSequence(1000, refSequence.getBytes());

        SAMRecord record = buildSamRecord("1M1D9M", readSequence);
        ReadContext victim = this.victim.createDelContext("GA", 1000, 0, record.getReadBases(), refBases);
        assertEquals("GTC", victim.coreString());
    }

    @Test
    public void testDeleteTwoBase()
    {
        String refSequence = "GATCATCTAGG";
        String readSequence = "GCATCTAGG";
        RefSequence refBases = new RefSequence(1000, refSequence.getBytes());

        SAMRecord record = buildSamRecord("1M2D8M", readSequence);
        ReadContext victim = this.victim.createDelContext("GAT", 1000, 0, record.getReadBases(), refBases);
        assertEquals("GCA", victim.coreString());
    }

    @NotNull
    static SAMRecord buildSamRecord(@NotNull final String cigar, @NotNull final String readString)
    {
        final StringBuilder qualityString = new StringBuilder();
        for(int i = 0; i < readString.length(); i++)
        {
            qualityString.append("A");
        }

        return buildSamRecord(cigar, readString, qualityString.toString());
    }

    @NotNull
    static SAMRecord buildSamRecord(@NotNull final String cigar, @NotNull final String readString, @NotNull final String qualities)
    {
        final SAMRecord record = new SAMRecord(null);
        record.setAlignmentStart(1000);
        record.setCigarString(cigar);
        record.setReadString(readString);
        record.setReadNegativeStrandFlag(false);
        record.setBaseQualityString(qualities);
        record.setMappingQuality(20);
        record.setDuplicateReadFlag(false);
        record.setReadUnmappedFlag(false);
        return record;
    }
}
