package com.hartwig.hmftools.amber;

import static org.junit.Assert.assertEquals;

import com.hartwig.hmftools.common.amber.BaseDepth;
import com.hartwig.hmftools.common.amber.ModifiableBaseDepth;
import com.hartwig.hmftools.common.amber.ModifiableTumorBAF;
import com.hartwig.hmftools.common.utils.sam.SAMRecords;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

import htsjdk.samtools.SAMRecord;

public class TumorBAFFactoryTest
{

    @Test
    public void useQualityOfBaseAfterDel()
    {
        int minQuality = SAMRecords.getBaseQuality('J');

        final SAMRecord lowQualDel = buildSamRecord(1000, "1M1D1M", "CT", "FI");
        final SAMRecord highQualDel = buildSamRecord(1000, "1M1D1M", "CT", "FJ");
        final ModifiableTumorBAF victim = createDefault("5", 1001);

        new TumorBAFFactory(minQuality).addEvidence(victim, lowQualDel);
        assertEquals(0, victim.tumorReadDepth());

        new TumorBAFFactory(minQuality).addEvidence(victim, highQualDel);
        assertEquals(1, victim.tumorReadDepth());
    }

    @NotNull
    private static ModifiableTumorBAF createDefault(@NotNull final String chromosome, final long position)
    {
        final ModifiableBaseDepth normal = ModifiableBaseDepth.create()
                .setChromosome(chromosome)
                .setPosition(position)
                .setRef(BaseDepth.Base.A)
                .setRefSupport(3)
                .setAlt(BaseDepth.Base.T)
                .setAltSupport(3)
                .setReadDepth(6)
                .setIndelCount(0);
        return TumorBAFFactory.create(normal);
    }

    private SAMRecord buildSamRecord(
            final int alignmentStart, @NotNull final String cigar, @NotNull final String readString, @NotNull final String qualities)
    {
        final SAMRecord record = new SAMRecord(null);
        record.setAlignmentStart(alignmentStart);
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
