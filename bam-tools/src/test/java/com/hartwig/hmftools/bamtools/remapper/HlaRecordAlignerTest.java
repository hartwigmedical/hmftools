package com.hartwig.hmftools.bamtools.remapper;

import static com.hartwig.hmftools.bamtools.remapper.HlaRecordAligner.mergeFlags;

import static htsjdk.samtools.SAMFlag.FIRST_OF_PAIR;
import static htsjdk.samtools.SAMFlag.MATE_REVERSE_STRAND;
import static htsjdk.samtools.SAMFlag.MATE_UNMAPPED;
import static htsjdk.samtools.SAMFlag.PROPER_PAIR;
import static htsjdk.samtools.SAMFlag.READ_PAIRED;
import static htsjdk.samtools.SAMFlag.READ_REVERSE_STRAND;
import static htsjdk.samtools.SAMFlag.READ_UNMAPPED;
import static htsjdk.samtools.SAMFlag.SECOND_OF_PAIR;
import static htsjdk.samtools.SAMFlag.SUPPLEMENTARY_ALIGNMENT;

import java.util.Set;

import javax.annotation.Nullable;

import org.broadinstitute.hellbender.utils.bwa.BwaMemAlignment;
import org.junit.Assert;
import org.junit.Test;

import htsjdk.samtools.SAMFlag;
import htsjdk.samtools.SAMRecord;

public class HlaRecordAlignerTest extends RemapperTestBase
{
    @Test
    public void mergeFlagsTest()
    {
        // These are the cases that are seen using COLO829v003R_AHHKYHDSXX_S13 as input
        // and GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.img as the reference genome for
        // the re-alignment of the HLAs.

        // 65, 16
        checkFlagMerge(
                Set.of(FIRST_OF_PAIR, READ_PAIRED),
                READ_REVERSE_STRAND,
                Set.of(FIRST_OF_PAIR, READ_REVERSE_STRAND, READ_PAIRED));

        // 73, 16
        checkFlagMerge(
                Set.of(FIRST_OF_PAIR, MATE_UNMAPPED, READ_PAIRED),
                READ_REVERSE_STRAND,
                Set.of(FIRST_OF_PAIR, READ_REVERSE_STRAND, MATE_UNMAPPED, READ_PAIRED));

        // 81, 16
        checkFlagMerge(
                Set.of(FIRST_OF_PAIR, READ_REVERSE_STRAND, READ_PAIRED),
                READ_REVERSE_STRAND,
                Set.of(FIRST_OF_PAIR, READ_REVERSE_STRAND, READ_PAIRED));

        // 83, 16
        checkFlagMerge(
                Set.of(FIRST_OF_PAIR, READ_REVERSE_STRAND, PROPER_PAIR, READ_PAIRED),
                READ_REVERSE_STRAND,
                Set.of(FIRST_OF_PAIR, READ_REVERSE_STRAND, PROPER_PAIR, READ_PAIRED));

        // 97, 16
        checkFlagMerge(
                Set.of(FIRST_OF_PAIR, MATE_REVERSE_STRAND, READ_PAIRED),
                READ_REVERSE_STRAND,
                Set.of(FIRST_OF_PAIR, MATE_REVERSE_STRAND, READ_REVERSE_STRAND, READ_PAIRED));

        // 99, 16
        checkFlagMerge(
                Set.of(FIRST_OF_PAIR, MATE_REVERSE_STRAND, PROPER_PAIR, READ_PAIRED),
                READ_REVERSE_STRAND,
                Set.of(FIRST_OF_PAIR, MATE_REVERSE_STRAND, READ_REVERSE_STRAND, PROPER_PAIR, READ_PAIRED));

        // 113, 16
        checkFlagMerge(
                Set.of(FIRST_OF_PAIR, MATE_REVERSE_STRAND, READ_REVERSE_STRAND, READ_PAIRED),
                READ_REVERSE_STRAND,
                Set.of(FIRST_OF_PAIR, MATE_REVERSE_STRAND, READ_REVERSE_STRAND, READ_PAIRED));

        // 121, 0
        checkFlagMerge(
                Set.of(FIRST_OF_PAIR, MATE_REVERSE_STRAND, READ_REVERSE_STRAND, MATE_UNMAPPED, READ_PAIRED),
                null,
                Set.of(FIRST_OF_PAIR, MATE_REVERSE_STRAND, MATE_UNMAPPED, READ_PAIRED));

        // 129, 0
        checkFlagMerge(
                Set.of(SECOND_OF_PAIR, READ_PAIRED),
                null,
                Set.of(SECOND_OF_PAIR, READ_PAIRED));

        // 133, 4
        checkFlagMerge(
                Set.of(SECOND_OF_PAIR, READ_UNMAPPED, READ_PAIRED),
                READ_UNMAPPED,
                Set.of(SECOND_OF_PAIR, READ_UNMAPPED, READ_PAIRED));

        // 145, 2048
        checkFlagMerge(
                Set.of(SECOND_OF_PAIR, READ_REVERSE_STRAND, READ_PAIRED),
                SUPPLEMENTARY_ALIGNMENT,
                Set.of(SUPPLEMENTARY_ALIGNMENT, SECOND_OF_PAIR, READ_PAIRED));

        // 147, 16
        checkFlagMerge(
                Set.of(SECOND_OF_PAIR, READ_REVERSE_STRAND, PROPER_PAIR, READ_PAIRED),
                READ_REVERSE_STRAND,
                Set.of(SECOND_OF_PAIR, READ_REVERSE_STRAND, PROPER_PAIR, READ_PAIRED));

        // 161, 0
        checkFlagMerge(
                Set.of(SECOND_OF_PAIR, MATE_REVERSE_STRAND, READ_PAIRED),
                null,
                Set.of(SECOND_OF_PAIR, MATE_REVERSE_STRAND, READ_PAIRED));

        // 163, 16
        checkFlagMerge(
                Set.of(SECOND_OF_PAIR, MATE_REVERSE_STRAND, PROPER_PAIR, READ_PAIRED),
                READ_REVERSE_STRAND,
                Set.of(SECOND_OF_PAIR, MATE_REVERSE_STRAND, READ_REVERSE_STRAND, PROPER_PAIR, READ_PAIRED));

        // 177, 0
        checkFlagMerge(
                Set.of(SECOND_OF_PAIR, MATE_REVERSE_STRAND, READ_REVERSE_STRAND, READ_PAIRED),
                null,
                Set.of(SECOND_OF_PAIR, MATE_REVERSE_STRAND, READ_PAIRED));

        // 181, 4
        checkFlagMerge(
                Set.of(SECOND_OF_PAIR, MATE_REVERSE_STRAND, READ_REVERSE_STRAND, READ_UNMAPPED, READ_PAIRED),
                READ_UNMAPPED,
                Set.of(SECOND_OF_PAIR, MATE_REVERSE_STRAND, READ_UNMAPPED, READ_PAIRED));
    }

    @Test
    public void mergeFlagsCornerCasesTest()
    {
        // Read unmapped flag.
        checkFlagMerge(
                Set.of(SECOND_OF_PAIR, READ_UNMAPPED, READ_PAIRED),
                null,
                Set.of(SECOND_OF_PAIR, READ_PAIRED));

        checkFlagMerge(
                Set.of(),
                READ_UNMAPPED,
                Set.of(READ_UNMAPPED));

        checkFlagMerge(
                Set.of(READ_UNMAPPED),
                READ_UNMAPPED,
                Set.of(READ_UNMAPPED));

        // Read reverse strand.
        checkFlagMerge(
                Set.of(FIRST_OF_PAIR, READ_REVERSE_STRAND, READ_PAIRED),
                null,
                Set.of(FIRST_OF_PAIR, READ_PAIRED));

        checkFlagMerge(
                Set.of(),
                READ_REVERSE_STRAND,
                Set.of(READ_REVERSE_STRAND));

        checkFlagMerge(
                Set.of(READ_REVERSE_STRAND),
                READ_REVERSE_STRAND,
                Set.of(READ_REVERSE_STRAND));

        // Mate reverse strand. TODO ??

        // Supplementary.
        checkFlagMerge(
                Set.of(SECOND_OF_PAIR, READ_PAIRED),
                SUPPLEMENTARY_ALIGNMENT,
                Set.of(SUPPLEMENTARY_ALIGNMENT, SECOND_OF_PAIR, READ_PAIRED));
    }

    private void checkFlagMerge(Set<SAMFlag> original, @Nullable SAMFlag alignerFlag, Set<SAMFlag> expected)
    {
        int originalFlag = original.stream().map(SAMFlag::intValue).reduce(0, Integer::sum);
        int alignerFlagValue = alignerFlag == null ? 0 : alignerFlag.intValue();
        int mergedFlag = mergeFlags(originalFlag, alignerFlagValue);
        int expectedFlag = expected.stream().map(SAMFlag::intValue).reduce(0, Integer::sum);
        Assert.assertEquals(expectedFlag, mergedFlag);
    }

    @Test
    public void createRemappedRecordTest()
    {
        SAMRecord record = records.get(12);
        BwaMemAlignment alignment =
                new BwaMemAlignment(16, 5, 31354375, 31354419, 0, 44, 60, 1, 39, 20, "44M107S", "22C1", null, -1, -1, 0);
        SAMRecord remappedRecord = HlaRecordAligner.createRemappedRecord(record, alignment);

        Assert.assertEquals(record.getReadName(), remappedRecord.getReadName());
        Assert.assertArrayEquals(record.getReadBases(), remappedRecord.getReadBases());
        Assert.assertArrayEquals(record.getBaseQualities(), remappedRecord.getBaseQualities());
        Assert.assertEquals("chr6", remappedRecord.getReferenceName());
        Assert.assertEquals(alignment.getRefStart(), remappedRecord.getAlignmentStart());
        Assert.assertEquals(alignment.getMapQual(), remappedRecord.getMappingQuality());
        //        Assert.assertEquals(alignment.getSamFlag(), remappedRecord.getFlags());
    }
}
