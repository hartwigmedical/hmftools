package com.hartwig.hmftools.markdups;

import static java.lang.Math.abs;

import static com.hartwig.hmftools.common.samtools.SamRecordUtils.MATE_CIGAR_ATTRIBUTE;
import static com.hartwig.hmftools.common.samtools.SamRecordUtils.SUPPLEMENTARY_ATTRIBUTE;
import static com.hartwig.hmftools.markdups.common.Constants.DEFAULT_PARTITION_SIZE;
import static com.hartwig.hmftools.markdups.common.Constants.DEFAULT_POS_BUFFER_SIZE;

import java.util.List;

import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeCoordinates;
import com.hartwig.hmftools.common.samtools.SupplementaryReadData;
import com.hartwig.hmftools.common.test.MockRefGenome;
import com.hartwig.hmftools.markdups.common.Fragment;
import com.hartwig.hmftools.markdups.common.FragmentStatus;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordSetBuilder;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;

public final class TestUtils
{
    public static final String TEST_READ_BASES = MockRefGenome.generateRandomBases(100);
    public static final String TEST_READ_ID = "READ_01";
    public static final String TEST_READ_CIGAR = "100M";

    public static final int DEFAULT_QUAL = 37;

    public static MarkDupsConfig createTestConfig()
    {
        return new MarkDupsConfig(DEFAULT_PARTITION_SIZE, DEFAULT_POS_BUFFER_SIZE, new MockRefGenome());
    }

    public static SAMSequenceDictionary SAM_DICTIONARY_V37;

    static
    {
        SAM_DICTIONARY_V37 = new SAMSequenceDictionary();

        RefGenomeCoordinates v37Coords = RefGenomeCoordinates.COORDS_37;

        for(HumanChromosome chromosome : HumanChromosome.values())
        {
            SAM_DICTIONARY_V37.addSequence(new SAMSequenceRecord(chromosome.toString(), v37Coords.Lengths.get(chromosome)));
        }
    }

    public static void resetFragment(final Fragment fragment) { fragment.setStatus(FragmentStatus.UNSET); }
    public static void resetFragments(final List<Fragment> fragments) { fragments.forEach(x -> x.setStatus(FragmentStatus.UNSET)); }

    public static Fragment createFragment(final String readId, final String chrStr, int readStart)
    {
        SAMRecord read = createSamRecord(readId, chrStr, readStart, TEST_READ_BASES, TEST_READ_CIGAR, chrStr, 200,
                false, false, null);
        return new Fragment(read);
    }

    public static Fragment createFragment(
            final String readId, final String chrStr, int readStart, final String readBases, final String cigar, final String mateChr,
            int mateStart, boolean isReversed, boolean isSupplementary, final SupplementaryReadData suppAlignment)
    {
        SAMRecord read = createSamRecord(readId, chrStr, readStart, readBases, cigar, mateChr, mateStart,
                isReversed, isSupplementary, suppAlignment);
        return new Fragment(read);
    }

    public static Fragment createFragment(
            final String readId, final String chrStr, int readStart, final String cigar, boolean isReversed,
            final String mateChr, int mateStart, boolean mateReversed, final String mateCigar)
    {
        SAMRecord read = createSamRecord(
                readId, chrStr, readStart, TEST_READ_BASES, cigar, mateChr, mateStart, isReversed, false, null);

        read.setAttribute(MATE_CIGAR_ATTRIBUTE, mateCigar);
        read.setMateNegativeStrandFlag(mateReversed);
        return new Fragment(read);
    }

    public static SAMRecord createSamRecord(
            final String readId, final String chrStr, int readStart, final String readBases, final String cigar, final String mateChr,
            int mateStart, boolean isReversed, boolean isSupplementary, final SupplementaryReadData suppAlignment)
    {
        SAMRecordSetBuilder recordBuilder = new SAMRecordSetBuilder();
        recordBuilder.getHeader().setSequenceDictionary(SAM_DICTIONARY_V37);
        recordBuilder.setUnmappedHasBasesAndQualities(false);

        HumanChromosome chromosome = HumanChromosome.fromString(chrStr);

        SAMRecord record = recordBuilder.addFrag(
                readId, chromosome.ordinal(), readStart, isReversed, false,
                cigar, readBases, DEFAULT_QUAL, false);

        record.setReadBases(readBases.getBytes());

        final byte[] qualities = new byte[readBases.length()];

        for(int i = 0; i < readBases.length(); ++i)
            qualities[i] = DEFAULT_QUAL;

        record.setBaseQualities(qualities);
        record.setReferenceName(chrStr);
        record.setReferenceIndex(chromosome.ordinal()); // need to override since no header is present

        if(!mateChr.isEmpty())
        {
            record.setMateReferenceName(mateChr);
            record.setMateAlignmentStart(mateStart);
            record.setMateReferenceIndex(HumanChromosome.fromString(mateChr).ordinal());
            record.setReadPairedFlag(true);
            record.setProperPairFlag(true);
        }
        else
        {
            record.setReadPairedFlag(false);
            record.setProperPairFlag(false);
        }

        // to be correct this should match the cigar element count
        record.setFirstOfPairFlag(true);

        record.setSupplementaryAlignmentFlag(isSupplementary);

        if(suppAlignment != null)
            record.setAttribute(SUPPLEMENTARY_ATTRIBUTE, suppAlignment.asCsv());

        if(chrStr.equals(mateChr))
            record.setInferredInsertSize(abs(readStart - mateStart));
        else
            record.setInferredInsertSize(0);

        return record;
    }

    public static void setBaseQualities(final Fragment fragment, int value)
    {
        fragment.reads().forEach(x -> setBaseQualities(x, value));
    }

    public static void setBaseQualities(final SAMRecord read, int value)
    {
        for(int i = 0; i < read.getBaseQualities().length; ++i)
            read.getBaseQualities()[i] = (byte)value;
    }

}
