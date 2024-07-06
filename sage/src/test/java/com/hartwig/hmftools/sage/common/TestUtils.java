package com.hartwig.hmftools.sage.common;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.bam.SamRecordUtils.MATE_CIGAR_ATTRIBUTE;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.NUM_MUTATONS_ATTRIBUTE;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.SamRecordTestUtils.DEFAULT_MAP_QUAL;
import static com.hartwig.hmftools.common.test.SamRecordTestUtils.buildDefaultBaseQuals;

import java.util.Collections;
import java.util.List;

import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.common.test.MockRefGenome;
import com.hartwig.hmftools.common.test.ReadIdGenerator;
import com.hartwig.hmftools.sage.SageConfig;
import com.hartwig.hmftools.sage.quality.MsiJitterCalcs;
import com.hartwig.hmftools.sage.quality.QualityCalculator;
import com.hartwig.hmftools.sage.bqr.BqrRecordMap;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordSetBuilder;

public class TestUtils
{
    public static final SageConfig TEST_CONFIG = createSageConfig();
    public static final SageConfig HIGH_QUAL_CONFIG = new SageConfig(true);

    public static final String TEST_SAMPLE = "SAMPLE";

    public static final BqrRecordMap RECALIBRATION = new BqrRecordMap(Collections.emptyList());

    public static final MsiJitterCalcs MSI_JITTER_CALCS = new MsiJitterCalcs();

    public static final MockRefGenome MOCK_REF_GENOME = new MockRefGenome();

    public static final ReadIdGenerator READ_ID_GENERATOR = new ReadIdGenerator();

    public static final String REF_BASES_200 =
        //             10        20        30        40        50        60        70        80        90
        //   0123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789
            "CGCAATATTCGGGTGGGAGTGACCCGATTTTCCAGGTGCGTTCGTCACCGCTGTCTGTGACTCGGAAAAAAAACTCCCTGACCCCTTGCGCTTCCCAGGT"
          + "GAGGCAATGCCTCGCCCTGCTTCGGCTCGCGCACAGTGCGCGCTACACACACTGGCCTGCGCCCACTGTCTGGCACTCCCTAGTGAGATGAACCCGGTAC";

    public static final RefSequence REF_SEQUENCE_200 = new RefSequence(0, REF_BASES_200.getBytes()); // note zero-based to line up with indices

    public static final QualityCalculator QUALITY_CALCULATOR = new QualityCalculator(
            TEST_CONFIG, RECALIBRATION, REF_SEQUENCE_200, MOCK_REF_GENOME, MSI_JITTER_CALCS);

    public static SageConfig createSageConfig()
    {
        // add input arguments as necessary or take the defaults
        return new SageConfig(false);
    }

    public static void setTumorQuality(final SageVariant variant, int count, int quality)
    {
        variant.tumorReadCounters().get(0).readSupportCounts().Full = count;
        variant.tumorReadCounters().get(0).readSupportQualityCounts().Full = quality;
    }

    public static String buildCigarString(int alignedLength) { return format("%dM", alignedLength); }

    public static String buildCigarString(int alignedLength, int leftSoftClip, int rightSoftClip)
    {
        StringBuilder sb = new StringBuilder();

        if(leftSoftClip > 0)
            sb.append(format("%dS", leftSoftClip));

        sb.append(format("%dM", alignedLength));

        if(rightSoftClip > 0)
            sb.append(format("%dS", rightSoftClip));

        return sb.toString();
    }

    public static void addLocalPhaseSet(final SageVariant variant, int lps, int readCount)
    {
        variant.tumorReadCounters().get(0).addLocalPhaseSet(lps, readCount, 0);
    }

    public static void clearFilters(final List<SageVariant> variants)
    {
        variants.forEach(x -> x.filters().clear());
    }

    public static void clearFilters(final SageVariant variant)
    {
        variant.filters().clear();
    }

    public static SAMRecord createSamRecord(
            final String readId, final String chrStr, int readStart, final String readBases, final String cigar)
    {
        SAMRecordSetBuilder recordBuilder = new SAMRecordSetBuilder();
        recordBuilder.setUnmappedHasBasesAndQualities(false);

        HumanChromosome chromosome = HumanChromosome.fromString(chrStr);

        SAMRecord record = recordBuilder.addFrag(
                readId, chromosome.ordinal(), readStart, false, false,
                cigar, readBases, 37, false);

        record.setReadBases(readBases.getBytes());

        final byte[] qualities = new byte[readBases.length()];

        for(int i = 0; i < readBases.length(); ++i)
            qualities[i] = 37;

        record.setBaseQualities(qualities);
        record.setReferenceName(chrStr);
        record.setReferenceIndex(chromosome.ordinal()); // need to override since no header is present

        record.setMateReferenceName(chrStr);
        record.setMateReferenceIndex(chromosome.ordinal());
        record.setMateAlignmentStart(readStart + 300);
        record.setMateNegativeStrandFlag(true);
        record.setMappingQuality(DEFAULT_MAP_QUAL);
        record.setAttribute(MATE_CIGAR_ATTRIBUTE, cigar);

        // to be correct this should match the cigar element count
        record.setAttribute(NUM_MUTATONS_ATTRIBUTE, 1);
        record.setFirstOfPairFlag(true);

        record.setReadPairedFlag(true);
        record.setProperPairFlag(true);

        return record;
    }

    public static SAMRecord buildSamRecord(final int alignmentStart, final String cigar, final String readBases, final String qualities)
    {
        return buildSamRecord(alignmentStart, cigar, readBases, qualities.getBytes());
    }

    public static SAMRecord buildSamRecord(final int alignmentStart, final String cigar, final String readString)
    {
        return buildSamRecord(alignmentStart, cigar, readString, buildDefaultBaseQuals(readString.length()));
    }

    public static SAMRecord buildSamRecord(final int alignmentStart, final String cigar, final String readBases, final byte[] qualities)
    {
        final SAMRecord record = new SAMRecord(null);
        record.setReadName(READ_ID_GENERATOR.nextId());
        record.setReferenceName(CHR_1);
        record.setAlignmentStart(alignmentStart);
        record.setCigarString(cigar);
        record.setReadString(readBases);
        record.setReadNegativeStrandFlag(false);
        record.setBaseQualities(qualities);
        record.setMappingQuality(DEFAULT_MAP_QUAL);
        record.setDuplicateReadFlag(false);
        record.setReadUnmappedFlag(false);
        record.setProperPairFlag(true);
        record.setReadPairedFlag(true);
        record.setInferredInsertSize(600);
        record.setAttribute(NUM_MUTATONS_ATTRIBUTE, 0);
        return record;
    }

    public static BaseRegion region(int start, int end)
    {
        return new BaseRegion(start, end);
    }
}
