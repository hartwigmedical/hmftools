package com.hartwig.hmftools.sage.common;

import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.MockRefGenome.generateRandomBases;
import static com.hartwig.hmftools.sage.SageConstants.DEFAULT_READ_CONTEXT_FLANK_SIZE;
import static com.hartwig.hmftools.sage.SageConstants.MIN_CORE_DISTANCE;
import static com.hartwig.hmftools.sage.evidence.ReadContextCounter.RC_FULL;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.hotspot.ImmutableVariantHotspotImpl;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.sage.candidate.Candidate;
import com.hartwig.hmftools.sage.evidence.ReadContextCounter;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordSetBuilder;

public class TestUtils
{
    public static SageVariant createVariant(int position, final String ref, final String alt)
    {
        String readBases = buildReadContextBases(alt);

        // LF          L   I   RC          RF
        // 0123456789  01  2  34  0123456789
        int leftCoreIndex = DEFAULT_READ_CONTEXT_FLANK_SIZE;
        int index = leftCoreIndex + MIN_CORE_DISTANCE;
        int rightCoreIndex = index + alt.length() - 1 + MIN_CORE_DISTANCE;

        IndexedBases indexBases = new IndexedBases(
                position, index, leftCoreIndex, rightCoreIndex, DEFAULT_READ_CONTEXT_FLANK_SIZE, readBases.getBytes());

        return createVariant(position, ref, alt, indexBases);
    }

    public static SageVariant createVariant(int position, final String ref, final String alt, final IndexedBases indexBases)
    {
        VariantHotspot variant = createVariantHotspot(position, ref, alt);

        ReadContext readContext = new ReadContext(position, "", 0, "", indexBases, false);

        ReadContextCounter readCounter =  new ReadContextCounter(0, variant, readContext, VariantTier.LOW_CONFIDENCE,
                100, 1, false);

        List<ReadContextCounter> tumorCounters = Lists.newArrayList(readCounter);

        Candidate candidate = new Candidate(
                VariantTier.HIGH_CONFIDENCE, variant, tumorCounters.get(0).readContext(),
                100, 1, 1, 0, 0);

        List<ReadContextCounter> normalCounters = Lists.newArrayList();

        return new SageVariant(candidate, normalCounters, tumorCounters);
    }

    public static void setTumorQuality(final SageVariant variant, int count, int quality)
    {
        variant.tumorReadCounters().get(0).counts()[RC_FULL] = count;
        variant.tumorReadCounters().get(0).quality()[RC_FULL] = quality;
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

    public static boolean hasFilter(final SageVariant variant, final String filter)
    {
        return variant.filters().contains(filter);
    }

    public static String buildReadContextBases(final String alt)
    {
        String flank = generateRandomBases(DEFAULT_READ_CONTEXT_FLANK_SIZE);
        String core = generateRandomBases(MIN_CORE_DISTANCE);
        return flank + core + alt + core + flank;
    }

    public static VariantHotspot createVariantHotspot(int position, final String ref, final String alt)
    {
        return ImmutableVariantHotspotImpl.builder().chromosome("1").position(position).ref(ref).alt(alt).build();
    }

    public static ReadContextCounter createReadCounter(
            int position, final String ref, final String alt,
            int index, int leftCoreIndex, int rightCoreIndex, int flankSize, final String readBases)
    {
        VariantHotspot variant = ImmutableVariantHotspotImpl.builder()
                .chromosome("1")
                .position(position)
                .ref(ref)
                .alt(alt).build();

        IndexedBases indexBases = new IndexedBases(position, index, leftCoreIndex, rightCoreIndex, flankSize, readBases.getBytes());
        ReadContext readContext = new ReadContext(position, "", 0, "", indexBases, false);

        return new ReadContextCounter(0, variant, readContext, VariantTier.LOW_CONFIDENCE,
                100, 1, false);
    }

    public static SAMRecord createSamRecord(
            final String readId, final String chromosome, int readStart, final String readBases, final String cigar)
    {
        SAMRecordSetBuilder recordBuilder = new SAMRecordSetBuilder();
        recordBuilder.setUnmappedHasBasesAndQualities(false);
        SAMRecord record = recordBuilder.addFrag(
                readId, 1, readStart, false, false, cigar, readBases, 37, false);

        record.setReadBases(readBases.getBytes());

        final byte[] qualities = new byte[readBases.length()];

        for(int i = 0; i < readBases.length(); ++i)
            qualities[i] = 37;

        record.setBaseQualities(qualities);
        record.setReferenceName(chromosome);

        // to be correct this should match the cigar element count
        record.setAttribute("NM", 1);
        record.setFirstOfPairFlag(true);

        record.setReadPairedFlag(true);
        record.setProperPairFlag(true);

        return record;
    }

}
