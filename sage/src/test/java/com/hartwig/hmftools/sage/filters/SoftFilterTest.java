package com.hartwig.hmftools.sage.filters;

import static com.hartwig.hmftools.common.samtools.SamRecordUtils.MATE_CIGAR_ATTRIBUTE;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.MockRefGenome.generateRandomBases;
import static com.hartwig.hmftools.common.test.MockRefGenome.getNextBase;
import static com.hartwig.hmftools.sage.common.TestUtils.RECALIBRATION;
import static com.hartwig.hmftools.sage.common.TestUtils.TEST_CONFIG;
import static com.hartwig.hmftools.sage.common.TestUtils.createReadContext;
import static com.hartwig.hmftools.sage.common.TestUtils.createSamRecord;
import static com.hartwig.hmftools.sage.common.TestUtils.createVariant;
import static com.hartwig.hmftools.sage.filter.SoftFilter.FRAGMENT_COORDS;
import static com.hartwig.hmftools.sage.filter.SoftFilter.MAX_EDGE_DISTANCE;
import static com.hartwig.hmftools.sage.filter.SoftFilter.MAX_GERMLINE_REL_RAW_BASE_QUAL;
import static com.hartwig.hmftools.sage.filter.SoftFilter.MAX_GERMLINE_VAF;
import static com.hartwig.hmftools.sage.filter.SoftFilter.MIN_TUMOR_QUAL;
import static com.hartwig.hmftools.sage.filter.SoftFilter.isGermlineAndNotTumorFiltered;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.sage.SageConfig;
import com.hartwig.hmftools.sage.common.IndexedBases;
import com.hartwig.hmftools.sage.common.ReadContext;
import com.hartwig.hmftools.sage.common.SageVariant;
import com.hartwig.hmftools.sage.common.SimpleVariant;
import com.hartwig.hmftools.sage.common.VariantTier;
import com.hartwig.hmftools.sage.evidence.ReadContextCounter;
import com.hartwig.hmftools.sage.filter.VariantFilters;
import com.hartwig.hmftools.sage.quality.QualityCalculator;
import com.hartwig.hmftools.sage.sync.FragmentData;

import org.apache.logging.log4j.util.Strings;
import org.junit.Test;

import htsjdk.samtools.SAMRecord;

public class SoftFilterTest
{
    private static final String REF_BASES = "X" + generateRandomBases(100);
    private static final IndexedBases REF_INDEXED_BASES = new IndexedBases(1, 0, REF_BASES.getBytes());
    private static final QualityCalculator QUALITY_CALCULATOR = new QualityCalculator(TEST_CONFIG.Quality, RECALIBRATION, REF_INDEXED_BASES);

    private static final String TEST_READ_ID = "READ_01";
    private static final String TEST_CIGAR = "30M";

    public static final SageConfig HIGH_QUAL_CONFIG = new SageConfig(true);

    private static final VariantFilters FILTERS = new VariantFilters(TEST_CONFIG);

    @Test
    public void testIsGermlineAndNotTumorFiltered()
    {
        assertTrue(isGermlineAndNotTumorFiltered(Sets.newHashSet(MAX_GERMLINE_REL_RAW_BASE_QUAL.filterName())));
        assertTrue(isGermlineAndNotTumorFiltered(Sets.newHashSet(MAX_GERMLINE_REL_RAW_BASE_QUAL.filterName(), MAX_GERMLINE_VAF.filterName())));

        assertFalse(isGermlineAndNotTumorFiltered(Sets.newHashSet()));
        assertFalse(isGermlineAndNotTumorFiltered(Sets.newHashSet(MAX_GERMLINE_REL_RAW_BASE_QUAL.filterName(), MIN_TUMOR_QUAL.filterName())));
        assertFalse(isGermlineAndNotTumorFiltered(Sets.newHashSet(MAX_GERMLINE_REL_RAW_BASE_QUAL.filterName(), "Random Filter")));
    }

    @Test
    public void testMaxReadEdgeDistanceFilter()
    {
        int position = 50;

        ReadContextCounter readContextCounter = createSnvReadContextCounter(position);

        // all reads have the variant near the end of the read, thereby failing the max edge distance filter

        String altBase = readContextCounter.alt();

        SAMRecord read = createSamRecord(
                TEST_READ_ID, CHR_1, 45, REF_BASES.substring(45, position) + altBase + REF_BASES.substring(position + 1, 75), TEST_CIGAR);

        readContextCounter.processRead(read, 1, null);
        readContextCounter.processRead(read, 1, null);
        readContextCounter.processRead(read, 1, null);

        SAMRecord read2 = createSamRecord(
                TEST_READ_ID, CHR_1, 45,  REF_BASES.substring(42, position) + altBase + REF_BASES.substring(position + 1, 70), TEST_CIGAR);

        readContextCounter.processRead(read2, 1, null);
        readContextCounter.processRead(read2, 1, null);

        // factor in soft-clipped bases
        SAMRecord read3 = createSamRecord(
                TEST_READ_ID, CHR_1, 46,
                REF_BASES.substring(43, position) + altBase + REF_BASES.substring(position + 1, 80), "3S30M4S");

        readContextCounter.processRead(read3, 1, null);

        // a read supporting the ref
        SAMRecord read4 = createSamRecord(TEST_READ_ID, CHR_1, 1, REF_BASES.substring(1, 99), "98M");
        readContextCounter.processRead(read4, 1, null);

        assertEquals(5, readContextCounter.readEdgeDistance().maxAltDistanceFromAlignedEdge());
        assertEquals(7, readContextCounter.readEdgeDistance().maxAltDistanceFromUnclippedEdge());
        assertEquals(48, readContextCounter.readEdgeDistance().maxDistanceFromUnclippedEdge());

        SageVariant variant = createVariant(readContextCounter);

        FILTERS.applySoftFilters(variant);

        assertTrue(variant.filters().contains(MAX_EDGE_DISTANCE.filterName()));
    }

    @Test
    public void testFragmentCoordsFilter()
    {
        // variants are filtered if the read coords supporting the alt are not sufficiently unique
        int position = 20;

        ReadContextCounter readContextCounter = createSnvReadContextCounter(position);
        String altBase = readContextCounter.alt();

        // 2 reads at the start, 3 at the end is insufficient
        String cigar1 = "2S25M3S";

        SAMRecord read1 = createSamRecord(
                TEST_READ_ID, CHR_1, 12, REF_BASES.substring(10, position) + altBase + REF_BASES.substring(position + 1, 40), cigar1);
        read1.setMateReferenceName(CHR_1);
        read1.setMateReferenceIndex(read1.getReferenceIndex());
        read1.setMateAlignmentStart(100);
        read1.setAttribute(MATE_CIGAR_ATTRIBUTE, cigar1);
        read1.setMateNegativeStrandFlag(true);

        // same frag coordinates processed twice
        readContextCounter.processRead(read1, 1, null);
        readContextCounter.processRead(read1, 1, null);

        assertEquals(1, readContextCounter.fragmentCoords().lowerCount());
        assertEquals(1, readContextCounter.fragmentCoords().upperCount());

        // same fragment start (unclipped), different mate end
        SAMRecord read2 = createSamRecord(
                TEST_READ_ID, CHR_1, 10,  REF_BASES.substring(10, position) + altBase + REF_BASES.substring(position + 1, 40), "30M");

        SAMRecord mate2 = createSamRecord(
                TEST_READ_ID, CHR_1, 14,  REF_BASES.substring(14, position) + altBase + REF_BASES.substring(position + 1, 44), "30M");

        mate2.setReadNegativeStrandFlag(true);

        readContextCounter.processRead(read2, 1, new FragmentData(read2, mate2));

        assertEquals(1, readContextCounter.fragmentCoords().lowerCount());
        assertEquals(2, readContextCounter.fragmentCoords().upperCount());

        SageVariant variant = createVariant(readContextCounter);

        FILTERS.applySoftFilters(variant);

        assertTrue(variant.filters().contains(FRAGMENT_COORDS.filterName()));

        // different frag coords pass the filter
        SAMRecord read3 = createSamRecord(
                TEST_READ_ID, CHR_1, 10,  REF_BASES.substring(8, position) + altBase + REF_BASES.substring(position + 1, 42), "2S30M");

        SAMRecord mate3 = createSamRecord(
                TEST_READ_ID, CHR_1, 14,  REF_BASES.substring(14, position) + altBase + REF_BASES.substring(position + 1, 46), "30M2S");

        mate3.setReadNegativeStrandFlag(true);

        readContextCounter.processRead(read3, 1, new FragmentData(read3, mate3));
        readContextCounter.processRead(read3, 1, new FragmentData(read3, mate3));

        assertEquals(2, readContextCounter.fragmentCoords().lowerCount());
        assertEquals(3, readContextCounter.fragmentCoords().upperCount());

        SAMRecord read4 = createSamRecord(
                TEST_READ_ID, CHR_1, 11, REF_BASES.substring(11, position) + altBase + REF_BASES.substring(position + 1, 40), "30M");
        read1.setMateReferenceName(CHR_1);
        read1.setMateReferenceIndex(read1.getReferenceIndex());
        read1.setMateAlignmentStart(100);
        read1.setAttribute(MATE_CIGAR_ATTRIBUTE, cigar1);
        read1.setMateNegativeStrandFlag(true);

        readContextCounter.processRead(read4, 1, null);

        assertTrue(readContextCounter.fragmentCoords().atCapacity());
    }

    private ReadContextCounter createSnvReadContextCounter(final int variantPosition)
    {
        String refBase = REF_BASES.substring(variantPosition, variantPosition + 1);
        String altBase = getNextBase(refBase);

        SimpleVariant variant = new SimpleVariant(CHR_1, variantPosition, refBase, altBase);

        // use minimum core with +/-2
        String readBases = REF_BASES.substring(
                variantPosition - 2, variantPosition) + altBase + REF_BASES.substring(variantPosition + 1, variantPosition + 3);

        final ReadContext readContext = createReadContext(variantPosition, 2, 0, 4, readBases, Strings.EMPTY);

        return new ReadContextCounter(
                1, variant, readContext, VariantTier.PANEL, 100, 0,
                HIGH_QUAL_CONFIG, QUALITY_CALCULATOR, null);
    }
}
