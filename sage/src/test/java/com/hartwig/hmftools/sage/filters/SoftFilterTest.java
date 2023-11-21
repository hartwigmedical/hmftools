package com.hartwig.hmftools.sage.filters;

import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.MockRefGenome.generateRandomBases;
import static com.hartwig.hmftools.common.test.MockRefGenome.getNextBase;
import static com.hartwig.hmftools.sage.common.TestUtils.RECALIBRATION;
import static com.hartwig.hmftools.sage.common.TestUtils.TEST_CONFIG;
import static com.hartwig.hmftools.sage.common.TestUtils.createReadContext;
import static com.hartwig.hmftools.sage.common.TestUtils.createSamRecord;
import static com.hartwig.hmftools.sage.common.TestUtils.createVariant;
import static com.hartwig.hmftools.sage.filter.SoftFilter.MAX_EDGE_DISTANCE;
import static com.hartwig.hmftools.sage.filter.SoftFilter.MAX_GERMLINE_REL_RAW_BASE_QUAL;
import static com.hartwig.hmftools.sage.filter.SoftFilter.MAX_GERMLINE_VAF;
import static com.hartwig.hmftools.sage.filter.SoftFilter.MIN_TUMOR_QUAL;
import static com.hartwig.hmftools.sage.filter.SoftFilter.isGermlineAndNotTumorFiltered;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.variant.hotspot.ImmutableVariantHotspotImpl;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.sage.common.IndexedBases;
import com.hartwig.hmftools.sage.common.ReadContext;
import com.hartwig.hmftools.sage.common.SageVariant;
import com.hartwig.hmftools.sage.common.VariantTier;
import com.hartwig.hmftools.sage.evidence.ReadContextCounter;
import com.hartwig.hmftools.sage.filter.VariantFilters;
import com.hartwig.hmftools.sage.quality.QualityCalculator;

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
        int position = 20;

        String refBase = REF_BASES.substring(position, position + 1);
        String altBase = getNextBase(refBase);

        final VariantHotspot hotspot = ImmutableVariantHotspotImpl.builder()
                .chromosome(CHR_1).ref(refBase).alt(altBase).position(position).build();

        String readBases = REF_BASES.substring(18, position) + altBase + REF_BASES.substring(position + 1, 23);
        final ReadContext readContext = createReadContext(position, 2, 0, 4, readBases, Strings.EMPTY);

        final ReadContextCounter readContextCounter = new ReadContextCounter(
                1, hotspot, readContext, VariantTier.PANEL, 100, 0,
                TEST_CONFIG, QUALITY_CALCULATOR, null);

        // all reads have the variant near the end of the read, thereby failing the max edge distance filter

        SAMRecord read = createSamRecord(
                TEST_READ_ID, CHR_1, 15, REF_BASES.substring(15, position) + altBase + REF_BASES.substring(position + 1, 45), TEST_CIGAR);

        readContextCounter.processRead(read, 1, null);
        readContextCounter.processRead(read, 1, null);
        readContextCounter.processRead(read, 1, null);

        SAMRecord read2 = createSamRecord(
                TEST_READ_ID, CHR_1, 5,  REF_BASES.substring(5, position) + altBase + REF_BASES.substring(position + 1, 35), TEST_CIGAR);

        readContextCounter.processRead(read2, 1, null);
        readContextCounter.processRead(read2, 1, null);

        assertEquals(14, readContextCounter.maxDistanceFromEdge());

        SageVariant variant = createVariant(readContextCounter);

        FILTERS.applySoftFilters(variant);

        assertTrue(variant.filters().contains(MAX_EDGE_DISTANCE.filterName()));
    }
}
