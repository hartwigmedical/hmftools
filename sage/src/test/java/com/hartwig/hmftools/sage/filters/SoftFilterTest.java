package com.hartwig.hmftools.sage.filters;

import static com.hartwig.hmftools.common.bam.SamRecordUtils.MATE_CIGAR_ATTRIBUTE;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.MockRefGenome.generateRandomBases;
import static com.hartwig.hmftools.common.test.MockRefGenome.getNextBase;
import static com.hartwig.hmftools.sage.SageConstants.MIN_CORE_DISTANCE;
import static com.hartwig.hmftools.sage.common.TestUtils.HIGH_QUAL_CONFIG;
import static com.hartwig.hmftools.sage.common.TestUtils.MOCK_REF_GENOME;
import static com.hartwig.hmftools.sage.common.TestUtils.MSI_JITTER_CALCS;
import static com.hartwig.hmftools.sage.common.TestUtils.READ_ID_GENERATOR;
import static com.hartwig.hmftools.sage.common.TestUtils.RECALIBRATION;
import static com.hartwig.hmftools.sage.common.TestUtils.TEST_CONFIG;
import static com.hartwig.hmftools.sage.common.TestUtils.buildCigarString;
import static com.hartwig.hmftools.sage.common.TestUtils.createSamRecord;
import static com.hartwig.hmftools.sage.common.VariantTier.HIGH_CONFIDENCE;
import static com.hartwig.hmftools.sage.common.VariantUtils.createReadContext;
import static com.hartwig.hmftools.sage.common.VariantUtils.createReadCounter;
import static com.hartwig.hmftools.sage.common.VariantUtils.createSageVariant;
import static com.hartwig.hmftools.sage.common.VariantUtils.createSimpleVariant;
import static com.hartwig.hmftools.sage.common.VariantUtils.sageVariantFromReadContextCounter;
import static com.hartwig.hmftools.sage.filter.SoftFilter.FRAGMENT_COORDS;
import static com.hartwig.hmftools.sage.filter.SoftFilter.MAX_EDGE_DISTANCE;
import static com.hartwig.hmftools.sage.filter.SoftFilter.MAX_GERMLINE_ALT_SUPPORT;
import static com.hartwig.hmftools.sage.filter.SoftFilter.MAX_GERMLINE_VAF;
import static com.hartwig.hmftools.sage.filter.SoftFilter.MIN_TUMOR_QUAL;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.Collections;
import java.util.List;

import com.hartwig.hmftools.sage.SageConfig;
import com.hartwig.hmftools.sage.candidate.Candidate;
import com.hartwig.hmftools.sage.common.VariantReadContext;
import com.hartwig.hmftools.sage.common.RefSequence;
import com.hartwig.hmftools.sage.common.SageVariant;
import com.hartwig.hmftools.sage.common.SimpleVariant;
import com.hartwig.hmftools.sage.common.VariantTier;
import com.hartwig.hmftools.sage.dedup.DedupMixedGermlineSomatic;
import com.hartwig.hmftools.sage.evidence.ReadContextCounter;
import com.hartwig.hmftools.sage.filter.VariantFilters;
import com.hartwig.hmftools.sage.quality.QualityCalculator;
import com.hartwig.hmftools.sage.sync.FragmentData;

import org.junit.Test;

import htsjdk.samtools.SAMRecord;

public class SoftFilterTest
{
    private static final String REF_BASES = "X" + generateRandomBases(100);

    private static final RefSequence REF_SEQUENCE = new RefSequence(1, REF_BASES.getBytes());

    private static final QualityCalculator QUALITY_CALCULATOR = new QualityCalculator(
            TEST_CONFIG, RECALIBRATION, REF_SEQUENCE, MOCK_REF_GENOME, MSI_JITTER_CALCS);

    private static final String TEST_READ_ID = "READ_01";
    private static final String TEST_CIGAR = "30M";

    private static final VariantFilters FILTERS = new VariantFilters(TEST_CONFIG);

    @Test
    public void testTumorQualFilter()
    {
        int position = 50;

        VariantReadContext readContext = createReadContext(
                createSimpleVariant(position, "A", "T"),
                REF_BASES.substring(48, 50), REF_BASES.substring(51, 53), REF_BASES.substring(38, 48), REF_BASES.substring(53, 63));

        ReadContextCounter readContextCounter = createReadCounter(readContext, HIGH_QUAL_CONFIG, VariantTier.PANEL);

        String altBase = readContextCounter.alt();

        String altReadBases = REF_BASES.substring(20, position) + altBase + REF_BASES.substring(position + 1, 80);
        String readCigar = buildCigarString(altReadBases.length());

        SAMRecord altRead = createSamRecord(READ_ID_GENERATOR.nextId(), CHR_1, 20, altReadBases, readCigar);

        SAMRecord refRead = createSamRecord(READ_ID_GENERATOR.nextId(), CHR_1, 20, REF_BASES.substring(20, 80), readCigar);

        readContextCounter.processRead(altRead, 1, null);
        readContextCounter.processRead(altRead, 1, null);
        readContextCounter.processRead(refRead, 1, null);
        readContextCounter.processRead(refRead, 1, null);

        SageVariant variant = sageVariantFromReadContextCounter(readContextCounter);
        FILTERS.applySoftFilters(variant);

        assertTrue(variant.filters().contains(MIN_TUMOR_QUAL));

        readContextCounter.processRead(altRead, 1, null);
        readContextCounter.processRead(altRead, 1, null);

        variant = sageVariantFromReadContextCounter(readContextCounter);
        FILTERS.applySoftFilters(variant);

        assertFalse(variant.filters().contains(MIN_TUMOR_QUAL));

        // TODO: add more variations to test quality site criteria
    }

    @Test
    public void testMaxReadEdgeDistanceFilter()
    {
        int position = 50;

        ReadContextCounter readContextCounter = createSnvReadContextCounter(position);

        // all reads have the variant near the end of the read, and so fail the max edge distance filter

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

        assertEquals(5, readContextCounter.readEdgeDistance().maxAltDistanceFromEdge());
        assertEquals(48, readContextCounter.readEdgeDistance().maxDistanceFromEdge());

        SageVariant variant = sageVariantFromReadContextCounter(readContextCounter);

        FILTERS.applySoftFilters(variant);

        assertTrue(variant.filters().contains(MAX_EDGE_DISTANCE));
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

        SageVariant variant = sageVariantFromReadContextCounter(readContextCounter);

        FILTERS.applySoftFilters(variant);

        assertTrue(variant.filters().contains(FRAGMENT_COORDS));

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

    @Test
    public void testSimpleAtGermlineVaf()
    {
        int position = 50;

        SimpleVariant variant = createSimpleVariant(position, "ACG", "TCA");
        VariantReadContext readContext = createReadContext(
                variant, REF_BASES.substring(48, position), REF_BASES.substring(53, 55), REF_BASES.substring(38, 48), REF_BASES.substring(55, 63));

        ReadContextCounter refCounter = createReadCounter(readContext, true);
        ReadContextCounter tumorCounter = createReadCounter(readContext);

        int readPosStart = 20;
        String altBases = variant.alt();

        // trigger alt in the germline and and full support tumor reads to observe the germline filter still being applied
        String altReadBases = REF_BASES.substring(readPosStart, position) + altBases + REF_BASES.substring(position + 3, 80);
        String readCigar = buildCigarString(altReadBases.length());

        SAMRecord altRead = createSamRecord(READ_ID_GENERATOR.nextId(), CHR_1, readPosStart, altReadBases, readCigar);
        tumorCounter.processRead(altRead, 1, null);

        String simpleAltReadBases = REF_BASES.substring(readPosStart, position - 1) + "G" + altBases + REF_BASES.substring(position + 3, 80);
        SAMRecord simpleAltRead = createSamRecord(READ_ID_GENERATOR.nextId(), CHR_1, 20, simpleAltReadBases, readCigar);
        refCounter.processRead(simpleAltRead, 1, null);

        assertEquals(1, refCounter.simpleAltMatches());
        assertEquals(0, refCounter.altSupport());

        Candidate candidate = new Candidate(HIGH_CONFIDENCE, readContext, 1, 1);
        SageVariant mnv = new SageVariant(candidate, List.of(refCounter), List.of(tumorCounter));

        FILTERS.applySoftFilters(mnv);

        assertTrue(mnv.filters().contains(MAX_GERMLINE_VAF));
    }

    private ReadContextCounter createSnvReadContextCounter(final int variantPosition)
    {
        String refBase = REF_BASES.substring(variantPosition, variantPosition + 1);
        String altBase = getNextBase(refBase);

        SimpleVariant variant = new SimpleVariant(CHR_1, variantPosition, refBase, altBase);

        String leftCore = REF_BASES.substring(variantPosition - MIN_CORE_DISTANCE, variantPosition);
        String rightCore = REF_BASES.substring(variantPosition + 1, variantPosition + 1 + MIN_CORE_DISTANCE);

        VariantReadContext readContext = createReadContext(variant, leftCore, rightCore);

        return createReadCounter(readContext, HIGH_QUAL_CONFIG, VariantTier.PANEL);
    }
}
