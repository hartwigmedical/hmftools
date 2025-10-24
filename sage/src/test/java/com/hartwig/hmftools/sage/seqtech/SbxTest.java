package com.hartwig.hmftools.sage.seqtech;

import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.GeneTestUtils.getCodingBases;
import static com.hartwig.hmftools.common.test.MockRefGenome.getNextBase;
import static com.hartwig.hmftools.common.test.SamRecordTestUtils.DEFAULT_BASE_QUAL;
import static com.hartwig.hmftools.common.variant.VariantTier.LOW_CONFIDENCE;
import static com.hartwig.hmftools.sage.SageConstants.DEFAULT_FLANK_LENGTH;
import static com.hartwig.hmftools.sage.SageConstants.DEFAULT_MAX_READ_DEPTH;
import static com.hartwig.hmftools.sage.common.TestUtils.QUALITY_CALCULATOR;
import static com.hartwig.hmftools.sage.common.TestUtils.READ_ID_GENERATOR;
import static com.hartwig.hmftools.sage.common.TestUtils.REF_BASES_200;
import static com.hartwig.hmftools.sage.common.TestUtils.REF_SEQUENCE_200;
import static com.hartwig.hmftools.sage.common.TestUtils.TEST_CONFIG;
import static com.hartwig.hmftools.sage.common.TestUtils.TEST_SAMPLE;
import static com.hartwig.hmftools.sage.common.TestUtils.buildCigarString;
import static com.hartwig.hmftools.sage.common.TestUtils.createSamRecord;
import static com.hartwig.hmftools.sage.common.TestUtils.setIlluminaSequencing;
import static com.hartwig.hmftools.sage.common.TestUtils.setSbxSequencing;
import static com.hartwig.hmftools.sage.common.VariantUtils.createSimpleVariant;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import com.hartwig.hmftools.common.redux.BaseQualAdjustment;
import com.hartwig.hmftools.common.variant.SimpleVariant;
import com.hartwig.hmftools.sage.common.VariantReadContext;
import com.hartwig.hmftools.sage.common.VariantReadContextBuilder;
import com.hartwig.hmftools.sage.evidence.ReadContextCounter;

import org.junit.After;
import org.junit.Test;

import htsjdk.samtools.SAMRecord;

public class SbxTest
{
    public SbxTest()
    {
        setSbxSequencing();
    }

    @After
    public void resetSequencingType() { setIlluminaSequencing(); }

    @Test
    public void testUncertainCoreBases()
    {
        int position = 40;

        String refBase = REF_BASES_200.substring(position, position + 1);
        String altBase = getNextBase(refBase);
        SimpleVariant variant = createSimpleVariant(position, refBase, altBase);

        String readBases = REF_BASES_200.substring(20, position) + altBase + REF_BASES_200.substring(position + 1, 60);
        String readCigar = buildCigarString(readBases.length());

        SAMRecord read = createSamRecord(READ_ID_GENERATOR.nextId(), CHR_1, 20, readBases, readCigar);

        VariantReadContextBuilder builder = new VariantReadContextBuilder(DEFAULT_FLANK_LENGTH);

        VariantReadContext readContext = builder.createContext(variant, read, 20, REF_SEQUENCE_200);

        assertEquals(12, readContext.VarIndex);
        assertTrue(readContext.isValid());
        // assertEquals(readBases.substring(33, 49), readContext.coreStr());

        ReadContextCounter readCounter = new ReadContextCounter(
                0, readContext, LOW_CONFIDENCE,
                DEFAULT_MAX_READ_DEPTH, 1, TEST_CONFIG, QUALITY_CALCULATOR, TEST_SAMPLE, true);

        readCounter.processRead(read, 1, null);

        assertEquals(1, readCounter.readCounts().Full);

        // now a read with uncertain bases in the core - these are tracked
        read.getBaseQualities()[19] = BaseQualAdjustment.BASE_QUAL_MINIMUM;

        readCounter.processRead(read, 1, null);

        assertEquals(2, readCounter.readCounts().Full);
        assertEquals(0, readCounter.uncertainCoreBaseCount());

        read.getBaseQualities()[19] = DEFAULT_BASE_QUAL;
        read.getBaseQualities()[20] = BaseQualAdjustment.BASE_QUAL_MINIMUM;

        readCounter.processRead(read, 1, null);

        assertEquals(2, readCounter.readCounts().Full);
        assertEquals(1, readCounter.uncertainCoreBaseCount());
    }

}
