package com.hartwig.hmftools.esvee.assembly;

import static com.hartwig.hmftools.esvee.assembly.AssemblyTestUtils.setIlluminaSequencing;
import static com.hartwig.hmftools.esvee.assembly.AssemblyTestUtils.setSbxSequencing;
import static com.hartwig.hmftools.esvee.common.CommonUtils.isHigherBaseQualCategory;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import com.hartwig.hmftools.common.redux.BaseQualAdjustment;
import com.hartwig.hmftools.common.sequencing.SbxBamUtils;

import org.junit.After;
import org.junit.Test;

public class SbxAssemblyTest
{
    public SbxAssemblyTest()
    {
        setSbxSequencing();
    }

    @After
    public void resetSequencingType() { setIlluminaSequencing(); }

    @Test
    public void testBaseQualComparisons()
    {
        byte qual1 = SbxBamUtils.SBX_SIMPLEX_QUAL;
        byte qual2 = BaseQualAdjustment.LOW_BASE_QUAL_THRESHOLD;
        assertTrue(isHigherBaseQualCategory(qual1, qual2));

        qual1 = SbxBamUtils.SBX_SIMPLEX_QUAL;
        qual2 = SbxBamUtils.SBX_SIMPLEX_QUAL;
        assertFalse(isHigherBaseQualCategory(qual1, qual2));

        qual2 = SbxBamUtils.SBX_DUPLEX_QUAL;
        assertFalse(isHigherBaseQualCategory(qual1, qual2));
        assertTrue(isHigherBaseQualCategory(qual2, qual1));

        // assertFalse(isHigherBaseQualCategory(qual1, qual2));
    }


}
