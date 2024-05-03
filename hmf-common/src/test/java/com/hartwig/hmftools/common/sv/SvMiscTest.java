package com.hartwig.hmftools.common.sv;

import static com.hartwig.hmftools.common.genome.region.Orientation.FORWARD;
import static com.hartwig.hmftools.common.genome.region.Orientation.ORIENT_FWD;
import static com.hartwig.hmftools.common.genome.region.Orientation.ORIENT_NEG_CHAR;
import static com.hartwig.hmftools.common.genome.region.Orientation.ORIENT_POS_CHAR;
import static com.hartwig.hmftools.common.genome.region.Orientation.ORIENT_REV;
import static com.hartwig.hmftools.common.genome.region.Orientation.REVERSE;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.genome.region.Orientation;

import org.jetbrains.annotations.NotNull;
import org.junit.Before;
import org.junit.Test;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.filter.CompoundFilter;
import htsjdk.variant.vcf.VCFCodec;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderVersion;

public class SvMiscTest
{
    @Test
    public void testOrientation()
    {
        assertEquals(FORWARD, REVERSE.opposite());
        assertEquals(REVERSE, FORWARD.opposite());

        assertTrue(FORWARD.equalsByte(ORIENT_FWD));
        assertTrue(REVERSE.equalsByte(ORIENT_REV));

        assertEquals(FORWARD, Orientation.fromByte(ORIENT_FWD));
        assertEquals(REVERSE, Orientation.fromByte(ORIENT_REV));

        assertEquals(FORWARD, Orientation.fromChar(ORIENT_POS_CHAR));
        assertEquals(REVERSE, Orientation.fromChar(ORIENT_NEG_CHAR));

        assertEquals(FORWARD, Orientation.fromByteStr(String.valueOf(ORIENT_FWD)));
        assertEquals(REVERSE, Orientation.fromByteStr(String.valueOf(ORIENT_REV)));
    }

    @Test
    public void testSvBreakendVcfFormation()
    {

    }

}
