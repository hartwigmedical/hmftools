package com.hartwig.hmftools.common.purple;

import static org.junit.Assert.assertNotNull;

import java.util.Collection;
import java.util.Optional;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.purple.copynumber.sv.ImmutableStructuralVariantLegPloidy;
import com.hartwig.hmftools.common.purple.region.ImmutableFittedRegion;
import com.hartwig.hmftools.common.purple.region.ObservedRegion;
import com.hartwig.hmftools.common.purple.segment.SegmentSupport;
import com.hartwig.hmftools.common.variant.structural.ImmutableStructuralVariantImpl;
import com.hartwig.hmftools.common.variant.structural.ImmutableStructuralVariantLegImpl;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantType;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.junit.Test;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;

public class PurpleDatamodelTest {

    private static final String CHROMOSOME = "1";

    @Test
    public void testDefaultFittedRegion() {
        assertNotNull(PurpleTestUtils.createDefaultFittedRegion(CHROMOSOME, 1, 100).build());
    }

    @Test
    public void testDefaultCopyNumber() {
        assertNotNull(PurpleTestUtils.createCopyNumber(CHROMOSOME, 1, 100, 2).build());
    }

}
