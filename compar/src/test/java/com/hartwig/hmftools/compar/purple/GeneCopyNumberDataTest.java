package com.hartwig.hmftools.compar.purple;

import static org.junit.Assert.assertEquals;

import java.util.Collections;
import java.util.Map;

import com.hartwig.hmftools.compar.ComparConfig;
import com.hartwig.hmftools.compar.ComparableItemTest;
import com.hartwig.hmftools.compar.common.MatchLevel;
import com.hartwig.hmftools.compar.common.MismatchType;

import org.junit.Before;
import org.junit.Test;

public class GeneCopyNumberDataTest extends ComparableItemTest<GeneCopyNumberData, GeneCopyNumberComparer, TestGeneCopyNumberDataBuilder>
{
    @Before
    public void setUp()
    {
        comparer = new GeneCopyNumberComparer(new ComparConfig(), Collections.emptyMap());
        builder = TestGeneCopyNumberDataBuilder.BUILDER;
        GeneCopyNumberData alternateValueSource = builder.createWithAlternateDefaults();

        fieldToAlternateValueInitializer = Map.of(
                GeneCopyNumberComparer.Fields.MinCopyNumber.toString(), b -> b.minCopyNumber = alternateValueSource.CopyNumber.minCopyNumber(),
                GeneCopyNumberComparer.Fields.MaxCopyNumber.toString(), b -> b.maxCopyNumber = alternateValueSource.CopyNumber.maxCopyNumber(),
                GeneCopyNumberComparer.Fields.DriverType.toString(), b -> b.driverType = alternateValueSource.CopyNumber.driverType()
        );
        nameToAlternateIndexInitializer = Map.of("Gene", b -> b.gene = alternateValueSource.CopyNumber.geneName());
        reportabilityFieldToFalseReportabilityInitializer = Collections.emptyMap();
        nameToNonPassInitializer = Collections.emptyMap();
    }

    @Test
    public void transcriptIsShownButNotCompared()
    {
        GeneCopyNumberData canonical = builder.create(b -> b.transcriptId = "ENST00000498124");
        GeneCopyNumberData otherTranscript = builder.create(b ->
        {
            b.transcriptId = "ENST00000579755";
            b.isCanonical = false;
        });

        assertEquals(
                "ENST00000579755",
                otherTranscript.fieldValues().get(GeneCopyNumberComparer.Fields.TranscriptId.toString()).displayValue());
        assertEquals("false", otherTranscript.fieldValues().get(GeneCopyNumberComparer.Fields.IsCanonical.toString()).displayValue());
        assertEquals(MismatchType.FULL_MATCH, canonical.findMismatch(comparer, otherTranscript, MatchLevel.DETAILED, true).Type);
    }
}
