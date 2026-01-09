package com.hartwig.hmftools.compar.purple;

import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_CHROMOSOME;
import static com.hartwig.hmftools.compar.common.CommonUtils.FLD_CHROMOSOME_BAND;
import static com.hartwig.hmftools.compar.common.CommonUtils.FLD_REPORTED;
import static com.hartwig.hmftools.compar.purple.GermlineDeletionData.FLD_GERMLINE_CN;
import static com.hartwig.hmftools.compar.purple.GermlineDeletionData.FLD_GERMLINE_STATUS;
import static com.hartwig.hmftools.compar.purple.GermlineDeletionData.FLD_TUMOR_CN;
import static com.hartwig.hmftools.compar.purple.GermlineDeletionData.FLD_TUMOR_STATUS;

import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import static junit.framework.TestCase.assertEquals;

import java.util.Collections;
import java.util.Map;

import com.hartwig.hmftools.compar.ComparConfig;
import com.hartwig.hmftools.compar.ComparableItemTest;
import com.hartwig.hmftools.compar.common.DiffThresholds;
import com.hartwig.hmftools.compar.common.MatchLevel;
import com.hartwig.hmftools.compar.common.Mismatch;
import com.hartwig.hmftools.compar.common.MismatchType;

import org.junit.Before;
import org.junit.Test;

public class GermlineAmpDelDataTest
        extends ComparableItemTest<GermlineDeletionData, GermlineDeletionComparer, TestGermlineDeletionDataBuilder>
{
    @Before
    public void setUp()
    {
        comparer = new GermlineDeletionComparer(new ComparConfig());
        builder = TestGermlineDeletionDataBuilder.BUILDER;
        GermlineDeletionData alternateValueSource = builder.createWithAlternateDefaults();

        fieldToAlternateValueInitializer = Map.of(
            FLD_GERMLINE_STATUS, b -> b.germlineStatus = alternateValueSource.Deletion.NormalStatus,
            FLD_TUMOR_STATUS, b -> b.tumorStatus = alternateValueSource.Deletion.TumorStatus,
            FLD_GERMLINE_CN, b -> b.germlineCopyNumber = alternateValueSource.Deletion.GermlineCopyNumber,
            FLD_TUMOR_CN, b -> b.tumorCopyNumber = alternateValueSource.Deletion.TumorCopyNumber,
            FLD_CHROMOSOME, b -> b.comparisonChromosome = alternateValueSource.mComparisonChromosome,
            FLD_CHROMOSOME_BAND, b -> b.chromosomeBand = alternateValueSource.Deletion.ChromosomeBand
        );
        nameToAlternateIndexInitializer = Map.of("Gene", b -> b.gene = alternateValueSource.Deletion.GeneName);
        reportabilityFieldToFalseReportabilityInitializer = Map.of(FLD_REPORTED, b -> b.reported = false);
    }

    @Test
    public void fullyMatchesSelfWithLiftover()
    {
        GermlineDeletionData victim = TestGermlineDeletionDataBuilder.BUILDER.create(b -> b.comparisonChromosome = "8");
        GermlineDeletionData liftoverVictim = TestGermlineDeletionDataBuilder.BUILDER.create(b -> b.comparisonChromosome = "8");
        DiffThresholds diffThresholds = createDefaultThresholds();

        assertTrue(victim.matches(liftoverVictim));
        assertTrue(liftoverVictim.matches(victim));
        assertNull(victim.findMismatch(liftoverVictim, MatchLevel.DETAILED, diffThresholds, false));
        assertNull(victim.findMismatch(liftoverVictim, MatchLevel.REPORTABLE, diffThresholds, false));

        Mismatch expectedMatch = new Mismatch(victim, liftoverVictim, MismatchType.FULL_MATCH, Collections.emptyList());
        assertEquals(expectedMatch, victim.findMismatch(liftoverVictim, MatchLevel.DETAILED, diffThresholds, true));
        assertEquals(expectedMatch, victim.findMismatch(liftoverVictim, MatchLevel.REPORTABLE, diffThresholds, true));
    }
}
