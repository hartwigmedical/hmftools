package com.hartwig.hmftools.compar.purple;

import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_CHROMOSOME;
import static com.hartwig.hmftools.compar.common.CommonUtils.FLD_CHROMOSOME_BAND;
import static com.hartwig.hmftools.compar.common.CommonUtils.FLD_REPORTED;
import static com.hartwig.hmftools.compar.purple.GermlineAmpDelData.FLD_GERMLINE_CN;
import static com.hartwig.hmftools.compar.purple.GermlineAmpDelData.FLD_GERMLINE_STATUS;
import static com.hartwig.hmftools.compar.purple.GermlineAmpDelData.FLD_TUMOR_CN;
import static com.hartwig.hmftools.compar.purple.GermlineAmpDelData.FLD_TUMOR_STATUS;

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
        extends ComparableItemTest<GermlineAmpDelData, GermlineAmpDelComparer, TestGermlineDeletionDataBuilder>
{
    @Before
    public void setUp()
    {
        comparer = new GermlineAmpDelComparer(new ComparConfig());
        builder = TestGermlineDeletionDataBuilder.BUILDER;
        GermlineAmpDelData alternateValueSource = builder.createWithAlternateDefaults();

        fieldToAlternateValueInitializer = Map.of(
            FLD_GERMLINE_STATUS, b -> b.germlineStatus = alternateValueSource.AmpDelData.NormalStatus,
            FLD_TUMOR_STATUS, b -> b.tumorStatus = alternateValueSource.AmpDelData.TumorStatus,
            FLD_GERMLINE_CN, b -> b.germlineCopyNumber = alternateValueSource.AmpDelData.GermlineCopyNumber,
            FLD_TUMOR_CN, b -> b.tumorCopyNumber = alternateValueSource.AmpDelData.TumorCopyNumber,
            FLD_CHROMOSOME, b -> b.comparisonChromosome = alternateValueSource.mComparisonChromosome,
            FLD_CHROMOSOME_BAND, b -> b.chromosomeBand = alternateValueSource.AmpDelData.ChromosomeBand
        );
        nameToAlternateIndexInitializer = Map.of("Gene", b -> b.gene = alternateValueSource.AmpDelData.GeneName);
        reportabilityFieldToFalseReportabilityInitializer = Map.of(FLD_REPORTED, b -> b.reported = false);
        nameToNonPassInitializer = Collections.emptyMap();
    }

    @Test
    public void fullyMatchesSelfWithLiftover()
    {
        GermlineAmpDelData victim = TestGermlineDeletionDataBuilder.BUILDER.create(b -> b.comparisonChromosome = "8");
        GermlineAmpDelData liftoverVictim = TestGermlineDeletionDataBuilder.BUILDER.create(b -> b.comparisonChromosome = "8");
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
