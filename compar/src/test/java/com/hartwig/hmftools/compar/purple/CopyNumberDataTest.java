package com.hartwig.hmftools.compar.purple;

import static com.hartwig.hmftools.compar.purple.CopyNumberData.FLD_COPY_NUMBER;
import static com.hartwig.hmftools.compar.purple.CopyNumberData.FLD_MAJOR_ALLELE_CN;
import static com.hartwig.hmftools.compar.purple.CopyNumberData.FLD_METHOD;

import static org.junit.Assert.assertFalse;
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

public class CopyNumberDataTest extends ComparableItemTest<CopyNumberData, CopyNumberComparer, TestCopyNumberDataBuilder>
{
    @Before
    public void setUp()
    {
        comparer = new CopyNumberComparer(new ComparConfig());
        builder = TestCopyNumberDataBuilder.BUILDER;
        CopyNumberData alternateValueSource = builder.createWithAlternateDefaults();

        fieldToAlternateValueInitializer = Map.of(
                FLD_COPY_NUMBER, b -> b.copyNumber = alternateValueSource.copyNumber(),
                FLD_MAJOR_ALLELE_CN, b -> b.majorAlleleCopyNumber = alternateValueSource.majorAlleleCopyNumber(),
                FLD_METHOD, b -> b.method = alternateValueSource.method()
        );
        nameToAlternateIndexInitializer = Map.of(
                "Chromosome", b -> {
                    b.chromosome = alternateValueSource.chromosome();
                    b.comparisonChromosomeStart = alternateValueSource.comparisonPositionStart().Chromosome;
                    b.comparisonChromosomeEnd = alternateValueSource.comparisonPositionEnd().Chromosome;
                },
                "PositionStart", b -> {
                    b.positionStart = alternateValueSource.positionStart();
                    b.comparisonPositionStart = alternateValueSource.comparisonPositionStart().Position;
                },
                "PositionEnd", b -> {
                    b.positionEnd = alternateValueSource.positionEnd();
                    b.comparisonPositionEnd = alternateValueSource.comparisonPositionEnd().Position;
                }
        );
        reportabilityFieldToFalseReportabilityInitializer = Collections.emptyMap();
        nameToNonPassInitializer = Collections.emptyMap();
    }

    @Override
    @Test
    public void fullyMatchesSelfInReportableMode()
    {
        // Override since copy numbers are never compared in reportable mode
    }

    @Override
    @Test
    public void singleFieldMismatchesAreRecognizedInReportableMode()
    {
        // Override since copy numbers are never compared in reportable mode
    }

    @Test
    public void fullyMatchesSelfWithLiftover()
    {
        CopyNumberData victim = TestCopyNumberDataBuilder.BUILDER.create(b ->
        {
            b.comparisonChromosomeStart = "8";
            b.comparisonChromosomeEnd = "8";
            b.comparisonPositionStart = 15000;
            b.comparisonPositionEnd = 25000;
        });
        CopyNumberData liftoverVictim = TestCopyNumberDataBuilder.BUILDER.create(b ->
        {
            b.chromosome = "8";
            b.positionStart = 15000;
            b.positionEnd = 25000;
        });
        DiffThresholds diffThresholds = createDefaultThresholds();

        assertTrue(victim.matches(liftoverVictim));
        assertTrue(liftoverVictim.matches(victim));
        assertNull(victim.findMismatch(liftoverVictim, MatchLevel.DETAILED, diffThresholds, false));

        Mismatch expectedMatch = new Mismatch(victim, liftoverVictim, MismatchType.FULL_MATCH, Collections.emptyList());
        assertEquals(expectedMatch, victim.findMismatch(liftoverVictim, MatchLevel.DETAILED, diffThresholds, true));
    }

    @Test
    public void keyNonEmptyForLiftover()
    {
        CopyNumberData victim = TestCopyNumberDataBuilder.BUILDER.create(b ->
        {
            b.comparisonChromosomeStart = "8";
            b.comparisonChromosomeEnd = "8";
            b.comparisonPositionStart = 15000;
            b.comparisonPositionEnd = 25000;
        });
        assertFalse(victim.key().isEmpty());
    }
}
