package com.hartwig.hmftools.compar.lilac;

import static com.hartwig.hmftools.common.hla.LilacQcData.FLD_DISC_ALIGN_FRAGS;
import static com.hartwig.hmftools.common.hla.LilacQcData.FLD_DISC_INDELS;
import static com.hartwig.hmftools.common.hla.LilacQcData.FLD_FIT_FRAGS;
import static com.hartwig.hmftools.common.hla.LilacQcData.FLD_HLA_Y;
import static com.hartwig.hmftools.common.hla.LilacQcData.FLD_QC_STATUS;
import static com.hartwig.hmftools.common.hla.LilacQcData.FLD_TOTAL_FRAGS;
import static com.hartwig.hmftools.compar.ComparTestUtil.assertSingleFieldMismatch;
import static com.hartwig.hmftools.compar.ComparTestUtil.assertValueDifferencesAsExpected;
import static com.hartwig.hmftools.compar.lilac.LilacData.FLD_ALLELES;

import static org.junit.Assert.assertNull;

import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.function.Consumer;

import com.hartwig.hmftools.common.hla.LilacAllele;
import com.hartwig.hmftools.compar.ComparConfig;
import com.hartwig.hmftools.compar.ComparableItemTest;
import com.hartwig.hmftools.compar.common.DiffThresholds;
import com.hartwig.hmftools.compar.common.MatchLevel;
import com.hartwig.hmftools.compar.common.MismatchType;

import org.junit.Before;
import org.junit.Test;

public class LilacDataTest extends ComparableItemTest<LilacData, LilacComparer, TestLilacDataBuilder>
{
    @Before
    public void setUp()
    {
        comparer = new LilacComparer(new ComparConfig());
        builder = TestLilacDataBuilder.BUILDER;

        List<LilacAllele> alleleSource = builder.create().Alleles;
        LilacData alternateValueSource = builder.createWithAlternateDefaults();

        // Does not include every field because field comparisons within alleles don't work well in generic tests
        fieldToAlternateValueInitializer = Map.of(
                FLD_QC_STATUS, b -> b.qcStatus = alternateValueSource.QcData.status(),
                FLD_TOTAL_FRAGS, b -> b.totalFragments = alternateValueSource.QcData.totalFragments(),
                FLD_FIT_FRAGS, b -> b.fittedFragments = alternateValueSource.QcData.fittedFragments(),
                FLD_DISC_ALIGN_FRAGS,
                b -> b.discardedAlignmentFragments = alternateValueSource.QcData.discardedAlignmentFragments(),
                FLD_DISC_INDELS, b -> b.discardedIndels = alternateValueSource.QcData.discardedIndels(),
                FLD_HLA_Y, b -> b.hlaYAllele = alternateValueSource.QcData.hlaYAllele(),
                FLD_ALLELES, b -> b.alleles = List.of(
                        alleleSource.get(0), alleleSource.get(1), alleleSource.get(2), alleleSource.get(3), alleleSource.get(4),
                        TestLilacAlleleBuilder.buildFrom(
                                alleleSource.get(5),
                                c -> c.allele = alternateValueSource.Alleles.get(5).allele()
                        )
                )
        );
        nameToAlternateIndexInitializer = Collections.emptyMap();
        reportabilityFieldToFalseReportabilityInitializer = Collections.emptyMap();
    }

    @Override
    @Test
    public void fullyDifferent()
    {
        // Overridden because field comparisons within alleles don't work well in generic test
        Set<String> expectedFieldNames = Set.of("A*01:01:RefTotal", "A*01:01:SomaticInframeIndel", "A*01:01:SomaticMissense",
                "A*01:01:SomaticNonsenseOrFrameshift", "A*01:01:SomaticSplice", "A*01:01:SomaticSynonymous", "A*01:01:TumorCopyNumber",
                "A*01:01:TumorTotal", "Alleles", "B*01:01:RefTotal", "B*01:01:SomaticInframeIndel", "B*01:01:SomaticMissense",
                "B*01:01:SomaticNonsenseOrFrameshift", "B*01:01:SomaticSplice", "B*01:01:SomaticSynonymous", "B*01:01:TumorCopyNumber",
                "B*01:01:TumorTotal", "B*01:02:RefTotal", "B*01:02:SomaticInframeIndel", "B*01:02:SomaticMissense",
                "B*01:02:SomaticNonsenseOrFrameshift", "B*01:02:SomaticSplice", "B*01:02:SomaticSynonymous", "B*01:02:TumorCopyNumber",
                "B*01:02:TumorTotal", "C*02:01:RefTotal", "C*02:01:SomaticInframeIndel", "C*02:01:SomaticMissense",
                "C*02:01:SomaticNonsenseOrFrameshift", "C*02:01:SomaticSplice", "C*02:01:SomaticSynonymous", "C*02:01:TumorCopyNumber",
                "C*02:01:TumorTotal", "DiscardedAlignmentFragments", "DiscardedIndels", "FittedFragments", "HlaYAllele", "Status",
                "TotalFragments");
        assertFullyDifferentAsExpected(expectedFieldNames, MatchLevel.DETAILED);
    }

    @Test
    public void fullyDifferentReportableLevel()
    {
        Set<String> expectedFieldNames = Set.of("A*01:01:SomaticInframeIndel", "A*01:01:SomaticMissense",
                "A*01:01:SomaticNonsenseOrFrameshift", "A*01:01:SomaticSplice", "A*01:01:TumorCopyNumber",
                "Alleles", "B*01:01:SomaticInframeIndel", "B*01:01:SomaticMissense",
                "B*01:01:SomaticNonsenseOrFrameshift", "B*01:01:SomaticSplice", "B*01:01:TumorCopyNumber",
                "B*01:02:SomaticInframeIndel", "B*01:02:SomaticMissense",
                "B*01:02:SomaticNonsenseOrFrameshift", "B*01:02:SomaticSplice", "B*01:02:TumorCopyNumber",
                "C*02:01:SomaticInframeIndel", "C*02:01:SomaticMissense",
                "C*02:01:SomaticNonsenseOrFrameshift", "C*02:01:SomaticSplice", "C*02:01:TumorCopyNumber",
                "DiscardedAlignmentFragments", "DiscardedIndels", "FittedFragments", "HlaYAllele", "Status",
                "TotalFragments");
        assertFullyDifferentAsExpected(expectedFieldNames, MatchLevel.REPORTABLE);
    }

    private void assertFullyDifferentAsExpected(final Set<String> expectedFieldNames, final MatchLevel matchLevel)
    {
        DiffThresholds diffThresholds = createDefaultThresholds();

        LilacData refVictim = builder.create();
        LilacData newVictim = builder.createWithAlternateDefaults();
        assertValueDifferencesAsExpected(refVictim, newVictim, matchLevel, diffThresholds, expectedFieldNames, true);
    }

    @Test
    public void singleFieldMismatchesInAlleleAreRecognized()
    {
        LilacAllele alternateValueSource = builder.createWithAlternateDefaults().Alleles.get(0);
        Map<String, Consumer<TestLilacAlleleBuilder>> fieldToAlternateAlleleValueInitializer = Map.of(
                "A*01:01:RefTotal", b -> b.refTotal = alternateValueSource.refFragments(),
                "A*01:01:SomaticInframeIndel", b -> b.inframeIndel = alternateValueSource.somaticInframeIndel(),
                "A*01:01:SomaticMissense", b -> b.missense = alternateValueSource.somaticMissense(),
                "A*01:01:SomaticNonsenseOrFrameshift", b -> b.nonsenseOrFrameshift = alternateValueSource.somaticNonsenseOrFrameshift(),
                "A*01:01:SomaticSplice", b -> b.splice = alternateValueSource.somaticSplice(),
                "A*01:01:SomaticSynonymous", b -> b.synonymous = alternateValueSource.somaticSynonymous(),
                "A*01:01:TumorCopyNumber",b -> b.tumorCopyNumber = alternateValueSource.tumorCopyNumber(),
                "A*01:01:TumorTotal", b -> b.tumorTotal = alternateValueSource.tumorFragments()
        );
        Set<String> detailedOnlyFields = Set.of("A*01:01:RefTotal", "A*01:01:SomaticSynonymous", "A*01:01:TumorTotal");

        DiffThresholds diffThresholds = createDefaultThresholds();
        for(Map.Entry<String, Consumer<TestLilacAlleleBuilder>> entry : fieldToAlternateAlleleValueInitializer.entrySet())
        {
            String field = entry.getKey();
            LilacData refVictim = builder.create();
            LilacData newVictim = builder.create(b -> b.alleles = List.of(
                    TestLilacAlleleBuilder.buildFrom(refVictim.Alleles.get(0), entry.getValue()),
                    refVictim.Alleles.get(1),
                    refVictim.Alleles.get(2),
                    refVictim.Alleles.get(3),
                    refVictim.Alleles.get(4),
                    refVictim.Alleles.get(5)
                )
            );

            assertSingleFieldMismatch(field, refVictim, newVictim, MatchLevel.DETAILED, diffThresholds, MismatchType.VALUE);
            if(detailedOnlyFields.contains(field))
            {
                assertNull("", refVictim.findMismatch(newVictim, MatchLevel.REPORTABLE, diffThresholds, false));
            }
            else
            {
                assertSingleFieldMismatch(field, refVictim, newVictim, MatchLevel.REPORTABLE, diffThresholds, MismatchType.VALUE);
            }
        }
    }


}
