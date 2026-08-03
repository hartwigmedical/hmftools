package com.hartwig.hmftools.compar.lilac;

import static com.hartwig.hmftools.common.hla.LilacAllele.FLD_REF_TOTAL;
import static com.hartwig.hmftools.common.hla.LilacAllele.FLD_SYNON;
import static com.hartwig.hmftools.common.hla.LilacAllele.FLD_TUMOR_TOTAL;
import static com.hartwig.hmftools.compar.ComparTestUtil.assertSingleFieldMismatch;

import java.util.Collections;
import java.util.Map;
import java.util.Set;
import java.util.function.Consumer;

import com.hartwig.hmftools.compar.ComparConfig;
import com.hartwig.hmftools.compar.ComparableItemTest;
import com.hartwig.hmftools.compar.common.FieldCheckCache;
import com.hartwig.hmftools.compar.common.MatchLevel;
import com.hartwig.hmftools.compar.common.MismatchType;

import org.junit.Before;
import org.junit.Test;

public class LilacAlleleDataTest extends ComparableItemTest<LilacAlleleData, LilacAlleleComparer, TestLilacAlleleDataBuilder>
{
    public Set<String> detailedOnlyFields;

    @Before
    public void setUp()
    {
        comparer = new LilacAlleleComparer(new ComparConfig(), Collections.emptyMap());
        builder = TestLilacAlleleDataBuilder.BUILDER;

        LilacAlleleData alternateValueSource = builder.createWithAlternateDefaults();

        fieldToAlternateValueInitializer = Map.of(
                LilacAlleleComparer.Fields.SomaticMissense.toString(), b -> b.missense = alternateValueSource.Allele.somaticMissense(),
                LilacAlleleComparer.Fields.SomaticNonsenseOrFrameshift.toString(), b -> b.nonsenseOrFrameshift = alternateValueSource.Allele.somaticNonsenseOrFrameshift(),
                LilacAlleleComparer.Fields.SomaticSplice.toString(), b -> b.splice = alternateValueSource.Allele.somaticSplice(),
                LilacAlleleComparer.Fields.SomaticInframeIndel.toString(), b -> b.inframeIndel = alternateValueSource.Allele.somaticInframeIndel(),
                LilacAlleleComparer.Fields.TumorCopyNumber.toString(), b -> b.tumorCopyNumber = alternateValueSource.Allele.tumorCopyNumber(),
                LilacAlleleComparer.Fields.RefTotal.toString(), b -> b.refTotal = alternateValueSource.Allele.refFragments(),
                LilacAlleleComparer.Fields.TumorTotal.toString(),b -> b.tumorTotal = alternateValueSource.Allele.tumorFragments(),
                LilacAlleleComparer.Fields.SomaticSynonymous.toString(), b -> b.synonymous = alternateValueSource.Allele.somaticSynonymous()
        );
        detailedOnlyFields = Set.of(FLD_REF_TOTAL, FLD_TUMOR_TOTAL, FLD_SYNON);

        nameToAlternateIndexInitializer = Map.of(
                "genes", b -> b.genes = alternateValueSource.Allele.genes(),
                "allele", b -> b.allele = alternateValueSource.Allele.allele(),
                "index", b -> b.index = alternateValueSource.Index
        );

        reportabilityFieldToFalseReportabilityInitializer = Collections.emptyMap();
        nameToNonPassInitializer = Collections.emptyMap();
    }

    @Test
    @Override
    public void singleFieldMismatchesAreRecognizedInReportableMode()
    {
        // Not all fields are compared in reportable mode
        for(Map.Entry<String, Consumer<TestLilacAlleleDataBuilder>> entry : fieldToAlternateValueInitializer.entrySet())
        {
            if(!detailedOnlyFields.contains(entry.getKey()))
            {
                String field = entry.getKey();
                LilacAlleleData refVictim = builder.create();
                LilacAlleleData newVictim = builder.create(entry.getValue());

                MatchLevel matchLevel = MatchLevel.REPORTABLE;
                FieldCheckCache fieldConfig = createDefaultThresholds(matchLevel);

                assertSingleFieldMismatch(comparer, field, refVictim, newVictim, matchLevel, MismatchType.VALUE);
            }
        }
    }
}
