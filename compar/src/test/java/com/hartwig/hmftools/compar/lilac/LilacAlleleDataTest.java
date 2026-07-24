package com.hartwig.hmftools.compar.lilac;

import static com.hartwig.hmftools.common.hla.LilacAllele.FLD_INDEL;
import static com.hartwig.hmftools.common.hla.LilacAllele.FLD_MISSENSE;
import static com.hartwig.hmftools.common.hla.LilacAllele.FLD_NFS;
import static com.hartwig.hmftools.common.hla.LilacAllele.FLD_REF_TOTAL;
import static com.hartwig.hmftools.common.hla.LilacAllele.FLD_SPLICE;
import static com.hartwig.hmftools.common.hla.LilacAllele.FLD_SYNON;
import static com.hartwig.hmftools.common.hla.LilacAllele.FLD_TUMOR_CN;
import static com.hartwig.hmftools.common.hla.LilacAllele.FLD_TUMOR_TOTAL;
import static com.hartwig.hmftools.compar.ComparTestUtil.assertSingleFieldMismatch;

import java.util.Collections;
import java.util.Map;
import java.util.Set;
import java.util.function.Consumer;

import com.hartwig.hmftools.compar.ComparConfig;
import com.hartwig.hmftools.compar.ComparableItemTest;
import com.hartwig.hmftools.compar.common.FieldConfig;
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
        comparer = new LilacAlleleComparer(new ComparConfig());
        builder = TestLilacAlleleDataBuilder.BUILDER;

        LilacAlleleData alternateValueSource = builder.createWithAlternateDefaults();

        fieldToAlternateValueInitializer = Map.of(
                FLD_MISSENSE, b -> b.missense = alternateValueSource.Allele.somaticMissense(),
                FLD_NFS, b -> b.nonsenseOrFrameshift = alternateValueSource.Allele.somaticNonsenseOrFrameshift(),
                FLD_SPLICE, b -> b.splice = alternateValueSource.Allele.somaticSplice(),
                FLD_INDEL, b -> b.inframeIndel = alternateValueSource.Allele.somaticInframeIndel(),
                FLD_TUMOR_CN, b -> b.tumorCopyNumber = alternateValueSource.Allele.tumorCopyNumber(),
                FLD_REF_TOTAL, b -> b.refTotal = alternateValueSource.Allele.refFragments(),
                FLD_TUMOR_TOTAL,b -> b.tumorTotal = alternateValueSource.Allele.tumorFragments(),
                FLD_SYNON, b -> b.synonymous = alternateValueSource.Allele.somaticSynonymous()
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
                FieldConfig fieldConfig = createDefaultThresholds(matchLevel);

                assertSingleFieldMismatch(field, refVictim, newVictim, matchLevel, fieldConfig, MismatchType.VALUE);
            }
        }
    }
}
