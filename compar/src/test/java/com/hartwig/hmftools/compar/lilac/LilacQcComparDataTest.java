package com.hartwig.hmftools.compar.lilac;

import java.util.Collections;
import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.common.hla.LilacAllele;
import com.hartwig.hmftools.compar.ComparConfig;
import com.hartwig.hmftools.compar.ComparableItemTest;

import org.junit.Before;

public class LilacQcComparDataTest extends ComparableItemTest<LilacQcComparData, LilacQcComparer, TestLilacQcDataBuilder>
{
    @Before
    public void setUp()
    {
        comparer = new LilacQcComparer(new ComparConfig(), Collections.emptyMap());
        builder = TestLilacQcDataBuilder.BUILDER;

        List<LilacAllele> alleleSource = builder.create().Alleles;
        LilacQcComparData alternateValueSource = builder.createWithAlternateDefaults();

        // Does not include every field because field comparisons within alleles don't work well in generic tests
        fieldToAlternateValueInitializer = Map.of(
                LilacQcComparer.Fields.Status.toString(), b -> b.qcStatus = alternateValueSource.QcData.status(),
                LilacQcComparer.Fields.TotalFragments.toString(), b -> b.totalFragments = alternateValueSource.QcData.totalFragments(),
                LilacQcComparer.Fields.FittedFragments.toString(), b -> b.fittedFragments = alternateValueSource.QcData.fittedFragments(),
                LilacQcComparer.Fields.DiscardedAlignmentFragments.toString(),
                b -> b.discardedAlignmentFragments = alternateValueSource.QcData.discardedAlignmentFragments(),
                LilacQcComparer.Fields.DiscardedIndels.toString(), b -> b.discardedIndels = alternateValueSource.QcData.discardedIndels(),
                LilacQcComparer.Fields.HlaYAllele.toString(), b -> b.hlaYAllele = alternateValueSource.QcData.hlaYAllele(),
                LilacQcComparer.Fields.Alleles.toString(), b -> b.alleles = List.of(
                        alleleSource.get(0), alleleSource.get(1), alleleSource.get(2), alleleSource.get(3), alleleSource.get(4),
                        TestLilacAlleleDataBuilder.buildFrom(
                                alleleSource.get(5),
                                c -> c.allele = alternateValueSource.Alleles.get(5).allele()
                        ).Allele
                )
        );

        nameToAlternateIndexInitializer = Map.of(
                "genes", b -> b.genes = alternateValueSource.QcData.genes());

        reportabilityFieldToFalseReportabilityInitializer = Collections.emptyMap();
        nameToNonPassInitializer = Collections.emptyMap();
    }
}
