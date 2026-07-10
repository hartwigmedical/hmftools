package com.hartwig.hmftools.compar.lilac;

import static com.hartwig.hmftools.common.hla.LilacQcData.FLD_DISC_ALIGN_FRAGS;
import static com.hartwig.hmftools.common.hla.LilacQcData.FLD_DISC_INDELS;
import static com.hartwig.hmftools.common.hla.LilacQcData.FLD_FIT_FRAGS;
import static com.hartwig.hmftools.common.hla.LilacQcData.FLD_HLA_Y;
import static com.hartwig.hmftools.common.hla.LilacQcData.FLD_QC_STATUS;
import static com.hartwig.hmftools.common.hla.LilacQcData.FLD_TOTAL_FRAGS;
import static com.hartwig.hmftools.compar.lilac.LilacQcComparer.FLD_ALLELES;

import java.util.Collections;
import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.common.hla.LilacAllele;
import com.hartwig.hmftools.compar.ComparConfig;
import com.hartwig.hmftools.compar.ComparableItemTest;

import org.junit.Before;

public class LilacQcDataTest extends ComparableItemTest<LilacQcData, LilacQcComparer, TestLilacQcDataBuilder>
{
    @Before
    public void setUp()
    {
        comparer = new LilacQcComparer(new ComparConfig());
        builder = TestLilacQcDataBuilder.BUILDER;

        List<LilacAllele> alleleSource = builder.create().Alleles;
        LilacQcData alternateValueSource = builder.createWithAlternateDefaults();

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
