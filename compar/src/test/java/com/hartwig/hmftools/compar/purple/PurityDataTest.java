package com.hartwig.hmftools.compar.purple;

import java.util.Collections;
import java.util.HashMap;

import com.hartwig.hmftools.compar.ComparConfig;
import com.hartwig.hmftools.compar.ComparableItemTest;

import org.junit.Before;

public class PurityDataTest extends ComparableItemTest<PurityData, PurityComparer, TestPurityDataBuilder>
{
    @Before
    public void setUp()
    {
        comparer = new PurityComparer(new ComparConfig(), Collections.emptyMap());
        builder = TestPurityDataBuilder.BUILDER;
        PurityData alternateValueSource = builder.createWithAlternateDefaults();

        fieldToAlternateValueInitializer = new HashMap<>();
        fieldToAlternateValueInitializer.put(PurityComparer.PurityFields.Purity.toString(), b -> b.purity = alternateValueSource.Purity.bestFit().purity());
        fieldToAlternateValueInitializer.put(PurityComparer.PurityFields.Ploidy.toString(), b -> b.ploidy = alternateValueSource.Purity.bestFit().ploidy());
        fieldToAlternateValueInitializer.put(PurityComparer.PurityFields.Contamination.toString(), b -> b.contamination = alternateValueSource.Purity.qc().contamination());
        fieldToAlternateValueInitializer.put(PurityComparer.PurityFields.TmbPerMb.toString(), b -> b.tmb = alternateValueSource.Purity.tumorMutationalBurdenPerMb());
        fieldToAlternateValueInitializer.put(PurityComparer.PurityFields.MsIndelsPerMb.toString(), b -> b.msIndels = alternateValueSource.Purity.microsatelliteIndelsPerMb());
        fieldToAlternateValueInitializer.put(PurityComparer.PurityFields.Tml.toString(), b -> b.tml = alternateValueSource.Purity.tumorMutationalLoad());
        fieldToAlternateValueInitializer.put(PurityComparer.PurityFields.CopyNumberSegments.toString(), b -> b.copyNumberSegments = alternateValueSource.Purity.qc().copyNumberSegments());
        fieldToAlternateValueInitializer.put(PurityComparer.PurityFields.UnsupportedCopyNumberSegments.toString(), b -> b.unsupportedCopyNumberSegments =
                alternateValueSource.Purity.qc().unsupportedCopyNumberSegments());
        fieldToAlternateValueInitializer.put(PurityComparer.PurityFields.SvTmb.toString(), b -> b.svTmb = alternateValueSource.Purity.svTumorMutationalBurden());
        fieldToAlternateValueInitializer.put(PurityComparer.PurityFields.QcStatus.toString(), b -> b.qcStatus = alternateValueSource.Purity.qc().status());
        fieldToAlternateValueInitializer.put(PurityComparer.PurityFields.Gender.toString(), b -> b.gender = alternateValueSource.Purity.gender());
        fieldToAlternateValueInitializer.put(PurityComparer.PurityFields.GermlineAberrations.toString(), b -> b.germlineAberrations = alternateValueSource.Purity.qc()
                .germlineAberrations());
        fieldToAlternateValueInitializer.put(PurityComparer.PurityFields.FitMethod.toString(), b -> b.fitMethod = alternateValueSource.Purity.method());
        fieldToAlternateValueInitializer.put(PurityComparer.PurityFields.MsStatus.toString(), b -> b.msStatus = alternateValueSource.Purity.microsatelliteStatus());
        fieldToAlternateValueInitializer.put(PurityComparer.PurityFields.TmbStatus.toString(), b -> b.tmbStatus = alternateValueSource.Purity.tumorMutationalBurdenStatus());
        fieldToAlternateValueInitializer.put(PurityComparer.PurityFields.TmlStatus.toString(), b -> b.tmlStatus = alternateValueSource.Purity.tumorMutationalLoadStatus());
        fieldToAlternateValueInitializer.put(PurityComparer.PurityFields.TincLevel.toString(), b -> b.tincLevel = alternateValueSource.Purity.qc().tincLevel());

        nameToAlternateIndexInitializer = Collections.emptyMap();
        reportabilityFieldToFalseReportabilityInitializer = Collections.emptyMap();
        nameToNonPassInitializer = Collections.emptyMap();
    }
}
