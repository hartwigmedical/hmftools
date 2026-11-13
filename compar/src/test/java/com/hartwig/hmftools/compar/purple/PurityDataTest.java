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
        fieldToAlternateValueInitializer.put(PurityComparer.Fields.Purity.toString(), b -> b.purity = alternateValueSource.Purity.bestFit().purity());
        fieldToAlternateValueInitializer.put(PurityComparer.Fields.Ploidy.toString(), b -> b.ploidy = alternateValueSource.Purity.bestFit().ploidy());
        fieldToAlternateValueInitializer.put(PurityComparer.Fields.Contamination.toString(), b -> b.contamination = alternateValueSource.Purity.qc().contamination());
        fieldToAlternateValueInitializer.put(PurityComparer.Fields.TmbPerMb.toString(), b -> b.tmb = alternateValueSource.Purity.tumorMutationalBurdenPerMb());
        fieldToAlternateValueInitializer.put(PurityComparer.Fields.MsIndelsPerMb.toString(), b -> b.msIndels = alternateValueSource.Purity.microsatelliteIndelsPerMb());
        fieldToAlternateValueInitializer.put(PurityComparer.Fields.Tml.toString(), b -> b.tml = alternateValueSource.Purity.tumorMutationalLoad());
        fieldToAlternateValueInitializer.put(PurityComparer.Fields.CopyNumberSegments.toString(), b -> b.copyNumberSegments = alternateValueSource.Purity.qc().copyNumberSegments());
        fieldToAlternateValueInitializer.put(PurityComparer.Fields.UnsupportedCopyNumberSegments.toString(), b -> b.unsupportedCopyNumberSegments =
                alternateValueSource.Purity.qc().unsupportedCopyNumberSegments());
        fieldToAlternateValueInitializer.put(PurityComparer.Fields.SvTmb.toString(), b -> b.svTmb = alternateValueSource.Purity.svTumorMutationalBurden());
        fieldToAlternateValueInitializer.put(PurityComparer.Fields.QcStatus.toString(), b -> b.qcStatus = alternateValueSource.Purity.qc().status());
        fieldToAlternateValueInitializer.put(PurityComparer.Fields.Gender.toString(), b -> b.gender = alternateValueSource.Purity.gender());
        fieldToAlternateValueInitializer.put(PurityComparer.Fields.GermlineAberrations.toString(), b -> b.germlineAberrations = alternateValueSource.Purity.qc()
                .germlineAberrations());
        fieldToAlternateValueInitializer.put(PurityComparer.Fields.FitMethod.toString(), b -> b.fitMethod = alternateValueSource.Purity.method());
        fieldToAlternateValueInitializer.put(PurityComparer.Fields.MsStatus.toString(), b -> b.msStatus = alternateValueSource.Purity.microsatelliteStatus());
        fieldToAlternateValueInitializer.put(PurityComparer.Fields.TmbStatus.toString(), b -> b.tmbStatus = alternateValueSource.Purity.tumorMutationalBurdenStatus());
        fieldToAlternateValueInitializer.put(PurityComparer.Fields.TmlStatus.toString(), b -> b.tmlStatus = alternateValueSource.Purity.tumorMutationalLoadStatus());
        fieldToAlternateValueInitializer.put(PurityComparer.Fields.TincLevel.toString(), b -> b.tincLevel = alternateValueSource.Purity.qc().tincLevel());

        nameToAlternateIndexInitializer = Collections.emptyMap();
        reportabilityFieldToFalseReportabilityInitializer = Collections.emptyMap();
        nameToNonPassInitializer = Collections.emptyMap();
    }
}
