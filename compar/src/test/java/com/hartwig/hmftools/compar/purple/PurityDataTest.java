package com.hartwig.hmftools.compar.purple;

import static com.hartwig.hmftools.compar.purple.PurityData.FLD_CN_SEGS;
import static com.hartwig.hmftools.compar.purple.PurityData.FLD_CONTAMINATION;
import static com.hartwig.hmftools.compar.purple.PurityData.FLD_FIT_METHOD;
import static com.hartwig.hmftools.compar.purple.PurityData.FLD_GENDER;
import static com.hartwig.hmftools.compar.purple.PurityData.FLD_GERM_ABS;
import static com.hartwig.hmftools.compar.purple.PurityData.FLD_MS_INDELS;
import static com.hartwig.hmftools.compar.purple.PurityData.FLD_MS_STATUS;
import static com.hartwig.hmftools.compar.purple.PurityData.FLD_PLOIDY;
import static com.hartwig.hmftools.compar.purple.PurityData.FLD_PURITY;
import static com.hartwig.hmftools.compar.purple.PurityData.FLD_QC_STATUS;
import static com.hartwig.hmftools.compar.purple.PurityData.FLD_SV_TMB;
import static com.hartwig.hmftools.compar.purple.PurityData.FLD_TINC_LEVEL;
import static com.hartwig.hmftools.compar.purple.PurityData.FLD_TMB;
import static com.hartwig.hmftools.compar.purple.PurityData.FLD_TMB_STATUS;
import static com.hartwig.hmftools.compar.purple.PurityData.FLD_TML;
import static com.hartwig.hmftools.compar.purple.PurityData.FLD_TML_STATUS;
import static com.hartwig.hmftools.compar.purple.PurityData.FLD_UNS_CN_SEGS;

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
        comparer = new PurityComparer(new ComparConfig());
        builder = TestPurityDataBuilder.BUILDER;
        PurityData alternateValueSource = builder.createWithAlternateDefaults();

        fieldToAlternateValueInitializer = new HashMap<>();
        fieldToAlternateValueInitializer.put(FLD_PURITY, b -> b.purity = alternateValueSource.Purity.bestFit().purity());
        fieldToAlternateValueInitializer.put(FLD_PLOIDY, b -> b.ploidy = alternateValueSource.Purity.bestFit().ploidy());
        fieldToAlternateValueInitializer.put(FLD_CONTAMINATION, b -> b.contamination = alternateValueSource.Purity.qc().contamination());
        fieldToAlternateValueInitializer.put(FLD_TMB, b -> b.tmb = alternateValueSource.Purity.tumorMutationalBurdenPerMb());
        fieldToAlternateValueInitializer.put(FLD_MS_INDELS, b -> b.msIndels = alternateValueSource.Purity.microsatelliteIndelsPerMb());
        fieldToAlternateValueInitializer.put(FLD_TML, b -> b.tml = alternateValueSource.Purity.tumorMutationalLoad());
        fieldToAlternateValueInitializer.put(FLD_CN_SEGS, b -> b.copyNumberSegments = alternateValueSource.Purity.qc().copyNumberSegments());
        fieldToAlternateValueInitializer.put(FLD_UNS_CN_SEGS, b -> b.unsupportedCopyNumberSegments =
                alternateValueSource.Purity.qc().unsupportedCopyNumberSegments());
        fieldToAlternateValueInitializer.put(FLD_SV_TMB, b -> b.svTmb = alternateValueSource.Purity.svTumorMutationalBurden());
        fieldToAlternateValueInitializer.put(FLD_QC_STATUS, b -> b.qcStatus = alternateValueSource.Purity.qc().status());
        fieldToAlternateValueInitializer.put(FLD_GENDER, b -> b.gender = alternateValueSource.Purity.gender());
        fieldToAlternateValueInitializer.put(FLD_GERM_ABS, b -> b.germlineAberrations = alternateValueSource.Purity.qc()
                .germlineAberrations());
        fieldToAlternateValueInitializer.put(FLD_FIT_METHOD, b -> b.fitMethod = alternateValueSource.Purity.method());
        fieldToAlternateValueInitializer.put(FLD_MS_STATUS, b -> b.msStatus = alternateValueSource.Purity.microsatelliteStatus());
        fieldToAlternateValueInitializer.put(FLD_TMB_STATUS, b -> b.tmbStatus = alternateValueSource.Purity.tumorMutationalBurdenStatus());
        fieldToAlternateValueInitializer.put(FLD_TML_STATUS, b -> b.tmlStatus = alternateValueSource.Purity.tumorMutationalLoadStatus());
        fieldToAlternateValueInitializer.put(FLD_TINC_LEVEL, b -> b.tincLevel = alternateValueSource.Purity.qc().tincLevel());

        nameToAlternateIndexInitializer = Collections.emptyMap();
        reportabilityFieldToFalseReportabilityInitializer = Collections.emptyMap();
    }
}
