package com.hartwig.hmftools.compar.vchord;

import static com.hartwig.hmftools.compar.vchord.VChordData.FLD_BREAST;
import static com.hartwig.hmftools.compar.vchord.VChordData.FLD_OTHER;
import static com.hartwig.hmftools.compar.vchord.VChordData.FLD_OVARIAN;
import static com.hartwig.hmftools.compar.vchord.VChordData.FLD_PANCREATIC;
import static com.hartwig.hmftools.compar.vchord.VChordData.FLD_PROSTATE;

import java.util.Collections;
import java.util.Map;

import com.hartwig.hmftools.compar.ComparConfig;
import com.hartwig.hmftools.compar.ComparableItemTest;

import org.junit.Before;

public class VChordDataTest extends ComparableItemTest<VChordData, VChordComparer, TestVChordDataBuilder>
{
    @Before
    public void setUp()
    {
        comparer = new VChordComparer(new ComparConfig());
        builder = TestVChordDataBuilder.BUILDER;
        VChordData alternateValueSource = builder.createWithAlternateDefaults();
        fieldToAlternateValueInitializer = Map.of(
                FLD_BREAST, b -> b.breast = alternateValueSource.VChord().breastCancerHrdScore(),
                FLD_OVARIAN, b -> b.ovarian = alternateValueSource.VChord().ovarianCancerHrdScore(),
                FLD_PANCREATIC, b -> b.pancreatic = alternateValueSource.VChord().pancreaticCancerScore(),
                FLD_PROSTATE, b -> b.prostate = alternateValueSource.VChord().prostateCancerScore(),
                FLD_OTHER, b -> b.other = alternateValueSource.VChord().otherCancerScore()
        );
        nameToAlternateIndexInitializer = Collections.emptyMap();
        reportabilityFieldToFalseReportabilityInitializer = Collections.emptyMap();
        nameToNonPassInitializer = Collections.emptyMap();
    }
}
