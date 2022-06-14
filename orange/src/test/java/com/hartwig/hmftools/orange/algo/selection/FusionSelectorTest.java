package com.hartwig.hmftools.orange.algo.selection;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.fusion.KnownFusionType;
import com.hartwig.hmftools.common.linx.LinxTestFactory;
import com.hartwig.hmftools.common.protect.ProtectEventGenerator;
import com.hartwig.hmftools.common.protect.ProtectEvidence;
import com.hartwig.hmftools.common.protect.ProtectTestFactory;
import com.hartwig.hmftools.common.sv.linx.LinxFusion;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;
import org.junit.Test;

public class FusionSelectorTest {

    @Test
    public void canSelectFusionsWithEvidence() {
        LinxFusion withEvidence = LinxTestFactory.builder().name("with evidence").geneStart("gene 1").geneEnd("gene 2").build();
        LinxFusion noEvidence = LinxTestFactory.builder().name("no evidence").geneStart("gene 3").geneEnd("gene 4").build();

        ProtectEvidence evidence = ProtectTestFactory.builder().event(ProtectEventGenerator.fusionEvent(withEvidence)).build();

        List<LinxFusion> fusions =
                FusionSelector.selectNonDriverFusions(Lists.newArrayList(withEvidence, noEvidence), Lists.newArrayList(evidence));

        assertEquals(1, fusions.size());
        assertNotNull(findByName(fusions, "with evidence"));
    }

    @Test
    public void canSelectFusionsWithReportedType() {
        LinxFusion withReportedType =
                LinxTestFactory.builder().name("with reported").reportedType(KnownFusionType.KNOWN_PAIR.toString()).build();
        LinxFusion withoutReportedType =
                LinxTestFactory.builder().name("without reported").reportedType(KnownFusionType.NONE.toString()).build();

        List<LinxFusion> fusions =
                FusionSelector.selectNonDriverFusions(Lists.newArrayList(withReportedType, withoutReportedType), Lists.newArrayList());

        assertEquals(1, fusions.size());
        assertNotNull(findByName(fusions, "with reported"));
    }

    @Nullable
    private static LinxFusion findByName(@NotNull List<LinxFusion> fusions, @NotNull String nameToFind) {
        for (LinxFusion fusion : fusions) {
            if (fusion.name().equals(nameToFind)) {
                return fusion;
            }
        }
        return null;
    }
}