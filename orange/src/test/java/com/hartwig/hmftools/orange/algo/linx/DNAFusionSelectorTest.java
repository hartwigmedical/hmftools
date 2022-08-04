package com.hartwig.hmftools.orange.algo.linx;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.drivercatalog.DriverCategory;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGene;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGeneTestFactory;
import com.hartwig.hmftools.common.fusion.KnownFusionType;
import com.hartwig.hmftools.common.linx.LinxTestFactory;
import com.hartwig.hmftools.common.protect.EventGenerator;
import com.hartwig.hmftools.common.protect.ProtectEvidence;
import com.hartwig.hmftools.common.protect.ProtectTestFactory;
import com.hartwig.hmftools.common.sv.linx.FusionPhasedType;
import com.hartwig.hmftools.common.sv.linx.LinxFusion;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;
import org.junit.Test;

public class DNAFusionSelectorTest {

    @Test
    public void canSelectFusionsWithEvidence() {
        LinxFusion withEvidence =
                LinxTestFactory.fusionBuilder().reported(false).name("with evidence").geneStart("gene 1").geneEnd("gene 2").build();
        LinxFusion noEvidence =
                LinxTestFactory.fusionBuilder().reported(false).name("no evidence").geneStart("gene 3").geneEnd("gene 4").build();

        ProtectEvidence evidence = ProtectTestFactory.builder().event(EventGenerator.fusionEvent(withEvidence)).build();

        List<LinxFusion> fusions = DNAFusionSelector.selectInterestingUnreportedFusions(Lists.newArrayList(withEvidence, noEvidence),
                Lists.newArrayList(evidence),
                Lists.newArrayList());

        assertEquals(1, fusions.size());
        assertNotNull(findByName(fusions, "with evidence"));
    }

    @Test
    public void canSelectFusionsWithReportedType() {
        LinxFusion withReportedType = LinxTestFactory.fusionBuilder()
                .reported(false)
                .name("with reported")
                .reportedType(KnownFusionType.KNOWN_PAIR.toString())
                .build();
        LinxFusion withoutReportedType = LinxTestFactory.fusionBuilder()
                .reported(false)
                .name("without reported")
                .reportedType(KnownFusionType.NONE.toString())
                .build();

        List<LinxFusion> fusions =
                DNAFusionSelector.selectInterestingUnreportedFusions(Lists.newArrayList(withReportedType, withoutReportedType),
                        Lists.newArrayList(),
                        Lists.newArrayList());

        assertEquals(1, fusions.size());
        assertNotNull(findByName(fusions, "with reported"));
    }

    @Test
    public void canSelectFusionOfOncogene() {
        DriverGene oncogene = DriverGeneTestFactory.builder().gene("onco").likelihoodType(DriverCategory.ONCO).build();
        DriverGene tsg = DriverGeneTestFactory.builder().gene("tsg").likelihoodType(DriverCategory.TSG).build();

        LinxFusion outOfFrameOnco = LinxTestFactory.fusionBuilder()
                .reported(false)
                .phased(FusionPhasedType.OUT_OF_FRAME)
                .name("out of frame five onco")
                .geneStart(oncogene.gene())
                .build();
        LinxFusion noDriver = LinxTestFactory.fusionBuilder()
                .reported(false)
                .phased(FusionPhasedType.INFRAME)
                .name("no driver")
                .geneStart("no driver")
                .build();
        LinxFusion withFiveOncogene = LinxTestFactory.fusionBuilder()
                .reported(false)
                .phased(FusionPhasedType.INFRAME)
                .name("with five onco")
                .geneStart(oncogene.gene())
                .build();
        LinxFusion withThreeOncogene = LinxTestFactory.fusionBuilder()
                .reported(false)
                .phased(FusionPhasedType.INFRAME)
                .name("with three onco")
                .geneEnd(oncogene.gene())
                .build();
        LinxFusion withFiveTSG = LinxTestFactory.fusionBuilder()
                .reported(false)
                .phased(FusionPhasedType.INFRAME)
                .name("with five tsg")
                .geneStart(tsg.gene())
                .build();
        LinxFusion withThreeTSG = LinxTestFactory.fusionBuilder()
                .reported(false)
                .phased(FusionPhasedType.INFRAME)
                .name("with five tsg")
                .geneEnd(tsg.gene())
                .build();

        List<LinxFusion> fusions = DNAFusionSelector.selectInterestingUnreportedFusions(Lists.newArrayList(outOfFrameOnco,
                noDriver,
                withFiveOncogene,
                withThreeOncogene,
                withFiveTSG,
                withThreeTSG), Lists.newArrayList(), Lists.newArrayList(oncogene, tsg));

        assertEquals(2, fusions.size());
        assertNotNull(findByName(fusions, "with five onco"));
        assertNotNull(findByName(fusions, "with three onco"));
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