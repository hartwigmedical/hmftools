package com.hartwig.hmftools.orange.algo.isofox;

import static org.junit.Assert.assertEquals;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.fusion.KnownFusionCache;
import com.hartwig.hmftools.common.fusion.KnownFusionCacheTestFactory;
import com.hartwig.hmftools.common.isofox.IsofoxTestFactory;
import com.hartwig.hmftools.common.rna.RnaFusion;
import com.hartwig.hmftools.common.sv.StructuralVariantType;
import com.hartwig.hmftools.datamodel.finding.Fusion;
import com.hartwig.hmftools.orange.algo.linx.LinxOrangeTestFactory;

import org.junit.Test;

public class RnaFusionSelectorTest
{
    @Test
    public void canSelectNovelKnownFusions()
    {
        RnaFusion match = IsofoxTestFactory.rnaFusionBuilder().name("A_B").build();
        RnaFusion hasLinxFusionAlready = IsofoxTestFactory.rnaFusionBuilder().name("C_D").build();
        RnaFusion noKnownFusion = IsofoxTestFactory.rnaFusionBuilder().name("E_F").build();
        List<RnaFusion> rnaFusions = Lists.newArrayList(match, hasLinxFusionAlready, noKnownFusion);

        KnownFusionCache knownFusionCache = new KnownFusionCache();
        knownFusionCache.addData(KnownFusionCacheTestFactory.createKnownPair("A", "B"));
        knownFusionCache.addData(KnownFusionCacheTestFactory.createKnownPair("C", "D"));

        List<Fusion> linxFusions = Lists.newArrayList(LinxOrangeTestFactory.fusionBuilder().geneStart("C").geneEnd("D").build());

        List<RnaFusion> novelFusions = RnaFusionSelector.selectNovelKnownFusions(rnaFusions, linxFusions, knownFusionCache);
        assertEquals(1, novelFusions.size());
        assertEquals(match, novelFusions.get(0));
    }

    @Test
    public void canSelectNovelPromiscuousFusions()
    {
        RnaFusion match =
                IsofoxTestFactory.rnaFusionBuilder().name("A_B").svType(StructuralVariantType.BND).positionUp(1).positionDown(2).build();

        RnaFusion tooShortDel =
                IsofoxTestFactory.rnaFusionBuilder().name("C_D").svType(StructuralVariantType.DEL).positionUp(1).positionDown(2).build();

        RnaFusion knownFusion = IsofoxTestFactory.rnaFusionBuilder().name("E_F").svType(StructuralVariantType.BND).build();

        RnaFusion noPromiscuous = IsofoxTestFactory.rnaFusionBuilder().name("G_H").svType(StructuralVariantType.BND).build();
        List<RnaFusion> rnaFusions = Lists.newArrayList(match, tooShortDel, knownFusion, noPromiscuous);

        KnownFusionCache knownFusionCache = new KnownFusionCache();
        knownFusionCache.addData(KnownFusionCacheTestFactory.createPromiscuousFive("A"));
        knownFusionCache.addData(KnownFusionCacheTestFactory.createPromiscuousThree("G"));
        knownFusionCache.addData(KnownFusionCacheTestFactory.createKnownPair("E", "F"));

        List<RnaFusion> novelPromiscuous =
                RnaFusionSelector.selectNovelPromiscuousFusions(rnaFusions, Lists.newArrayList(), knownFusionCache);

        assertEquals(1, novelPromiscuous.size());
        assertEquals(match, novelPromiscuous.get(0));
    }
}