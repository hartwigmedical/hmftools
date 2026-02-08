package com.hartwig.hmftools.orange.algo.isofox;

import static com.hartwig.hmftools.common.rna.AltSpliceJunctionType.SKIPPED_EXONS;
import static com.hartwig.hmftools.orange.algo.isofox.FusionNameUtil.geneDown;
import static com.hartwig.hmftools.orange.algo.isofox.FusionNameUtil.geneUp;

import java.util.List;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.driver.panel.DriverGene;
import com.hartwig.hmftools.common.fusion.KnownFusionCache;
import com.hartwig.hmftools.common.rna.AltSpliceJunctionType;
import com.hartwig.hmftools.common.rna.KnownFusionType;
import com.hartwig.hmftools.common.rna.NovelSpliceJunction;
import com.hartwig.hmftools.common.rna.RnaFusion;
import com.hartwig.hmftools.datamodel.linx.LinxFusion;
import com.hartwig.hmftools.orange.algo.linx.DnaFusionEvaluator;

import org.jetbrains.annotations.NotNull;

final class NovelSpliceJunctionSelector
{
    public static List<NovelSpliceJunction> selectSkippedExons(
            final List<RnaFusion> rnaFusions, final List<NovelSpliceJunction> junctions, final List<LinxFusion> linxFusions)
    {
        List<NovelSpliceJunction> result = Lists.newArrayList();

        for(RnaFusion rnaFusion : rnaFusions)
        {
            String geneUp = geneUp(rnaFusion);
            String geneDown = geneDown(rnaFusion);
            if(geneUp != null && geneDown != null && geneUp.equals(geneDown)
            && rnaFusion.knownType() == KnownFusionType.KNOWN_PAIR)
            {
                for(NovelSpliceJunction junction : junctions)
                {
                    if(!junction.geneName().equals(geneUp))
                        continue;
                    if(junction.type() != SKIPPED_EXONS)
                        continue;

                    boolean isTypeMatch = junction.type() == SKIPPED_EXONS;
                    boolean hasSufficientFragments = junction.fragmentCount() > 5;
                    boolean hasLimitedCohortFreq = junction.cohortFrequency() < 30;
                    boolean hasReportedLinxFusion = DnaFusionEvaluator.hasFusion(linxFusions, junction.geneName(), junction.geneName());
                    if(isTypeMatch && hasSufficientFragments && hasLimitedCohortFreq && !hasReportedLinxFusion)
                    {
                        result.add(junction);
                    }
                }
            }
        }

        return result;
    }

    public static List<NovelSpliceJunction> selectNovelExonsIntrons(
            final List<NovelSpliceJunction> junctions, final List<DriverGene> driverGenes)
    {
        List<NovelSpliceJunction> result = Lists.newArrayList();

        Set<String> drivers = Sets.newHashSet();
        for(DriverGene driverGene : driverGenes)
        {
            drivers.add(driverGene.gene());
        }

        for(NovelSpliceJunction junction : junctions)
        {
            if(drivers.contains(junction.geneName()))
            {
                boolean isTypeMatch =
                        junction.type() == AltSpliceJunctionType.NOVEL_INTRON || junction.type() == AltSpliceJunctionType.NOVEL_EXON;
                boolean hasSufficientFragments = junction.fragmentCount() > 5;
                boolean hasLimitedCohortFreq = junction.cohortFrequency() < 10;
                if(isTypeMatch && hasSufficientFragments && hasLimitedCohortFreq)
                {
                    result.add(junction);
                }
            }
        }

        return result;
    }
}
