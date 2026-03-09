package com.hartwig.hmftools.orange.conversion;

import static com.hartwig.hmftools.datamodel.isofox.RnaFusionType.EXON_DEL_DUP;
import static com.hartwig.hmftools.datamodel.isofox.RnaFusionType.KNOWN_PAIR;
import static com.hartwig.hmftools.datamodel.isofox.RnaFusionType.NONE;
import static com.hartwig.hmftools.datamodel.isofox.RnaFusionType.PROMISCUOUS_3;
import static com.hartwig.hmftools.datamodel.isofox.RnaFusionType.PROMISCUOUS_5;
import static com.hartwig.hmftools.datamodel.isofox.RnaFusionType.PROMISCUOUS_BOTH;
import static com.hartwig.hmftools.orange.algo.isofox.RnaFusionSelector.geneDown;
import static com.hartwig.hmftools.orange.algo.isofox.RnaFusionSelector.geneUp;

import com.hartwig.hmftools.common.rna.NovelSpliceJunction;
import com.hartwig.hmftools.common.rna.RnaFusion;
import com.hartwig.hmftools.common.rna.RnaFusionFile;
import com.hartwig.hmftools.common.rna.RnaQcFilter;
import com.hartwig.hmftools.datamodel.isofox.AltSpliceJunctionContext;
import com.hartwig.hmftools.datamodel.isofox.AltSpliceJunctionType;
import com.hartwig.hmftools.datamodel.isofox.GeneExpression;
import com.hartwig.hmftools.datamodel.isofox.ImmutableGeneExpression;
import com.hartwig.hmftools.datamodel.isofox.ImmutableNovelSpliceJunction;
import com.hartwig.hmftools.datamodel.isofox.ImmutableRnaFusion;
import com.hartwig.hmftools.datamodel.isofox.ImmutableRnaStatistics;
import com.hartwig.hmftools.datamodel.isofox.RnaStatistics;
import com.hartwig.hmftools.datamodel.isofox.RnaQCStatus;
import com.hartwig.hmftools.datamodel.isofox.StructuralVariantType;

public final class IsofoxConversion
{
    public static RnaStatistics convert(final com.hartwig.hmftools.common.rna.RnaStatistics rnaStatistics)
    {
        return ImmutableRnaStatistics.builder()
                .qcStatus(ConversionUtil.mapToIterable(rnaStatistics.qcStatus(), IsofoxConversion::convert))
                .totalFragments(rnaStatistics.totalFragments())
                .duplicateFragments(rnaStatistics.duplicateFragments())
                .splicedFragmentPerc(rnaStatistics.splicedFragmentPerc())
                .unsplicedFragmentPerc(rnaStatistics.unsplicedFragmentPerc())
                .altFragmentPerc(rnaStatistics.altFragmentPerc())
                .chimericFragmentPerc(rnaStatistics.chimericFragmentPerc())
                .build();
    }

    public static RnaQCStatus convert(final RnaQcFilter filter)
    {
        return RnaQCStatus.valueOf(filter.name());
    }

    public static GeneExpression convert(final com.hartwig.hmftools.common.rna.GeneExpression geneExpression)
    {
        boolean hasCancerCohortData = geneExpression.medianTpmCancer() > 0;
        return ImmutableGeneExpression.builder()
                .gene(geneExpression.geneName())
                .tpm(geneExpression.tpm())
                .medianTpmCancer(hasCancerCohortData ? geneExpression.medianTpmCancer() : null)
                .percentileCancer(hasCancerCohortData ? geneExpression.percentileCancer() : null)
                .medianTpmCohort(geneExpression.medianTpmCohort())
                .percentileCohort(geneExpression.percentileCohort())
                .build();
    }

    public static com.hartwig.hmftools.datamodel.isofox.RnaFusion convert(final RnaFusion rnaFusion)
    {
        return ImmutableRnaFusion.builder()
                .geneStart(geneUp(rnaFusion))
                .geneEnd(geneDown(rnaFusion))
                .chromosomeStart(rnaFusion.chromosomeUp())
                .chromosomeEnd(rnaFusion.chromosomeDown())
                .positionStart(rnaFusion.positionUp())
                .positionEnd(rnaFusion.positionDown())
                .junctionTypeStart(rnaFusion.junctionTypeUp())
                .junctionTypeEnd(rnaFusion.junctionTypeDown())
                .svType(StructuralVariantType.valueOf(rnaFusion.svType().name()))
                .knownType(convertKnownType(rnaFusion))
                .splitFragments(rnaFusion.splitFragments())
                .realignedFrags(rnaFusion.realignedFrags())
                .discordantFrags(rnaFusion.discordantFrags())
                .depthStart(rnaFusion.depthUp())
                .depthEnd(rnaFusion.depthDown())
                .cohortFrequency(rnaFusion.cohortFrequency())
                .build();
    }

    private static com.hartwig.hmftools.datamodel.isofox.RnaFusionType convertKnownType(final RnaFusion rnaFusion)
    {
        String[] geneNames = RnaFusionFile.geneNames(rnaFusion);

        if(geneNames.length == 2 && geneNames[0].equals(geneNames[1]))
            return EXON_DEL_DUP;

        switch(rnaFusion.knownType())
        {
            case KNOWN_PAIR: return KNOWN_PAIR;
            case PROM5_PROM3: return PROMISCUOUS_BOTH;

            case KNOWN_PROM3:
            case OTHER_PROM3:
                return PROMISCUOUS_3;

            case PROM5_KNOWN:
            case PROM5_OTHER:
                return PROMISCUOUS_5;

            default: return NONE;
        }
    }

    public static com.hartwig.hmftools.datamodel.isofox.NovelSpliceJunction convert(final NovelSpliceJunction novelSpliceJunction)
    {
        return ImmutableNovelSpliceJunction.builder()
                .gene(novelSpliceJunction.geneName())
                .chromosome(novelSpliceJunction.chromosome())
                .junctionStart(novelSpliceJunction.junctionStart())
                .junctionEnd(novelSpliceJunction.junctionEnd())
                .type(AltSpliceJunctionType.valueOf(novelSpliceJunction.type().name()))
                .fragmentCount(novelSpliceJunction.fragmentCount())
                .depthStart(novelSpliceJunction.depthStart())
                .depthEnd(novelSpliceJunction.depthEnd())
                .regionStart(AltSpliceJunctionContext.valueOf(novelSpliceJunction.regionStart().name()))
                .regionEnd(AltSpliceJunctionContext.valueOf(novelSpliceJunction.regionEnd().name()))
                .cohortFrequency(novelSpliceJunction.cohortFrequency())
                .build();
    }
}
