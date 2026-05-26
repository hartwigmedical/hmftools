package com.hartwig.hmftools.finding.util;

import static com.hartwig.hmftools.finding.util.TransformUtil.transformDriverFindingList;

import com.hartwig.hmftools.finding.datamodel.Disruption;
import com.hartwig.hmftools.finding.datamodel.DisruptionBuilder;
import com.hartwig.hmftools.finding.datamodel.FindingRecord;
import com.hartwig.hmftools.finding.datamodel.FindingRecordBuilder;
import com.hartwig.hmftools.finding.datamodel.GainDeletion;
import com.hartwig.hmftools.finding.datamodel.GainDeletionBuilder;
import com.hartwig.hmftools.finding.datamodel.SmallVariant;
import com.hartwig.hmftools.finding.datamodel.SmallVariantBuilder;

import jakarta.validation.constraints.NotNull;

public class ChromosomePrefixTransformer
{
    @NotNull
    public static FindingRecord transform(@NotNull FindingRecord record)
    {
        return FindingRecordBuilder.builder(record)
                .somaticSmallVariants(transformDriverFindingList(record.somaticSmallVariants(), ChromosomePrefixTransformer::transform))
                .germlineSmallVariants(transformDriverFindingList(record.germlineSmallVariants(), ChromosomePrefixTransformer::transform))
                .somaticDisruptions(transformDriverFindingList(record.somaticDisruptions(), ChromosomePrefixTransformer::transform))
                .germlineDisruptions(transformDriverFindingList(record.germlineDisruptions(), ChromosomePrefixTransformer::transform))
                .somaticGainDeletions(transformDriverFindingList(record.somaticGainDeletions(), ChromosomePrefixTransformer::transform))
                .germlineGainDeletions(transformDriverFindingList(record.germlineGainDeletions(), ChromosomePrefixTransformer::transform))
                .build();
    }

    private static SmallVariant transform(SmallVariant smallVariant)
    {
        return SmallVariantBuilder.builder(smallVariant)
                .chromosome(chromosome(smallVariant.chromosome()))
                .build();
    }

    private static Disruption transform(Disruption smallVariant)
    {
        return DisruptionBuilder.builder(smallVariant)
                .chromosome(chromosome(smallVariant.chromosome()))
                .build();
    }

    private static GainDeletion transform(GainDeletion gainDeletion)
    {
        return GainDeletionBuilder.builder(gainDeletion)
                .chromosome(chromosome(gainDeletion.chromosome()))
                .build();
    }

    @NotNull
    private static String chromosome(@NotNull String chromosome) {
        return chromosome.replace("chr", "");
    }
}
