package com.hartwig.hmftools.finding.util;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.hartwig.hmftools.finding.FindingKeys;
import com.hartwig.hmftools.finding.FindingUtil;
import com.hartwig.hmftools.finding.datamodel.FindingRecord;
import com.hartwig.hmftools.finding.datamodel.FindingRecordBuilder;
import com.hartwig.hmftools.finding.datamodel.SmallVariant;
import com.hartwig.hmftools.finding.datamodel.SmallVariantBuilder;
import com.hartwig.hmftools.finding.datamodel.driver.DriverFields;
import com.hartwig.hmftools.finding.datamodel.driver.DriverFieldsBuilder;
import com.hartwig.hmftools.finding.datamodel.driver.DriverFindingList;
import com.hartwig.hmftools.finding.datamodel.driver.DriverFindingListBuilder;
import com.hartwig.hmftools.finding.datamodel.driver.DriverSource;

import jakarta.validation.constraints.NotNull;

public class NoGermlineConverter
{
    public static FindingRecord convert(FindingRecord record)
    {
        return FindingRecordBuilder.builder(record)
                .somaticSmallVariants(combineVariants(record))
                .germlineSmallVariants(FindingUtil.notAvailableDriverFindingList(Set.of()))
                .build();
    }

    @NotNull
    private static DriverFindingList<SmallVariant> combineVariants(@NotNull FindingRecord record)
    {
        DriverFindingList<SmallVariant> somaticFindings = record.somaticSmallVariants();
        DriverFindingList<SmallVariant> germlineFindings = record.germlineSmallVariants();
        List<SmallVariant> allVariants = new ArrayList<>(record.somaticSmallVariants().findings());
        if(germlineFindings.status().isOK())
        {
            allVariants.addAll(germlineFindings.findings().stream().map(v -> SmallVariantBuilder.builder(v)
                    .driver(DriverFieldsBuilder.builder(v.driver())
                            .findingKey(FindingKeys.smallVariant(DriverSource.SOMATIC, v))
                            .driverSource(DriverSource.SOMATIC)
                            .build())
                    .build()).toList());
        }
        return DriverFindingListBuilder.<SmallVariant>builder()
                .status(somaticFindings.status())
                .findings(setMaxDriverLikelihoods(allVariants))
                .build();
    }

    @NotNull
    static List<SmallVariant> setMaxDriverLikelihoods(@NotNull List<SmallVariant> variants)
    {
        Map<String, DriverFields> maximumDriverLikelihoodPerGene = determineMaximumDriverLikelihoodPerGene(variants);
        return variants.stream().map(v ->
        {
            DriverFields driverFields = maximumDriverLikelihoodPerGene.get(v.gene());
            return SmallVariantBuilder.builder(v)
                    .driver(DriverFieldsBuilder.builder(driverFields).build())
                    .build();
        }).toList();
    }

    @NotNull
    static Map<String, DriverFields> determineMaximumDriverLikelihoodPerGene(
            @NotNull Collection<SmallVariant> variants)
    {
        return variants.stream()
                .collect(Collectors.groupingBy(SmallVariant::gene,
                        Collectors.collectingAndThen(Collectors.minBy(Comparator.comparing(SmallVariant::driverLikelihood)), optional -> optional.map(SmallVariant::driver)
                                .orElse(null))));
    }
}
