package com.hartwig.hmftools.finding.util;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.function.Function;
import java.util.stream.Collectors;

import com.hartwig.hmftools.finding.datamodel.Disruption;
import com.hartwig.hmftools.finding.datamodel.DisruptionBuilder;
import com.hartwig.hmftools.finding.datamodel.FindingRecord;
import com.hartwig.hmftools.finding.datamodel.FindingRecordBuilder;
import com.hartwig.hmftools.finding.datamodel.GainDeletion;
import com.hartwig.hmftools.finding.datamodel.GainDeletionBuilder;
import com.hartwig.hmftools.finding.datamodel.SmallVariant;
import com.hartwig.hmftools.finding.datamodel.SmallVariantBuilder;
import com.hartwig.hmftools.finding.datamodel.driver.Driver;
import com.hartwig.hmftools.finding.datamodel.driver.DriverFields;
import com.hartwig.hmftools.finding.datamodel.driver.DriverFieldsBuilder;
import com.hartwig.hmftools.finding.datamodel.driver.DriverFindingList;
import com.hartwig.hmftools.finding.datamodel.driver.DriverFindingListBuilder;
import com.hartwig.hmftools.finding.datamodel.driver.DriverSource;

import jakarta.validation.constraints.NotNull;

// TODO: Finding keys are still recognizable as germline
//  however if we would convert them we have to make sure these remain unique
public class NoGermlineTransformer
{
    public static FindingRecord transform(FindingRecord record)
    {
        return FindingRecordBuilder.builder(record)
                .somaticSmallVariants(combineVariants(record))
                .germlineSmallVariants(FindingUtil.notAvailableDriverFindingList(Set.of()))
                .somaticDisruptions(combineDisruptions(record))
                .germlineDisruptions(FindingUtil.notAvailableDriverFindingList(Set.of()))
                .somaticGainDeletions(combineGainDeletions(record))
                .germlineGainDeletions(FindingUtil.notAvailableDriverFindingList(Set.of()))
                .build();
    }

    @NotNull
    private static DriverFindingList<SmallVariant> combineVariants(@NotNull FindingRecord record)
    {
        DriverFindingList<SmallVariant> somaticFindings = record.somaticSmallVariants();
        DriverFindingList<SmallVariant> germlineFindings = record.germlineSmallVariants();
        return DriverFindingListBuilder.builder(somaticFindings)
                .findings(setMaxDriverLikelihoods(combineFindings(somaticFindings, germlineFindings, NoGermlineTransformer::transformSmallVariant)))
                .build();
    }

    private static SmallVariant transformSmallVariant(SmallVariant smallVariant)
    {
        // Germline only variant, return as is
        if(smallVariant.driverSource() == DriverSource.GERMLINE && smallVariant.variantCopyNumber() < 0.5)
        {
            return smallVariant;
        }
        else
        {
            return SmallVariantBuilder.builder(smallVariant)
                    .driver(DriverFieldsBuilder.builder(smallVariant.driver())
                            .driverSource(DriverSource.SOMATIC)
                            .build())
                    .build();
        }
    }

    @NotNull
    static List<SmallVariant> setMaxDriverLikelihoods(@NotNull List<SmallVariant> variants)
    {
        Map<String, DriverFields> maximumDriverLikelihoodPerGene = determineMaximumDriverLikelihoodPerGene(variants);
        return variants.stream().map(v ->
        {
            DriverFields driverFields = maximumDriverLikelihoodPerGene.get(v.gene());
            if (v.driverSource() != DriverSource.GERMLINE)
            {
                return SmallVariantBuilder.builder(v)
                        .driver(DriverFieldsBuilder.builder(driverFields).build())
                        .build();
            } else {
                return v;
            }
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

    private static DriverFindingList<Disruption> combineDisruptions(@NotNull FindingRecord record)
    {
        DriverFindingList<Disruption> somaticFindings = record.somaticDisruptions();
        DriverFindingList<Disruption> germlineFindings = record.germlineDisruptions();
        return combineDriverFindingLists(somaticFindings, germlineFindings, NoGermlineTransformer::transformDisruption);
    }

    private static Disruption transformDisruption(Disruption disruption)
    {
        return DisruptionBuilder.builder(disruption)
                .driver(DriverFieldsBuilder.builder(disruption.driver())
                        .driverSource(DriverSource.SOMATIC)
                        .build())
                .build();
    }

    private static DriverFindingList<GainDeletion> combineGainDeletions(@NotNull FindingRecord record)
    {
        //TODO: Review behavior. This was the old behavior (ConclusionAlgo)
        // purple.reportableSomaticGainsLosses();
        // purple.reportableGermlineFullLosses();
        // purple.reportableGermlineLossOfHeterozygosities();
        DriverFindingList<GainDeletion> somaticFindings = record.somaticGainDeletions();
        DriverFindingList<GainDeletion> germlineFindings = record.germlineGainDeletions();
        return combineDriverFindingLists(somaticFindings, germlineFindings, NoGermlineTransformer::transformGainDeletion);
    }

    private static GainDeletion transformGainDeletion(GainDeletion gainDeletion)
    {
        return GainDeletionBuilder.builder(gainDeletion)
                .driver(DriverFieldsBuilder.builder(gainDeletion.driver())
                        .driverSource(DriverSource.SOMATIC)
                        .build())
                .build();
    }

    @NotNull
    private static <T extends Driver> DriverFindingList<T> combineDriverFindingLists(DriverFindingList<T> somaticFindings,
            DriverFindingList<T> germlineFindings,
            Function<T, T> buildFinding)
    {
        return DriverFindingListBuilder.builder(somaticFindings)
                .findings(combineFindings(somaticFindings, germlineFindings, buildFinding))
                .build();
    }

    @NotNull
    private static <T extends Driver> List<T> combineFindings(DriverFindingList<T> somaticFindings,
            DriverFindingList<T> germlineFindings, Function<T, T> buildFinding)
    {
        List<T> combinedFindings = new ArrayList<>(somaticFindings.findings());
        // Only include germline findings if the somatic findings are ok, otherwise this is a germline only report.
        if(somaticFindings.status().isOK() && germlineFindings.status().isOK())
        {
            combinedFindings.addAll(germlineFindings.findings().stream().map(buildFinding).toList());
        }
        return combinedFindings;
    }
}
