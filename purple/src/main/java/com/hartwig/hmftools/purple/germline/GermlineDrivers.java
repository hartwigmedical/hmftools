package com.hartwig.hmftools.purple.germline;

import static java.util.stream.Collectors.toList;
import static java.util.stream.Collectors.toSet;

import static com.hartwig.hmftools.common.driver.LikelihoodMethod.GERMLINE;
import static com.hartwig.hmftools.common.driver.LikelihoodMethod.SPLICE_REGION;
import static com.hartwig.hmftools.common.driver.panel.DriverGeneGermlineReporting.NONE;
import static com.hartwig.hmftools.purple.drivers.SomaticVariantDrivers.getWorstReportableCodingEffect;

import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.driver.DriverCatalog;
import com.hartwig.hmftools.common.driver.DriverImpact;
import com.hartwig.hmftools.common.driver.DriverType;
import com.hartwig.hmftools.common.driver.ImmutableDriverCatalog;
import com.hartwig.hmftools.common.driver.LikelihoodMethod;
import com.hartwig.hmftools.common.driver.panel.DriverGene;
import com.hartwig.hmftools.common.purple.GeneCopyNumber;
import com.hartwig.hmftools.common.purple.ReportedStatus;
import com.hartwig.hmftools.common.variant.CodingEffect;

public class GermlineDrivers
{
    private final Map<String,DriverGene> mDriverGeneMap;

    public GermlineDrivers(final Map<String,DriverGene> driverGenes)
    {
        mDriverGeneMap = driverGenes;
    }

    public List<DriverCatalog> findDrivers(final List<GermlineVariant> variants, final Map<String,GeneCopyNumber> geneCopyNumberMap)
    {
        Set<String> genes = variants.stream().map(GermlineVariant::gene).collect(toSet());

        List<DriverCatalog> driverCatalog = Lists.newArrayList();

        for(String gene : genes)
        {
            DriverGene driverGene = mDriverGeneMap.get(gene);

            if(driverGene == null)
                continue;

            List<GermlineVariant> geneVariants = variants.stream().filter(x -> x.gene().equals(gene)).collect(toList());

            GeneCopyNumber geneCopyNumber = geneCopyNumberMap.get(gene);

            if(geneCopyNumber == null)
                continue;

            driverCatalog.add(germlineDriver(driverGene, gene, geneVariants, geneCopyNumber));
        }

        return driverCatalog;
    }

    private static DriverCatalog germlineDriver(
            final DriverGene driverGene, final String gene, final List<GermlineVariant> geneVariants, final GeneCopyNumber geneCopyNumber)
    {
        Map<DriverImpact,Integer> variantCounts = Maps.newHashMap();

        boolean hasCodingImpact = false;

        for(GermlineVariant variant : geneVariants)
        {
            CodingEffect codingEffect = getWorstReportableCodingEffect(variant.variantImpact());
            DriverImpact driverImpact = DriverImpact.select(variant.type(), codingEffect);

            hasCodingImpact |= variant.isHotspot() || driverImpact != DriverImpact.UNKNOWN;

            Integer count = variantCounts.get(driverImpact);
            variantCounts.put(driverImpact, count != null ? count + 1 : 1);
        }

        int missenseVariants = variantCounts.getOrDefault(DriverImpact.MISSENSE, 0);
        int nonsenseVariants = variantCounts.getOrDefault(DriverImpact.NONSENSE, 0);
        int spliceVariants = variantCounts.getOrDefault(DriverImpact.SPLICE, 0);
        int inframeVariants = variantCounts.getOrDefault(DriverImpact.INFRAME, 0);
        int frameshiftVariants = variantCounts.getOrDefault(DriverImpact.FRAMESHIFT, 0);

        ReportedStatus reportedStatus = ReportedStatus.NOT_REPORTED;

        if(driverGene.reportGermlineVariant() != NONE || driverGene.reportGermlineHotspot() != NONE)
        {
            for(GermlineVariant variant : geneVariants)
            {
                if(variant.reported())
                    reportedStatus = ReportedStatus.REPORTED;
            }
        }

        double impactLikelihood = hasCodingImpact ? 1 : 0;
        LikelihoodMethod likelihoodMethod = hasCodingImpact ? GERMLINE : SPLICE_REGION;

        final ImmutableDriverCatalog.Builder builder = ImmutableDriverCatalog.builder()
                .chromosome(geneVariants.get(0).chromosome())
                .chromosomeBand(geneCopyNumber.ChromosomeBand)
                .gene(gene)
                .transcript(geneCopyNumber.TransName)
                .isCanonical(geneCopyNumber.IsCanonical)
                .driver(DriverType.GERMLINE_MUTATION)
                .category(driverGene.likelihoodType())
                .driverLikelihood(impactLikelihood)
                .missense(missenseVariants)
                .nonsense(nonsenseVariants)
                .splice(spliceVariants)
                .inframe(inframeVariants)
                .frameshift(frameshiftVariants)
                .biallelic(geneVariants.stream().anyMatch(GermlineVariant::biallelic))
                .minCopyNumber(geneCopyNumber.minCopyNumber())
                .maxCopyNumber(geneCopyNumber.maxCopyNumber())
                .likelihoodMethod(likelihoodMethod)
                .reportedStatus(reportedStatus);

        return builder.build();
    }

    @VisibleForTesting
    public Map<String,DriverGene> driverGeneMap() { return mDriverGeneMap; }
}
