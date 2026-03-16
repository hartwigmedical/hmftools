package com.hartwig.hmftools.orange.algo.purple;

import static com.hartwig.hmftools.common.driver.DriverType.AMP;
import static com.hartwig.hmftools.common.driver.DriverType.GERMLINE_MUTATION;
import static com.hartwig.hmftools.common.driver.DriverType.MUTATION;
import static com.hartwig.hmftools.common.driver.DriverType.PARTIAL_AMP;
import static com.hartwig.hmftools.common.linx.DriverEventType.DEL;
import static com.hartwig.hmftools.common.linx.DriverEventType.GAIN;
import static com.hartwig.hmftools.common.linx.DriverEventType.GAIN_ARM;
import static com.hartwig.hmftools.common.linx.DriverEventType.GAIN_CHR;
import static com.hartwig.hmftools.common.linx.DriverEventType.LOH;
import static com.hartwig.hmftools.common.linx.DriverEventType.LOH_ARM;
import static com.hartwig.hmftools.common.linx.DriverEventType.LOH_CHR;
import static com.hartwig.hmftools.common.linx.DriverEventType.LOH_SV_CENTRO;
import static com.hartwig.hmftools.common.linx.DriverEventType.LOH_SV_TELO;
import static com.hartwig.hmftools.orange.OrangeApplication.LOGGER;
import static com.hartwig.hmftools.orange.algo.OrangeConstants.PURPLE_AMP_DEL_FULL;
import static com.hartwig.hmftools.orange.algo.OrangeConstants.PURPLE_AMP_DEL_PARTIAL;
import static com.hartwig.hmftools.orange.algo.OrangeConstants.PURPLE_ARM_CN_DIPLOID;
import static com.hartwig.hmftools.orange.algo.linx.LinxInterpreter.findReportableLinxPlot;

import java.util.EnumSet;
import java.util.List;
import java.util.Objects;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.driver.DriverCatalog;
import com.hartwig.hmftools.common.driver.DriverType;
import com.hartwig.hmftools.common.linx.DriverEventType;
import com.hartwig.hmftools.common.linx.LinxDriver;
import com.hartwig.hmftools.common.purple.GeneCopyNumber;
import com.hartwig.hmftools.common.rna.GeneExpression;
import com.hartwig.hmftools.common.variant.SmallVariant;
import com.hartwig.hmftools.datamodel.purple.ImmutablePurpleCharacteristics;
import com.hartwig.hmftools.datamodel.purple.ImmutablePurpleFit;
import com.hartwig.hmftools.datamodel.purple.ImmutablePurpleGainDeletion;
import com.hartwig.hmftools.datamodel.purple.ImmutablePurpleRecord;
import com.hartwig.hmftools.datamodel.purple.PurpleCharacteristics;
import com.hartwig.hmftools.datamodel.purple.PurpleChrArmCopyNumber;
import com.hartwig.hmftools.datamodel.purple.PurpleDriver;
import com.hartwig.hmftools.datamodel.purple.PurpleFit;
import com.hartwig.hmftools.datamodel.purple.PurpleFittedPurityMethod;
import com.hartwig.hmftools.datamodel.purple.PurpleGainDeletion;
import com.hartwig.hmftools.datamodel.purple.PurpleMicrosatelliteStatus;
import com.hartwig.hmftools.datamodel.purple.PurpleRecord;
import com.hartwig.hmftools.datamodel.purple.PurpleTumorMutationalStatus;
import com.hartwig.hmftools.datamodel.purple.PurpleVariant;
import com.hartwig.hmftools.orange.algo.isofox.IsofoxData;
import com.hartwig.hmftools.orange.algo.linx.LinxData;
import com.hartwig.hmftools.orange.conversion.ConversionUtil;
import com.hartwig.hmftools.orange.conversion.PurpleConversion;

import org.jetbrains.annotations.Nullable;

public class PurpleInterpreter
{
    public PurpleInterpreter() {}

    public PurpleRecord interpret(final PurpleData purpleData, final LinxData linxData, @Nullable final IsofoxData isofoxData)
    {
        List<PurpleVariant> somaticVariants = buildPurpleVariants(purpleData.somaticVariants(), purpleData.somaticDrivers(), false);

        List<PurpleVariant> germlineVariants = buildPurpleVariants(purpleData.germlineVariants(), purpleData.germlineDrivers(), true);

        List<PurpleDriver> germlineDrivers = ConversionUtil.mapToNullableList(purpleData.germlineDrivers(), PurpleConversion::convert);

        List<PurpleGainDeletion> driverSomaticGainsDels = somaticGainsDelsFromDrivers(
                purpleData.somaticDrivers(), purpleData.somaticGeneCopyNumbers(), linxData, isofoxData);

        List<PurpleGainDeletion> driverGermlineAmpDels = null;

        if(purpleData.germlineAmpDels() != null)
        {
            driverGermlineAmpDels = GermlineGainDeletionFactory.createGermlineGainDeletions(
                    purpleData.germlineAmpDels(), Objects.requireNonNull(germlineDrivers), purpleData.somaticGeneCopyNumbers(), isofoxData);
        }

        double ploidy = purpleData.purityContext().bestFit().ploidy();

        List<PurpleChrArmCopyNumber> armCopyNumberAbberations = purpleData.chrArmCopyNumbers().stream()
                .map(x -> PurpleConversion.convert(x, ploidy))
                .filter(x -> !x.type().equals(PURPLE_ARM_CN_DIPLOID))
                .collect(Collectors.toList());

        return ImmutablePurpleRecord.builder()
                .fit(createFit(purpleData))
                .characteristics(createCharacteristics(purpleData))
                .somaticDrivers(ConversionUtil.mapToIterable(purpleData.somaticDrivers(), PurpleConversion::convert))
                .germlineDrivers(germlineDrivers)
                .somaticVariants(somaticVariants)
                .germlineVariants(germlineVariants)
                .armCopyNumberAbberations(armCopyNumberAbberations)
                .somaticGainsDels(driverSomaticGainsDels)
                .germlineGainsDels(driverGermlineAmpDels)
                .build();
    }

    private static List<PurpleVariant> buildPurpleVariants(
            final List<SmallVariant> variants, final List<DriverCatalog> drivers, boolean isGermline)
    {
        if(variants == null)
            return null;

        List<PurpleVariant> purpleVariants = Lists.newArrayListWithCapacity(variants.size());

        DriverType requiredDriverType = isGermline ? GERMLINE_MUTATION : MUTATION;

        for(SmallVariant variant : variants)
        {
            DriverCatalog driver = drivers.stream().filter(x -> x.driver() == requiredDriverType).findFirst().orElse(null);
            PurpleVariant purpleVariant = PurpleVariantFactory.buildPurpleVariant(variant, driver, isGermline);
            purpleVariants.add(purpleVariant);
        }

        return purpleVariants;
    }

    private static final Set<DriverType> AMP_DEL_TYPES = EnumSet.of(
            DriverType.AMP, DriverType.PARTIAL_AMP, DriverType.DEL, DriverType.HET_DEL, DriverType.LOH);

    private static List<PurpleGainDeletion> somaticGainsDelsFromDrivers(
            final List<DriverCatalog> drivers, final List<GeneCopyNumber> geneCopyNumbers, final LinxData linxData,
            @Nullable final IsofoxData isofoxData)
    {
        List<PurpleGainDeletion> gainDeletions = Lists.newArrayList();

        for(DriverCatalog driver : drivers)
        {
            if(AMP_DEL_TYPES.contains(driver.driver()))
            {
                GeneCopyNumber geneCopyNumber = geneCopyNumbers.stream()
                        .filter(x -> x.GeneName.equals(driver.gene()) && driver.transcript().equals(x.TransName))
                        .findFirst().orElse(null);

                GeneExpression geneExpression = isofoxData != null ? isofoxData.geneExpressions().stream()
                        .filter(x -> x.geneName().equals(driver.gene())).findFirst().orElse(null) : null;

                gainDeletions.add(toGainDel(driver, geneCopyNumber, linxData, geneExpression));
            }
        }

        return gainDeletions;
    }

    private static PurpleGainDeletion toGainDel(
            final DriverCatalog driver, final GeneCopyNumber geneCopyNumber, final LinxData linxData,
            @Nullable final GeneExpression geneExpression)
    {
        Double tpm = null;
        Double tpmPercentile = null;
        Double tpmFoldChange = null;

        if(geneExpression != null)
        {
            tpm = geneExpression.tpm();
            tpmPercentile = geneExpression.percentileCohort();
            tpmFoldChange = geneExpression.medianTpmCohort() > 0 ? tpm / geneExpression.medianTpmCohort() : 0;
        }

        String geneRange = PURPLE_AMP_DEL_PARTIAL;

        if(driver.driver() == DriverType.AMP)
        {
            geneRange = PURPLE_AMP_DEL_FULL;
        }

        String plotFilename = findLinxDriverPlot(driver, linxData);

        return ImmutablePurpleGainDeletion.builder()
                .driver(PurpleConversion.convert(driver))
                .chromosome(driver.chromosome())
                .chromosomeBand(driver.chromosomeBand())
                .minCopyNumber(Math.max(0, driver.minCopyNumber()))
                .maxCopyNumber(Math.max(0, driver.maxCopyNumber()))
                .geneRange(geneRange)
                .exonStart(null)
                .exonEnd(null)
                .relativeCopyNumber(geneCopyNumber.RelativeMinCopyNumber)
                .minMinorAlleleCopies(geneCopyNumber.MinMinorAlleleCopyNumber)
                .tpm(tpm)
                .tpmPercentile(tpmPercentile)
                .tpmFoldChange(tpmFoldChange)
                .plotFilename(plotFilename)
                .build();
    }

    private static final Set<DriverEventType> AMP_LINX_DRIVER_TYPES = EnumSet.of(GAIN, GAIN_ARM, GAIN_CHR);
    private static final Set<DriverEventType> DEL_LINX_DRIVER_TYPES = EnumSet.of(DEL, LOH, LOH_SV_CENTRO, LOH_SV_TELO, LOH_ARM, LOH_CHR);

    private static String findLinxDriverPlot(final DriverCatalog driver, final LinxData linxData)
    {
        Integer clusterId = null;

        for(LinxDriver linxDriver : linxData.somaticDriverData())
        {
            if(!linxDriver.gene().equals(driver.gene()))
                continue;

            if(driver.driver() == AMP || driver.driver() == PARTIAL_AMP)
            {
                if(AMP_LINX_DRIVER_TYPES.contains(linxDriver.eventType()))
                {
                    clusterId = linxDriver.clusterId();
                    break;
                }
            }
            else if(driver.driver() == DriverType.DEL)
            {
                if(DEL_LINX_DRIVER_TYPES.contains(linxDriver.eventType()))
                {
                    clusterId = linxDriver.clusterId();
                    break;
                }
            }
        }

        return clusterId != null ? findReportableLinxPlot(linxData.reportableEventPlots(), clusterId) : null;
    }

    private static PurpleFit createFit(final PurpleData purple)
    {
        return ImmutablePurpleFit.builder()
                .qc(PurpleConversion.convert(purple.purityContext().qc()))
                .fittedPurityMethod(PurpleFittedPurityMethod.valueOf(purple.purityContext().method().name()))
                .purity(purple.purityContext().bestFit().purity())
                .minPurity(purple.purityContext().score().minPurity())
                .maxPurity(purple.purityContext().score().maxPurity())
                .ploidy(purple.purityContext().bestFit().ploidy())
                .minPloidy(purple.purityContext().score().minPloidy())
                .maxPloidy(purple.purityContext().score().maxPloidy())
                .build();
    }

    private static PurpleCharacteristics createCharacteristics(final PurpleData purple)
    {
        return ImmutablePurpleCharacteristics.builder()
                .wholeGenomeDuplication(purple.purityContext().wholeGenomeDuplication())
                .microsatelliteIndelsPerMb(purple.purityContext().microsatelliteIndelsPerMb())
                .microsatelliteStatus(PurpleMicrosatelliteStatus.valueOf(purple.purityContext().microsatelliteStatus().name()))
                .tumorMutationalBurdenPerMb(purple.purityContext().tumorMutationalBurdenPerMb())
                .tumorMutationalBurdenStatus(PurpleTumorMutationalStatus.valueOf(purple.purityContext()
                        .tumorMutationalBurdenStatus()
                        .name()))
                .tumorMutationalLoad(purple.purityContext().tumorMutationalLoad())
                .tumorMutationalLoadStatus(PurpleTumorMutationalStatus.valueOf(purple.purityContext().tumorMutationalLoadStatus().name()))
                .svTumorMutationalBurden(purple.purityContext().svTumorMutationalBurden())
                .lohPercentage(purple.purityContext().qc().lohPercent())
                .build();
    }
}
