package com.hartwig.hmftools.orange.algo.purple;

import static com.hartwig.hmftools.orange.OrangeApplication.LOGGER;

import java.util.EnumSet;
import java.util.List;
import java.util.Objects;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.driver.DriverCatalog;
import com.hartwig.hmftools.common.driver.DriverType;
import com.hartwig.hmftools.common.purple.GeneCopyNumber;
import com.hartwig.hmftools.common.rna.GeneExpression;
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
import com.hartwig.hmftools.orange.conversion.ConversionUtil;
import com.hartwig.hmftools.orange.conversion.PurpleConversion;

import org.jetbrains.annotations.Nullable;

public class PurpleInterpreter
{
    public PurpleInterpreter() {}

    public PurpleRecord interpret(final PurpleData purple, @Nullable final IsofoxData isofoxData)
    {
        LOGGER.info("Analysing Purple data");

        List<PurpleVariant> somaticVariants = PurpleVariantFactory.fromPurpleVariants(purple.somaticVariants(), purple.somaticDrivers());

        List<PurpleVariant> germlineVariants = PurpleVariantFactory.fromPurpleVariants(purple.germlineVariants(), purple.germlineDrivers());

        List<PurpleDriver> germlineDrivers = ConversionUtil.mapToNullableList(purple.germlineDrivers(), PurpleConversion::convert);

        List<PurpleGainDeletion> driverSomaticGainsDels = somaticGainsDelsFromDrivers(
                purple.somaticDrivers(), purple.somaticGeneCopyNumbers(), isofoxData);

        List<PurpleGainDeletion> driverGermlineAmpDels = null;

        if(purple.germlineDeletions() != null)
        {
            driverGermlineAmpDels = GermlineGainDeletionFactory.createGermlineGainDeletions(
                    purple.germlineDeletions(), Objects.requireNonNull(germlineDrivers), purple.somaticGeneCopyNumbers(), isofoxData);
        }

        double ploidy = purple.purityContext().bestFit().ploidy();

        List<PurpleChrArmCopyNumber> chrArmCopyNumbers = purple.chrArmCopyNumbers().stream()
                .map(x -> PurpleConversion.convert(x, ploidy)).collect(Collectors.toList());

        return ImmutablePurpleRecord.builder()
                .fit(createFit(purple))
                .tumorStats(TumorStatsFactory.compute(purple))
                .characteristics(createCharacteristics(purple))
                .somaticDrivers(ConversionUtil.mapToIterable(purple.somaticDrivers(), PurpleConversion::convert))
                .germlineDrivers(germlineDrivers)
                .somaticVariants(somaticVariants)
                .germlineVariants(germlineVariants)
                .somaticCopyNumbers(ConversionUtil.mapToIterable(purple.somaticCopyNumbers(), PurpleConversion::convert))
                .somaticGeneCopyNumbers(ConversionUtil.mapToIterable(purple.somaticGeneCopyNumbers(), PurpleConversion::convert))
                .chrArmCopyNumbers(chrArmCopyNumbers)
                .somaticGainsDels(driverSomaticGainsDels)
                .germlineGainsDels(driverGermlineAmpDels)
                .build();
    }

    private static final Set<DriverType> AMP_DEL_TYPES = EnumSet.of(
            DriverType.AMP, DriverType.PARTIAL_AMP, DriverType.DEL, DriverType.HET_DEL, DriverType.LOH);

    private static List<PurpleGainDeletion> somaticGainsDelsFromDrivers(
            final List<DriverCatalog> drivers, final List<GeneCopyNumber> geneCopyNumbers, @Nullable final IsofoxData isofoxData)
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

                gainDeletions.add(toGainDel(driver, geneCopyNumber, geneExpression));
            }
        }

        return gainDeletions;
    }

    private static PurpleGainDeletion toGainDel(
            final DriverCatalog driver, final GeneCopyNumber geneCopyNumber, @Nullable final GeneExpression geneExpression)
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

        String exonRange = "PARTIAL";

        if(driver.driver() == DriverType.AMP)
        {
            exonRange = "FULL";
        }

        return ImmutablePurpleGainDeletion.builder()
                .driver(PurpleConversion.convert(driver))
                .chromosome(driver.chromosome())
                .chromosomeBand(driver.chromosomeBand())
                .minCopyNumber(Math.max(0, driver.minCopyNumber()))
                .maxCopyNumber(Math.max(0, driver.maxCopyNumber()))
                .exonRange(exonRange)
                .relativeCopyNumber(geneCopyNumber.RelativeMinCopyNumber)
                .minMinorAlleleCopies(geneCopyNumber.MinMinorAlleleCopyNumber)
                .tpm(tpm)
                .tpmPercentile(tpmPercentile)
                .tpmFoldChange(tpmFoldChange)
                .build();
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
                .build();
    }
}
