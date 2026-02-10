package com.hartwig.hmftools.orange.algo.purple;

import static com.hartwig.hmftools.orange.OrangeApplication.LOGGER;
import static com.hartwig.hmftools.orange.algo.purple.GermlineGainDeletionUtil.findGeneCopyNumberForGeneTranscript;

import java.util.EnumSet;
import java.util.List;
import java.util.Objects;
import java.util.Set;

import com.hartwig.hmftools.common.driver.DriverCatalog;
import com.hartwig.hmftools.common.driver.DriverType;
import com.hartwig.hmftools.common.purple.GeneCopyNumber;
import com.hartwig.hmftools.datamodel.purple.ImmutablePurpleCharacteristics;
import com.hartwig.hmftools.datamodel.purple.ImmutablePurpleFit;
import com.hartwig.hmftools.datamodel.purple.ImmutablePurpleGainDeletion;
import com.hartwig.hmftools.datamodel.purple.ImmutablePurpleRecord;
import com.hartwig.hmftools.datamodel.purple.PurpleCharacteristics;
import com.hartwig.hmftools.datamodel.purple.PurpleDriver;
import com.hartwig.hmftools.datamodel.purple.PurpleFit;
import com.hartwig.hmftools.datamodel.purple.PurpleFittedPurityMethod;
import com.hartwig.hmftools.datamodel.purple.PurpleGainDeletion;
import com.hartwig.hmftools.datamodel.purple.PurpleMicrosatelliteStatus;
import com.hartwig.hmftools.datamodel.purple.PurpleRecord;
import com.hartwig.hmftools.datamodel.purple.PurpleTumorMutationalStatus;
import com.hartwig.hmftools.datamodel.purple.PurpleVariant;
import com.hartwig.hmftools.orange.conversion.ConversionUtil;
import com.hartwig.hmftools.orange.conversion.PurpleConversion;

import org.jetbrains.annotations.Nullable;

public class PurpleInterpreter
{
    private final GermlineGainDeletionFactory mGermlineGainDelFactory;

    public PurpleInterpreter(final GermlineGainDeletionFactory germlineGainDeletionFactory)
    {
        mGermlineGainDelFactory = germlineGainDeletionFactory;
    }

    public PurpleRecord interpret(final PurpleData purple)
    {
        LOGGER.info("Analysing purple data");

        List<PurpleVariant> somaticVariants = PurpleVariantFactory.fromPurpleVariants(purple.somaticVariants());

        List<PurpleVariant> germlineVariants = PurpleVariantFactory.fromPurpleVariants(purple.germlineVariants());

        @Nullable List<PurpleDriver> germlineDrivers = ConversionUtil.mapToNullableList(purple.germlineDrivers(), PurpleConversion::convert);

        List<PurpleGainDeletion> driverSomaticGainsDels = somaticGainsDelsFromDrivers(purple.somaticDrivers(), purple.somaticGeneCopyNumbers());

        List<PurpleGainDeletion> driverGermlineAmpDels = null;

        if(purple.germlineDeletions() != null)
        {
            driverGermlineAmpDels = mGermlineGainDelFactory.createGermlineGainDeletions(
                    purple.germlineDeletions(), Objects.requireNonNull(germlineDrivers), purple.somaticGeneCopyNumbers());
        }

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
                .somaticGainsDels(driverSomaticGainsDels)
                .germlineGainsDels(driverGermlineAmpDels)
                .build();
    }

    private static final Set<DriverType> AMP_DEL_TYPES = EnumSet.of(
            DriverType.AMP, DriverType.PARTIAL_AMP, DriverType.DEL, DriverType.HET_DEL, DriverType.LOH);

    private static List<PurpleGainDeletion> somaticGainsDelsFromDrivers(
            final List<DriverCatalog> drivers, final List<GeneCopyNumber> allSomaticGeneCopyNumbers)
    {
        return drivers.stream()
                .filter(o -> AMP_DEL_TYPES.contains(o.driver()))
                .map(o -> toGainDel(
                        o, findGeneCopyNumberForGeneTranscript(o.gene(), o.transcript(), allSomaticGeneCopyNumbers)))
                .toList();
    }

    private static PurpleGainDeletion toGainDel(final DriverCatalog driver, final GeneCopyNumber somaticGeneCopyNumber)
    {
        return ImmutablePurpleGainDeletion.builder()
                .driver(PurpleConversion.convert(driver))
                .interpretation(CopyNumberInterpretationUtil.fromCNADriver(driver))
                .chromosome(driver.chromosome())
                .chromosomeBand(driver.chromosomeBand())
                .minCopies(Math.max(0, driver.minCopyNumber()))
                .maxCopies(Math.max(0, driver.maxCopyNumber()))
                .minMinorAlleleCopies(somaticGeneCopyNumber.MinMinorAlleleCopyNumber)
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
