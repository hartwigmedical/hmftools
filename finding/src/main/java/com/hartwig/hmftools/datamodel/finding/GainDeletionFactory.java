package com.hartwig.hmftools.datamodel.finding;

import java.util.ArrayList;
import java.util.List;

import com.hartwig.hmftools.datamodel.driver.DriverInterpretation;
import com.hartwig.hmftools.datamodel.driver.DriverSource;
import com.hartwig.hmftools.datamodel.driver.ReportedStatus;
import com.hartwig.hmftools.datamodel.purple.CopyNumberInterpretation;
import com.hartwig.hmftools.datamodel.purple.PurpleDriver;
import com.hartwig.hmftools.datamodel.purple.PurpleDriverType;
import com.hartwig.hmftools.datamodel.purple.PurpleGainDeletion;
import com.hartwig.hmftools.datamodel.purple.PurpleLossOfHeterozygosity;
import com.hartwig.hmftools.datamodel.purple.PurpleRecord;

import org.jetbrains.annotations.NotNull;

final class GainDeletionFactory {

    public static DriverFindings<GainDeletion> gainDeletionFindings(@NotNull PurpleRecord purple, @NotNull FindingsStatus findingsStatus) {
        List<GainDeletion> allGainDels = new ArrayList<>();
        List<PurpleGainDeletion> germlineFullDels = purple.reportableGermlineFullDels();
        List<PurpleLossOfHeterozygosity> germlineLohs = purple.reportableGermlineLossOfHeterozygosities();
        List<PurpleDriver> purpleGermlineDrivers = purple.germlineDrivers();
        if (germlineFullDels != null && germlineLohs != null && purpleGermlineDrivers != null) {
            allGainDels.addAll(germlineDriverGainDels(germlineFullDels, germlineLohs, purple.germlineDrivers()));
        }
        allGainDels.addAll(somaticDriverGainDels(purple.reportableSomaticGainsDels(), purple.somaticDrivers()));
        return ImmutableDriverFindings.<GainDeletion>builder()
                .status(findingsStatus)
                .all(allGainDels)
                .build();
    }

    // in orange data, HOM_DELS are stored as germline full dels, HET_DELS are stored in LOH, they do not overlap.
    // all the reportable ones are in purple drivers. Other types are not reportable, we can ignore them
    public static List<GainDeletion> germlineDriverGainDels(List<PurpleGainDeletion> reportableGermlineFullDels,
            List<PurpleLossOfHeterozygosity> reportableGermlineLossOfHeterozygosities,
            final List<PurpleDriver> germlineDrivers)
    {
        List<GainDeletion> driverGainDels = new ArrayList<>();

        for(PurpleGainDeletion fullDels : reportableGermlineFullDels)
        {
            // find the purple driver object, it should be there
            PurpleDriver driver = findDriver(germlineDrivers, fullDels.gene(), fullDels.transcript(), PurpleDriverType.GERMLINE_DELETION);

            driverGainDels.add(toGainDel(fullDels, driver, GainDeletion.Type.GERMLINE_DEL_HOM_IN_TUMOR, DriverSource.GERMLINE));
        }

        for(PurpleLossOfHeterozygosity loh : reportableGermlineLossOfHeterozygosities)
        {
            // find the purple driver object, it should be there
            PurpleDriver driver = findDriver(germlineDrivers, loh.gene(), loh.transcript(), PurpleDriverType.GERMLINE_DELETION);

            driverGainDels.add(toGainDel(loh, driver));
        }

        // should we sort this?
        return driverGainDels;
    }

    public static List<GainDeletion> somaticDriverGainDels(List<PurpleGainDeletion> gainDeletions, final List<PurpleDriver> drivers)
    {
        List<GainDeletion> somaticGainsDels = new ArrayList<>();
        for(PurpleGainDeletion gainDeletion : gainDeletions)
        {
            // we have to reverse the copy number interpretation logic to get back the purple driver type
            final PurpleDriverType purpleDriverType = switch(gainDeletion.interpretation())
            {
                case FULL_GAIN -> PurpleDriverType.AMP;
                case PARTIAL_GAIN -> PurpleDriverType.PARTIAL_AMP;
                case FULL_DEL, PARTIAL_DEL -> PurpleDriverType.DEL;
            };

            PurpleDriver driver = findDriver(drivers, gainDeletion.gene(), gainDeletion.transcript(), purpleDriverType);

            final GainDeletion.Type type = switch(gainDeletion.interpretation())
            {
                case FULL_GAIN, PARTIAL_GAIN -> GainDeletion.Type.SOMATIC_GAIN;
                case FULL_DEL, PARTIAL_DEL -> GainDeletion.Type.SOMATIC_DEL;
            };

            somaticGainsDels.add(toGainDel(gainDeletion, driver, type, DriverSource.SOMATIC));
        }
        return somaticGainsDels;
    }

    private static PurpleDriver findDriver(final List<PurpleDriver> drivers, final String gene, final String transcript,
            final PurpleDriverType purpleDriverType) {
        return drivers.stream()
                .filter(o -> o.gene().equals(gene) && o.transcript().equals(transcript) && o.type().equals(purpleDriverType))
                .findFirst()
                .orElseThrow(() -> new IllegalStateException("No driver found for " + gene + " transcript " + transcript + " type " + purpleDriverType));
    }

    private static GainDeletion toGainDel(PurpleGainDeletion purpleGainDeletion,
            final PurpleDriver driver,
            GainDeletion.Type type,
            DriverSource sourceSample) {
        return ImmutableGainDeletion.builder()
                .findingKey(FindingKeys.gainDeletion(sourceSample,
                        purpleGainDeletion.gene(),
                        purpleGainDeletion.interpretation(),
                        driver.isCanonical(),
                        purpleGainDeletion.transcript()))
                .driverSource(sourceSample)
                .reportedStatus(ReportedStatus.REPORTED)
                .driverInterpretation(DriverInterpretation.interpret(driver.driverLikelihood()))
                .type(type)
                .chromosome(purpleGainDeletion.chromosome())
                .chromosomeBand(purpleGainDeletion.chromosomeBand())
                .gene(purpleGainDeletion.gene())
                .transcript(purpleGainDeletion.transcript())
                .isCanonical(driver.isCanonical())
                .interpretation(purpleGainDeletion.interpretation())
                .minCopies(purpleGainDeletion.minCopies())
                .maxCopies(purpleGainDeletion.maxCopies())
                .build();
    }

    private static GainDeletion toGainDel(PurpleLossOfHeterozygosity loh, final PurpleDriver driver) {

        CopyNumberInterpretation copyNumberInterpretation = switch (loh.geneProportion())
        {
            case FULL_GENE -> CopyNumberInterpretation.FULL_GAIN;
            case PARTIAL_GENE -> CopyNumberInterpretation.PARTIAL_DEL;
        };

        return ImmutableGainDeletion.builder()
                .findingKey(FindingKeys.gainDeletion(DriverSource.GERMLINE,
                        loh.gene(),
                        copyNumberInterpretation,
                        driver.isCanonical(),
                        loh.transcript()))
                .driverSource(DriverSource.GERMLINE)
                .reportedStatus(ReportedStatus.REPORTED)
                .driverInterpretation(DriverInterpretation.interpret(driver.driverLikelihood()))
                .type(GainDeletion.Type.GERMLINE_DEL_HET_IN_TUMOR)
                .chromosome(loh.chromosome())
                .chromosomeBand(loh.chromosomeBand())
                .gene(loh.gene())
                .transcript(loh.transcript())
                .isCanonical(driver.isCanonical())
                .interpretation(copyNumberInterpretation)
                .minCopies(loh.minCopies())
                .maxCopies(loh.maxCopies())
                .build();
    }
}