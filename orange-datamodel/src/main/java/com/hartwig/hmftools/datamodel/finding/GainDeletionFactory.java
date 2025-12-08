package com.hartwig.hmftools.datamodel.finding;

import java.util.List;
import java.util.Optional;

import com.hartwig.hmftools.datamodel.driver.DriverInterpretation;
import com.hartwig.hmftools.datamodel.driver.ReportedStatus;
import com.hartwig.hmftools.datamodel.purple.PurpleDriver;
import com.hartwig.hmftools.datamodel.purple.PurpleDriverType;
import com.hartwig.hmftools.datamodel.purple.PurpleGainDeletion;

public class GainDeletionFactory {

    public static List<GainDeletion> somaticGainsDelsFromDrivers(List<PurpleGainDeletion> gainDeletions, final List<PurpleDriver> drivers)
    {
        List<GainDeletion> somaticGainsDels = new java.util.ArrayList<>();
        for(PurpleGainDeletion gainDeletion : gainDeletions)
        {
            // we have to reverse the copy number interpretation logic to get back the purple driver type
            final PurpleDriverType type = switch(gainDeletion.interpretation())
            {
                case FULL_GAIN -> PurpleDriverType.AMP;
                case PARTIAL_GAIN -> PurpleDriverType.PARTIAL_AMP;
                case FULL_DEL, PARTIAL_DEL -> PurpleDriverType.DEL;
            };

            Optional<PurpleDriver> driver = drivers.stream()
                    .filter(o -> o.gene().equals(gainDeletion.gene()) && o.transcript().equals(gainDeletion.transcript()) && o.type().equals(type))
                    .findFirst();

            driver.ifPresent(purpleDriver -> somaticGainsDels.add(toGainDel(gainDeletion, purpleDriver, FindingKeys.SampleType.SOMATIC)));
        }
        return somaticGainsDels;
    }

    public static List<GainDeletion> convertGermlineFullDels(List<PurpleGainDeletion> gainDeletions, final List<PurpleDriver> drivers)
    {
        List<GainDeletion> germlineDels = new java.util.ArrayList<>();
        for(PurpleGainDeletion gainDeletion : gainDeletions)
        {
            Optional<PurpleDriver> driver = drivers.stream()
                    .filter(o -> o.gene().equals(gainDeletion.gene()) && o.transcript().equals(gainDeletion.transcript()) && o.type().equals(PurpleDriverType.GERMLINE_DELETION))
                    .findFirst();

            driver.ifPresent(purpleDriver -> germlineDels.add(toGainDel(gainDeletion, purpleDriver, FindingKeys.SampleType.GERMLINE)));
        }
        return germlineDels;
    }

    private static GainDeletion toGainDel(PurpleGainDeletion purpleGainDeletion, final PurpleDriver driver, FindingKeys.SampleType sampleType) {
        return ImmutableGainDeletion.builder()
                .findingKey(FindingKeys.gainDeletion(sampleType,
                        driver.gene(),
                        purpleGainDeletion.interpretation(),
                        driver.isCanonical(),
                        driver.transcript()))
                .reportedStatus(ReportedStatus.REPORTED) // fix this later
                .driverInterpretation(DriverInterpretation.interpret(driver.driverLikelihood()))
                .driver(driver)
                .chromosome(purpleGainDeletion.chromosome())
                .chromosomeBand(purpleGainDeletion.chromosomeBand())
                .gene(driver.gene())
                .transcript(driver.transcript())
                .isCanonical(driver.isCanonical())
                .interpretation(purpleGainDeletion.interpretation())
                .minCopies(purpleGainDeletion.minCopies())
                .maxCopies(purpleGainDeletion.maxCopies())
                .build();
    }
}