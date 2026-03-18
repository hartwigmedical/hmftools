package com.hartwig.hmftools.purple.plot;

import static java.lang.String.format;
import static java.util.Comparator.comparingInt;
import static java.util.stream.Collectors.groupingBy;
import static java.util.stream.Collectors.minBy;

import static com.hartwig.hmftools.common.driver.DriverInterpretation.LOW;
import static com.hartwig.hmftools.common.driver.DriverType.DRIVERS_PURPLE_SOMATIC;
import static com.hartwig.hmftools.common.variant.VariantType.INDEL;
import static com.hartwig.hmftools.purple.PurpleUtils.PPL_LOGGER;
import static com.hartwig.hmftools.purple.plot.CircosFileWriter.circosContig;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.Set;
import java.util.function.Function;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.driver.DriverInterpretation;
import com.hartwig.hmftools.common.driver.DriverType;
import com.hartwig.hmftools.common.purple.GeneCopyNumber;
import com.hartwig.hmftools.common.purple.ReportedStatus;
import com.hartwig.hmftools.common.utils.Doubles;
import com.hartwig.hmftools.purple.DriverSourceData;
import com.hartwig.hmftools.purple.somatic.SomaticVariant;

public final class CircosDriverWriter
{
    public static void writeDrivers(
            final String driverTextFilePath, final String driverLineFilePath, final List<DriverSourceData> driverDataList)
            throws IOException
    {
        List<DriverSourceData> filteredDriverSourceData = filterAndPrioritise(driverDataList);

        Files.write(new File(driverTextFilePath).toPath(), drivers(filteredDriverSourceData, CircosDriverWriter::driverText));
        Files.write(new File(driverLineFilePath).toPath(), drivers(filteredDriverSourceData, CircosDriverWriter::driverPointer));
    }

    private static List<DriverSourceData> filterAndPrioritise(final List<DriverSourceData> driverDataList)
    {
        // show reported events and only select one event per gene, prioritising copy-number events
        List<DriverSourceData> filteredDriverData = driverDataList.stream()
                .filter(x -> x.DriverData.reportedStatus() == ReportedStatus.REPORTED)
                .filter(x -> DriverInterpretation.interpret(x.DriverData.driverLikelihood()) != LOW)
                .filter(x -> DRIVERS_PURPLE_SOMATIC.contains(x.DriverData.driver()))
                .collect(Collectors.toList());

        Collections.sort(filteredDriverData, new DriverSourceDataSorter());

        Set<String> processedGenes = Sets.newHashSet();

        int index = 0;

        while(index < filteredDriverData.size())
        {
            DriverSourceData driver = filteredDriverData.get(index);

            if(processedGenes.contains(driver.DriverData.gene()))
            {
                filteredDriverData.remove(index);
                continue;
            }

            processedGenes.add(driver.DriverData.gene());
            ++index;
        }

        return filteredDriverData;
    }

    private static class DriverSourceDataSorter implements Comparator<DriverSourceData>
    {
        public int compare(final DriverSourceData first, final DriverSourceData second)
        {
            double likelihood1 = first.DriverData.driverLikelihood();
            double likelihood2 = second.DriverData.driverLikelihood();

            if(likelihood1 != likelihood2)
                return likelihood1 > likelihood2 ? -1 : 1;

            int driverRank1 = driverTypeRank(first.DriverData.driver());
            int driverRank2 = driverTypeRank(second.DriverData.driver());

            return Integer.compare(driverRank1, driverRank2);
        }
    }

    private static int driverTypeRank(final DriverType driverType)
    {
        switch(driverType)
        {
            case AMP:
            case PARTIAL_AMP:
            case DEL:
                return 0;

            case LOH:
            case HET_DEL:
                return 1;

            case MUTATION:
                return 2;

            default:
                return 3;
        }
    }

    private static List<String> drivers(
            final List<DriverSourceData> driverDataList, final Function<DriverSourceData, String> driverSourceDataToLine)
    {
        return driverDataList.stream()
                .map(driverSourceDataToLine)
                .toList();
    }

    private static String driverText(final DriverSourceData driverData)
    {
        if(driverData.SourceObject instanceof final GeneCopyNumber geneCopyNumber)
        {
            //PPL_LOGGER.info("CNV gene: {}, driver: {}, biallellic: {}", geneCopyNumber.geneName(), o.DriverData.driver(),
            //        o.DriverData.biallelic());

            int geneRegionMid = (geneCopyNumber.start() + geneCopyNumber.end()) / 2;
            return String.join("\t",
                    circosContig(geneCopyNumber.chromosome()),
                    String.valueOf(geneRegionMid),
                    String.valueOf(geneRegionMid),
                    geneCopyNumber.geneName(),
                    format("color=%s", driverData.DriverData.driver() == DriverType.DEL ? "vdred" : "vdgreen"));
        }
        else if(driverData.SourceObject instanceof final SomaticVariant somaticVariant)
        {
            return String.join("\t",
                    circosContig(somaticVariant.chromosome()),
                    String.valueOf(somaticVariant.position()),
                    String.valueOf(somaticVariant.position()),
                    somaticVariant.gene(),
                    format("color=black"));
        }
        throw new IllegalStateException("Unexpected value: " + driverData.SourceObject);
    }

    // from circos doc:
    // Each highlight line has five dimensions:
    // - outer padding
    // - outer line length (drawn at new position)
    // - connecting line length (connects old to new position)
    // - inner line length (drawn at old position)
    // - inner padding
    private static String driverPointer(final DriverSourceData driverData)
    {
        if(driverData.SourceObject instanceof final GeneCopyNumber geneCopyNumber)
        {
            double r0 = switch(driverData.DriverData.driver())
            {
                case AMP, PARTIAL_AMP -> 0.7;
                case DEL -> 0.55;
                default -> 0.5;
            };

            int geneRegionMid = (geneCopyNumber.start() + geneCopyNumber.end()) / 2;

            return String.join("\t",
                    circosContig(geneCopyNumber.chromosome()),
                    String.valueOf(geneRegionMid),
                    String.valueOf(geneRegionMid),
                    geneCopyNumber.geneName(),
                    format("r0=%.2fr", r0));
        }
        else if(driverData.SourceObject instanceof final SomaticVariant somaticVariant)
        {
            double r0;

            if(somaticVariant.type() == INDEL)
            {
                r0 = 0.775;
            }
            else
            {
                r0 = 0.775 + Doubles.clamp(somaticVariant.decorator().adjustedVaf(), 0, 1.0) * (0.975 - 0.775);
            }
            return String.join("\t",
                    circosContig(somaticVariant.chromosome()),
                    String.valueOf(somaticVariant.position()),
                    String.valueOf(somaticVariant.position()),
                    somaticVariant.gene(),
                    format("r0=%.2fr", r0));
        }
        throw new IllegalStateException("Unexpected value: " + driverData.SourceObject);
    }
}
