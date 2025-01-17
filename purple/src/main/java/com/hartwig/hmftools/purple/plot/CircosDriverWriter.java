package com.hartwig.hmftools.purple.plot;

import static java.lang.String.format;
import static java.util.Comparator.comparingInt;
import static java.util.stream.Collectors.groupingBy;
import static java.util.stream.Collectors.minBy;

import static com.hartwig.hmftools.common.variant.VariantType.INDEL;
import static com.hartwig.hmftools.purple.PurpleUtils.PPL_LOGGER;
import static com.hartwig.hmftools.purple.plot.CircosFileWriter.circosContig;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.function.Function;

import com.hartwig.hmftools.common.driver.DriverType;
import com.hartwig.hmftools.common.purple.GeneCopyNumber;
import com.hartwig.hmftools.common.utils.Doubles;
import com.hartwig.hmftools.purple.DriverSourceData;
import com.hartwig.hmftools.purple.somatic.SomaticVariant;

import org.jetbrains.annotations.NotNull;

public final class CircosDriverWriter
{
    public static void writeDrivers(final String driverTextFilePath, final String driverLineFilePath,
            final List<DriverSourceData> driverSourceData)
            throws IOException
    {
        List<DriverSourceData> filteredDriverSourceData = filterAndDropDup(driverSourceData);

        Files.write(new File(driverTextFilePath).toPath(), drivers(filteredDriverSourceData, CircosDriverWriter::driverText));
        Files.write(new File(driverLineFilePath).toPath(), drivers(filteredDriverSourceData, CircosDriverWriter::driverPointer));
    }

    // rank them
    private static List<DriverSourceData> filterAndDropDup(final List<DriverSourceData> driverSourceData)
    {
        // First filter out germline ones, and driverlikelyhood < 0.8
        // groupby gene and find the one with highest rank
        Map<String, Optional<DriverSourceData>> m = driverSourceData.stream()
                //.filter(o -> !DriverType.isGermline(o.DriverData.driver()) && o.DriverData.driverLikelihood() >= 0.8)
                .collect(groupingBy(o -> o.DriverData.gene(), minBy(comparingInt(CircosDriverWriter::driverRank))));

        // get all the resulting driver source data object into list
        return m.values().stream().flatMap(Optional::stream).toList();
    }

    private static int driverRank(DriverSourceData driverSourceData)
    {
        if(driverSourceData.SourceObject instanceof GeneCopyNumber)
        {
            // copy number changes are ranked highest
            return 1;
        }
        if(driverSourceData.SourceObject instanceof final SomaticVariant somaticVariant)
        {
            if(somaticVariant.type() == INDEL)
            {
                // indel is ranked higher than SNP
                return 2;
            }
            else
            {
                return 3;
            }
        }
        return 4;
    }

    private static List<String> drivers(final List<DriverSourceData> driverSourceData,
            final Function<DriverSourceData, String> driverSourceDataToLine)
    {
        return driverSourceData.stream()
                .map(driverSourceDataToLine)
                .toList();
    }

    @NotNull
    private static String driverText(final DriverSourceData o)
    {
        if(o.SourceObject instanceof final GeneCopyNumber geneCopyNumber)
        {
            PPL_LOGGER.info("CNV gene: {}, driver: {}, biallellic: {}", geneCopyNumber.geneName(), o.DriverData.driver(),
                    o.DriverData.biallelic());

            int geneRegionMid = (geneCopyNumber.start() + geneCopyNumber.end()) / 2;
            return String.join("\t",
                    circosContig(geneCopyNumber.chromosome()),
                    String.valueOf(geneRegionMid),
                    String.valueOf(geneRegionMid),
                    geneCopyNumber.geneName(),
                    format("color=%s", o.DriverData.driver() == DriverType.DEL ? "vdred" : "vdgreen"));
        }
        else if(o.SourceObject instanceof final SomaticVariant somaticVariant)
        {
            PPL_LOGGER.info("SV gene: {}, driver: {}, biallellic: {}", somaticVariant.gene(), o.DriverData.driver(),
                    o.DriverData.biallelic());

            return String.join("\t",
                    circosContig(somaticVariant.chromosome()),
                    String.valueOf(somaticVariant.position()),
                    String.valueOf(somaticVariant.position()),
                    somaticVariant.gene(),
                    format("color=black"));
        }
        throw new IllegalStateException("Unexpected value: " + o.SourceObject);
    }

    // from circos doc:
    // Each highlight line has five dimensions:
    // - outer padding
    // - outer line length (drawn at new position)
    // - connecting line length (connects old to new position)
    // - inner line length (drawn at old position)
    // - inner padding
    private static String driverPointer(final DriverSourceData o)
    {
        if(o.SourceObject instanceof final GeneCopyNumber geneCopyNumber)
        {
            double r0 = switch(o.DriverData.driver())
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
        else if(o.SourceObject instanceof final SomaticVariant somaticVariant)
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
        throw new IllegalStateException("Unexpected value: " + o.SourceObject);
    }
}
