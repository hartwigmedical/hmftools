package com.hartwig.hmftools.common.driver.panel;

import static com.hartwig.hmftools.common.driver.panel.DriverGeneGermlineReporting.ANY;
import static com.hartwig.hmftools.common.driver.panel.DriverGeneGermlineReporting.NONE;
import static com.hartwig.hmftools.common.purple.PurpleCommon.DEFAULT_DRIVER_AMPLIFICATION_PLOIDY_RATIO;
import static com.hartwig.hmftools.common.purple.PurpleCommon.DEFAULT_DRIVER_HET_DELETION_THRESHOLD;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.createFieldsIndexMap;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.driver.DriverCategory;

public final class DriverGeneFile
{
    private static final String OTHER_TRANS_DELIM = ";";

    private DriverGeneFile() {}

    public static void write(final String filename, final List<DriverGene> driverGenes) throws IOException
    {
        List<DriverGene> sorted = Lists.newArrayList(driverGenes);
        Files.write(new File(filename).toPath(), toLines(sorted));
    }

    public static List<DriverGene> read(final String filename) throws IOException
    {
        return fromLines(Files.readAllLines(new File(filename).toPath()));
    }

    private enum Columns
    {
        gene,
        reportMissense,
        reportNonsense,
        reportSplice,
        reportDeletion,
        reportHetDeletion,
        reportLoh,
        hetDeletionThreshold,
        reportDisruption,
        reportAmplification,
        amplificationRatio,
        reportSomaticHotspot,
        likelihoodType,
        reportGermlineVariant,
        reportGermlineHotspot,
        reportGermlineDisruption,
        reportGermlineDeletion,
        reportGermlineAmplification,
        additionalReportedTranscripts,
        reportPGX;
    }

    public static List<DriverGene> fromLines(final List<String> lines)
    {
        List<DriverGene> driverGenes = Lists.newArrayList();

        String header = lines.get(0);
        Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(header, TSV_DELIM);
        lines.remove(0);

        int geneIndex = fieldsIndexMap.get(Columns.gene.toString());
        int missenseIndex = fieldsIndexMap.get(Columns.reportMissense.toString());
        int nonsenseIndex = fieldsIndexMap.get(Columns.reportNonsense.toString());
        int spliceIndex = fieldsIndexMap.get(Columns.reportSplice.toString());
        int deletionIndex = fieldsIndexMap.get(Columns.reportDeletion.toString());
        Integer hetDelIndex = fieldsIndexMap.get(Columns.reportHetDeletion.toString());
        Integer lohIndex = fieldsIndexMap.get(Columns.reportLoh.toString());
        Integer hetDelThresholdIndex = fieldsIndexMap.get(Columns.hetDeletionThreshold.toString());
        int disruptionIndex = fieldsIndexMap.get(Columns.reportDisruption.toString());
        int amplificationIndex = fieldsIndexMap.get(Columns.reportAmplification.toString());
        Integer ampRatioIndex = fieldsIndexMap.get(Columns.amplificationRatio.toString());

        int somaticHotspotIndex = fieldsIndexMap.get(Columns.reportSomaticHotspot.toString());

        int likelihoodTypeIndex = fieldsIndexMap.get(Columns.likelihoodType.toString());
        int germlineVariantIndex = fieldsIndexMap.get(Columns.reportGermlineVariant.toString());
        int germlineHotspotIndex = fieldsIndexMap.get(Columns.reportGermlineHotspot.toString());
        int germlineDisruptionIndex = fieldsIndexMap.get(Columns.reportGermlineDisruption.toString());
        int germlineDeletionIndex = fieldsIndexMap.get(Columns.reportGermlineDeletion.toString());
        Integer germlineAmpIndex = fieldsIndexMap.get(Columns.reportGermlineAmplification.toString());
        int altTransIndex = fieldsIndexMap.get(Columns.additionalReportedTranscripts.toString());
        int reportPGXIndex = fieldsIndexMap.get(Columns.reportPGX.toString());

        ImmutableDriverGene.Builder builder = ImmutableDriverGene.builder();

        for(String line : lines)
        {
            String[] values = line.split(TSV_DELIM, -1);

            List<String> otherReportableTrans = Arrays.stream(values[altTransIndex].split(OTHER_TRANS_DELIM))
                    .filter(x -> !x.isEmpty())
                    .collect(Collectors.toList());

            // backwards compatibility prior to pipeline v5.32
            String reportGermlineDisruptionStr = values[germlineDisruptionIndex];
            DriverGeneGermlineReporting reportGermlineDisruption;

            if(reportGermlineDisruptionStr.toLowerCase().equals(Boolean.TRUE.toString()))
                reportGermlineDisruption = ANY;
            else if(reportGermlineDisruptionStr.toLowerCase().equals(Boolean.FALSE.toString()))
                reportGermlineDisruption = NONE;
            else
                reportGermlineDisruption = DriverGeneGermlineReporting.valueOf(reportGermlineDisruptionStr);

            DriverGeneGermlineReporting reportGermlineDeletion = DriverGeneGermlineReporting.valueOf(values[germlineDeletionIndex]);

            builder.gene(values[geneIndex])
                    .reportMissenseAndInframe(Boolean.parseBoolean(values[missenseIndex]))
                    .reportNonsenseAndFrameshift(Boolean.parseBoolean(values[nonsenseIndex]))
                    .reportSplice(Boolean.parseBoolean(values[spliceIndex]))
                    .reportDeletion(Boolean.parseBoolean(values[deletionIndex]))
                    .reportHetDeletion(hetDelIndex != null ? Boolean.parseBoolean(values[hetDelIndex]) : false)
                    .reportLoh(lohIndex != null ? Boolean.parseBoolean(values[lohIndex]) : false)
                    .hetDeletionThreshold(hetDelThresholdIndex != null ?
                            Double.parseDouble(values[hetDelThresholdIndex]) : DEFAULT_DRIVER_HET_DELETION_THRESHOLD)
                    .reportDisruption(Boolean.parseBoolean(values[disruptionIndex]))
                    .reportAmplification(Boolean.parseBoolean(values[amplificationIndex]))
                    .amplificationRatio(ampRatioIndex != null ?
                            Double.parseDouble(values[ampRatioIndex]) : DEFAULT_DRIVER_AMPLIFICATION_PLOIDY_RATIO)
                    .reportSomaticHotspot(Boolean.parseBoolean(values[somaticHotspotIndex]))
                    .likelihoodType(DriverCategory.valueOf(values[likelihoodTypeIndex]))
                    .reportGermlineVariant(DriverGeneGermlineReporting.valueOf(values[germlineVariantIndex].toUpperCase()))
                    .reportGermlineHotspot(DriverGeneGermlineReporting.valueOf(values[germlineHotspotIndex].toUpperCase()))
                    .reportGermlineDisruption(reportGermlineDisruption)
                    .reportGermlineDeletion(reportGermlineDeletion)
                    .reportGermlineAmplification(germlineAmpIndex != null ? Boolean.parseBoolean(values[germlineAmpIndex]) : false)
                    .additionalReportedTranscripts(otherReportableTrans)
                    .reportPGX(Boolean.parseBoolean(values[reportPGXIndex]));

            driverGenes.add(builder.build());
        }

        return driverGenes;
    }

    private static String otherReportableTransStr(final List<String> otherTrans)
    {
        if(otherTrans.isEmpty())
            return "";

        StringJoiner sj = new StringJoiner(OTHER_TRANS_DELIM);
        otherTrans.forEach(x -> sj.add(x));
        return sj.toString();
    }

    private static String header()
    {
        StringJoiner sj = new StringJoiner(TSV_DELIM);

        for(Columns column : Columns.values())
        {
            sj.add(column.toString());
        }

        return sj.toString();
    }

    private static String toString(final DriverGene gene)
    {
        return new StringJoiner(TSV_DELIM)
                .add(gene.gene())
                .add(String.valueOf(gene.reportMissenseAndInframe()))
                .add(String.valueOf(gene.reportNonsenseAndFrameshift()))
                .add(String.valueOf(gene.reportSplice()))
                .add(String.valueOf(gene.reportDeletion()))
                .add(String.valueOf(gene.reportHetDeletion()))
                .add(String.valueOf(gene.reportLoh()))
                .add(String.valueOf(gene.hetDeletionThreshold()))
                .add(String.valueOf(gene.reportDisruption()))
                .add(String.valueOf(gene.reportAmplification()))
                .add(String.valueOf(gene.amplificationRatio()))
                .add(String.valueOf(gene.reportSomaticHotspot()))
                .add(String.valueOf(gene.likelihoodType()))
                .add(String.valueOf(gene.reportGermlineVariant()))
                .add(String.valueOf(gene.reportGermlineHotspot()))
                .add(String.valueOf(gene.reportGermlineDisruption()))
                .add(String.valueOf(gene.reportGermlineDeletion()))
                .add(String.valueOf(gene.reportGermlineAmplification()))
                .add(otherReportableTransStr(gene.additionalReportedTranscripts()))
                .add(String.valueOf(gene.reportPGX()))
                .toString();
    }

    private static List<String> toLines(final List<DriverGene> driverGenes)
    {
        final List<String> lines = Lists.newArrayList();
        lines.add(header());
        driverGenes.stream().map(DriverGeneFile::toString).forEach(lines::add);
        return lines;
    }
}
