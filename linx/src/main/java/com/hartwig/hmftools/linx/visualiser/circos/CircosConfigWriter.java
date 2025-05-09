package com.hartwig.hmftools.linx.visualiser.circos;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.nio.charset.Charset;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.util.StringJoiner;
import java.util.function.Function;

import com.hartwig.hmftools.linx.visualiser.CircosConfig;

import org.apache.logging.log4j.core.util.IOUtils;
import org.jetbrains.annotations.NotNull;

public class CircosConfigWriter
{
    private final String sample;
    private final String outputDir;
    private final CircosData circosData;
    private final CircosConfig config;

    private final double exonOuterRadius;
    private final double exonInnerRadius;
    private final double geneOuterRadius;
    private final double geneInnerRadius;

    private final double segmentOuterRadius;
    private final double segmentInnerRadius;

    private final double copyNumberOuterRadius;
    private final double copyNumberMiddleRadius;
    private final double copyNumberInnerRadius;

    private final double cobaltOuterRadius;
    private final double cobaltInnerRadius;

    private final double mapOuterRadius;
    private final double mapMiddleRadius;
    private final double mapInnerRadius;

    private final double amberOuterRadius;
    private final double amberInnerRadius;

    private final double labelSize;
    private final double labelPosition;

    public CircosConfigWriter(final String filename, final String outputDir, final CircosData data, final CircosConfig config)
    {
        this.sample = filename;
        this.circosData = data;
        this.config = config;
        this.outputDir = outputDir;
        this.labelSize = data.labelSize();

        double gapSize = config.GapRadius;
        boolean displayGenes = data.displayGenes();

        double geneRelativeSize = displayGenes ? config.GeneRelativeSize : 0;
        double segmentRelativeSize = config.SegmentRelativeSize;
        double copyNumberRelativeSize = config.CopyNumberRelativeSize;

        double totalRelativeSize = geneRelativeSize + segmentRelativeSize + copyNumberRelativeSize;

        int numberOfGaps = displayGenes ? 5 : 3;
        if(!data.AmberBAFs.isEmpty()) numberOfGaps += 1;
        if(!data.CobaltRatios.isEmpty()) numberOfGaps += 1;

        double totalSpaceAvailable = 1 - numberOfGaps * gapSize - config.InnerRadius;
        double purpleSpaceAvailable = copyNumberRelativeSize / totalRelativeSize * totalSpaceAvailable;
        int cnaGainTracks = Math.max(2, (int) Math.round(Math.ceil(data.MaxCopyNumber - 2)));
        int mapGainTracks = Math.max(1, (int) Math.round(Math.ceil(data.MaxMinorAllelePloidy - 1)));
        int amberTracks = !data.AmberBAFs.isEmpty() ? 1 : 0;
        int cobaltTracks = !data.CobaltRatios.isEmpty() ? 1 : 0;
        double purpleTrackSize = purpleSpaceAvailable / (1 + 2 + cnaGainTracks + mapGainTracks + amberTracks + cobaltTracks);

        if(displayGenes)
        {
            exonOuterRadius = 1 - gapSize - config.ExonRankRadius;
            exonInnerRadius = exonOuterRadius - geneRelativeSize / totalRelativeSize * totalSpaceAvailable;
            segmentOuterRadius = exonInnerRadius - gapSize;
        }
        else
        {
            exonOuterRadius = 0;
            exonInnerRadius = 0;
            segmentOuterRadius = 1 - gapSize;
        }

        // Radius positions for each track are calculated starting from the outermost track
        double exonDistance = exonOuterRadius - exonInnerRadius;
        geneOuterRadius = exonOuterRadius - 9d / 20d * exonDistance;
        geneInnerRadius = exonInnerRadius + 9d / 20d * exonDistance;

        segmentInnerRadius = segmentOuterRadius - segmentRelativeSize / totalRelativeSize * totalSpaceAvailable;

        copyNumberOuterRadius = segmentInnerRadius - gapSize;
        copyNumberMiddleRadius = copyNumberOuterRadius - cnaGainTracks * purpleTrackSize;
        copyNumberInnerRadius = copyNumberMiddleRadius - 2 * purpleTrackSize;

        cobaltOuterRadius = copyNumberInnerRadius - gapSize;
        cobaltInnerRadius = cobaltOuterRadius - cobaltTracks * 2 * purpleTrackSize;

        mapOuterRadius = cobaltInnerRadius - gapSize;
        mapMiddleRadius = mapOuterRadius - mapGainTracks * purpleTrackSize;
        mapInnerRadius = mapMiddleRadius - 1 * purpleTrackSize;

        amberOuterRadius = mapInnerRadius - gapSize;
        amberInnerRadius = amberOuterRadius - amberTracks * 2 * purpleTrackSize;

        labelPosition = copyNumberMiddleRadius + purpleTrackSize*0.05; // Add slight offset so that the label isn't flush with the axis
    }

    public double svTrackRelative(int track)
    {
        double difference = segmentOuterRadius() - segmentInnerRadius;
        double singleTrack = difference / circosData.MaxTracks;

        return segmentInnerRadius + track * singleTrack;
    }

    private double segmentOuterRadius()
    {
        return segmentOuterRadius;
    }

    public String writeConfig(int frame) throws IOException
    {
        final String fileName = sample + ".circos." + String.format("%03d", frame) + ".conf";
        final String configPath = outputDir + File.separator + fileName;

        int chromosomeCount = circosData.ContigLengths.size();
        int totalContigLength = circosData.totalContigLength();

        int cnaMaxTracks = Math.max(2, (int) Math.round(Math.ceil(circosData.MaxCopyNumber - 2)));

        int mapMaxTracks = Math.max(1, (int) Math.round(Math.ceil(circosData.MaxMinorAllelePloidy - 1)));

        final Charset charset = StandardCharsets.UTF_8;
        final String template =
                readResource("/visualisation/cluster.template")
                        .replaceAll("SUBSTITUTE_CURRENT_FRAME", String.valueOf(frame))
                        .replaceAll("SUBSTITUTE_IDEOGRAM_RADIUS", String.valueOf(config.OuterRadius))
                        .replaceAll("SUBSTITUTE_IDEOGRAM_SPACING", chromosomeCount > 1 ? "0.005r" : 0.005 * totalContigLength + "u")

                        .replaceAll("SUBSTITUTE_EXON_INNER_RADIUS", String.valueOf(exonInnerRadius))
                        .replaceAll("SUBSTITUTE_EXON_OUTER_RADIUS", String.valueOf(exonOuterRadius))
                        .replaceAll("SUBSTITUTE_GENE_INNER_RADIUS", String.valueOf(geneInnerRadius))
                        .replaceAll("SUBSTITUTE_GENE_OUTER_RADIUS", String.valueOf(geneOuterRadius))

                        .replaceAll("SUBSTITUTE_SV_INNER_RADIUS", String.valueOf(segmentInnerRadius))
                        .replaceAll("SUBSTITUTE_SV_OUTER_RADIUS", String.valueOf(segmentOuterRadius))

                        .replaceAll("SUBSTITUTE_MAP_INNER_RADIUS", String.valueOf(mapInnerRadius))
                        .replaceAll("SUBSTITUTE_MAP_OUTER_RADIUS", String.valueOf(mapOuterRadius))
                        .replaceAll("SUBSTITUTE_MAP_MIDDLE_RADIUS", String.valueOf(mapMiddleRadius))
                        .replaceAll("SUBSTITUTE_MAP_GAIN_MAX", String.valueOf(mapMaxTracks))
                        .replaceAll("SUBSTITUTE_MAP_GAIN_SPACING", String.valueOf(1d / mapMaxTracks))

                        .replaceAll("SUBSTITUTE_COBALT_INNER_RADIUS", String.valueOf(cobaltInnerRadius))
                        .replaceAll("SUBSTITUTE_COBALT_OUTER_RADIUS", String.valueOf(cobaltOuterRadius))

                        .replaceAll("SUBSTITUTE_CNA_INNER_RADIUS", String.valueOf(copyNumberInnerRadius))
                        .replaceAll("SUBSTITUTE_CNA_OUTER_RADIUS", String.valueOf(copyNumberOuterRadius))
                        .replaceAll("SUBSTITUTE_CNA_MIDDLE_RADIUS", String.valueOf(copyNumberMiddleRadius))
                        .replaceAll("SUBSTITUTE_CNA_GAIN_MAX", String.valueOf(cnaMaxTracks))
                        .replaceAll("SUBSTITUTE_CNA_GAIN_AXIS_POSITION", cnaAxisPositions(cnaMaxTracks))
                        .replaceAll("SUBSTITUTE_CNA_DISTANCE_RADIUS", labelPosition + "r")

                        .replaceAll("SUBSTITUTE_AMBER_INNER_RADIUS", String.valueOf(amberInnerRadius))
                        .replaceAll("SUBSTITUTE_AMBER_OUTER_RADIUS", String.valueOf(amberOuterRadius))

                        .replaceAll("SUBSTITUTE_SV_SPACING", String.valueOf(1d / circosData.MaxTracks))
                        .replaceAll("SUBSTITUTE_SV_MAX", String.valueOf(circosData.MaxTracks))

                        .replaceAll("SUBSTITUTE_LABEL_SIZE", String.valueOf(labelSize))

                        .replaceAll("SUBSTITUTE_SAMPLE", sample);

        Files.write(new File(configPath).toPath(), template.getBytes(charset));
        return fileName;
    }

    private String readResource(final String resource) throws IOException
    {
        InputStream in = getClass().getResourceAsStream(resource);
        BufferedReader reader = new BufferedReader(new InputStreamReader(in));
        return IOUtils.toString(reader);
    }

    static String cnaAxisPositions(int maxTracks)
    {
        StringJoiner builder = new StringJoiner(",");

        final double rel = 1d / maxTracks;
        final Function<Integer, String> relString = i -> (Math.round(i * rel * 10000) / 10000D) + "r";

        for(int i = 1; i <= Math.min(7, maxTracks); i++)
        {
            builder.add(relString.apply(i));
        }

        for(int i = 8; i <= maxTracks; i += 10)
        {
            builder.add(relString.apply(i));
        }

        return builder.toString();
    }

}
