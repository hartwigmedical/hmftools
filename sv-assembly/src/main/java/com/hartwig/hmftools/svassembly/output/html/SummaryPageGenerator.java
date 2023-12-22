package com.hartwig.hmftools.svassembly.output.html;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.genome.position.GenomePosition;
import com.hartwig.hmftools.common.genome.position.GenomePositions;
import com.hartwig.hmftools.svassembly.assembly.AssemblyExtenderCounters;
import com.hartwig.hmftools.svassembly.assembly.JunctionMetrics;
import com.hartwig.hmftools.svassembly.assembly.PrimaryAssemblerCounters;
import com.hartwig.hmftools.svassembly.processor.OverallCounters;
import com.hartwig.hmftools.svassembly.processor.VariantCall;
import com.hartwig.hmftools.svassembly.util.Counter;
import com.hartwig.hmftools.svassembly.util.StringUtils;

import org.apache.commons.compress.utils.IOUtils;
import org.apache.commons.lang3.tuple.Pair;
import org.jetbrains.annotations.Nullable;

public enum SummaryPageGenerator
{
    ;

    private static void appendResource(final StringBuilder page, final String resourceName)
    {
        try
        {
            final String resource = new String(IOUtils.toByteArray(Objects.requireNonNull(
                    SummaryPageGenerator.class.getResourceAsStream(resourceName))));
            page.append(resource);
        }
        catch(final Exception e)
        {
            throw new RuntimeException(e);
        }
    }

    private static String truncate(@Nullable final String s, final int maxLength)
    {
        if (s == null)
            return "";

        return s.length() < maxLength
                ? s
                : s.substring(0, maxLength - 4) + " ...";
    }

    public static void generatePage(final String folder, final OverallCounters counters, final List<VariantCall> results)
    {
        final List<Pair<GenomePosition, VariantCall>> sortedResults = new ArrayList<>();
        for(final VariantCall call : results)
        {
            if(call.LeftChromosome != null)
                sortedResults.add(Pair.of(GenomePositions.create(call.LeftChromosome, call.LeftPosition), call));
            else
                sortedResults.add(Pair.of(GenomePositions.create(call.RightChromosome, call.RightPosition), call));
        }
        sortedResults.sort(Map.Entry.comparingByKey());

        //noinspection ResultOfMethodCallIgnored
        new File(folder).mkdirs();

        final HTMLBuilder builder = new HTMLBuilder();
        builder.append("<!DOCTYPE html>\n");
        builder.appendStartTag("<html>");
        builder.appendStartTag("<head>");
        builder.append("<title>Results Summary</title>\n");
        builder.appendStartTag("<style>");
        appendResource(builder.getStringBuilder(), "/summary.css");
        builder.appendEndTag("</style>");
        builder.appendEndTag("</head>");

        builder.appendStartTag("<body>");

        builder.appendStartTag("<div class=\"summaryInfo\">");

        builder.appendStartTag("<div class=\"summaryBox\">");
        builder.appendCountersTable(counters.all(), 1);
        builder.appendEndTag("</div>");
        builder.appendStartTag("<div class=\"summaryBox\">");
        builder.appendCountersTable(counters.PrimaryAssemblerCounters.all(), 1);
        builder.appendEndTag("</div>");
        builder.appendStartTag("<div class=\"summaryBox\">");
        builder.appendCountersTable(counters.AssemblyExtenderCounters.all(), 1);
        builder.appendEndTag("</div>");

        builder.appendEndTag("</div>");

        appendPieCharts(counters, builder);

        builder.appendStartTag("<div class=\"junctionList\">");
        builder.appendStartTag("<table class=\"smart-table\">");
        builder.appendStartTag("<thead>");
        builder.append("<tr>")
                .append("<th>Left Position</th>")
                .append("<th>Right Position</th>")
                .append("<th>Assemblies</th>")
                .append("<th>Support</th>")
                .append("<th>Qual</th>")
                .append("<th>Classification</th>")
                .append("<th>Length</th>")
                .append("<th>Assembly Length</th>")
                .append("<th>Process Time</th>")
                .append("<th>Left Variant</th>")
                .append("<th>Right Variant</th>")
                .append("</tr>\n");
        builder.appendEndTag("</thead>");
        builder.appendStartTag("<tbody>");
        for(final Pair<GenomePosition, VariantCall> pair : sortedResults)
        {
            final VariantCall call = pair.getValue();

            builder.append("<tr>");
            final String junctionURL = String.format("./Variant_%s.html", call.name());

            final String junctionLeftLabel = call.LeftChromosome == null ? "" : call.LeftChromosome + ":" + call.LeftPosition;
            final String junctionLeftLink = String.format("<a href=\"%s\">%s</a>", junctionURL, junctionLeftLabel);
            builder.append("<td>").append(junctionLeftLink).append("</td>");

            final String junctionRightLabel = call.RightChromosome == null ? "" : call.RightChromosome + ":" + call.RightPosition;
            final String junctionRightLink = String.format("<a href=\"%s\">%s</a>", junctionURL, junctionRightLabel);
            builder.append("<td>").append(junctionRightLink).append("</td>");

            builder.append("<td>").append(call.associatedAssemblies().size()).append("</td>");
            builder.append("<td>").append("%s/%s", call.somaticSupport(), call.germlineSupport()).append("</td>");
            builder.append("<td>").append(call.quality()).append("</td>");

            builder.append("<td>").append(call.Classification.Type).append("</td>");
            builder.append("<td>").append(call.Classification.Length).append("</td>");

            final String assemblyLengths = call.associatedAssemblies().stream()
                    .mapToInt(assembly -> assembly.Assembly.length())
                    .distinct()
                    .sorted()
                    .mapToObj(String::valueOf)
                    .collect(Collectors.joining(" | "));
            builder.append("<td>").append(assemblyLengths).append("</td>");

            final long processTimeNanos = call.associatedAssemblies().stream()
                    .mapToLong(assembly ->
                    {
                        final PrimaryAssemblerCounters primaryCounters = assembly.getAllErrata(JunctionMetrics.class).stream()
                                .map(metrics -> metrics.Counters)
                                .reduce(new PrimaryAssemblerCounters(), (l, r) -> (PrimaryAssemblerCounters) l.add(r));

                        final AssemblyExtenderCounters extension = assembly.getAllErrata(AssemblyExtenderCounters.class).stream()
                                .reduce(new AssemblyExtenderCounters(), (l, r) -> (AssemblyExtenderCounters) l.add(r));

                        return primaryCounters.ProcessTimeNanos.getValue() + extension.ProcessTimeNanos.getValue();
                    }).sum();
            builder.append("<td>").append(StringUtils.formatNanos(processTimeNanos)).append("</td>");

            builder.append("<td><code>").append(truncate(call.LeftDescriptor, 80)).append("</code></td>");
            builder.append("<td><code>").append(truncate(call.RightDescriptor, 80)).append("</code></td>");
            builder.append("</tr>\n");
        }
        builder.appendEndTag("</tbody>");
        builder.appendEndTag("</table>");
        builder.appendEndTag("</div>");

        builder.appendEndTag("</body>");
        builder.appendEndTag("</html>");

        try
        {
            Files.writeString(Path.of(folder, "Summary.html"), builder.toString());
        }
        catch(final IOException e)
        {
            throw new RuntimeException(e);
        }
    }

    private static void appendPieCharts(final OverallCounters overall, final HTMLBuilder builder)
    {
        builder.appendStartTag("<div class=\"performanceSummary\">");

        final PieChart mainPieChart = new PieChart("Overall");
        append(mainPieChart, overall.PrimaryAssemblyTime);
        append(mainPieChart, overall.InterJunctionDeduplicationTime);
        append(mainPieChart, overall.ExtensionTime);
        append(mainPieChart, overall.PrimaryPhasingTime);
        append(mainPieChart, overall.PhasedAssemblyMergingTime);
        append(mainPieChart, overall.SecondaryPhasingTime);
        append(mainPieChart, overall.MergeSecondaryTime);
        append(mainPieChart, overall.AlignmentTime);
        append(mainPieChart, overall.HomologyTime);
        append(mainPieChart, overall.VariantCallingTime);
        append(mainPieChart, overall.VariantDeduplicationCounters.OverallTimeNanos);
        append(mainPieChart, overall.SupportScanTime);
        mainPieChart.appendAsSVG(builder, "width: 400px;");

        final PieChart ioPieChart = new PieChart("IO Breakdown");
        append(ioPieChart, overall.PrimaryAssemblerCounters.InitialReadTimeNanos);
        append(ioPieChart, overall.AssemblyExtenderCounters.SubsequentReadTimeNanos);
        append(ioPieChart, overall.AssemblyExtenderCounters.DiscordantSearchTimeNanos);
        final long processingTimeNanos = overall.PrimaryAssemblerCounters.ProcessTimeNanos.getValue()
                + overall.AssemblyExtenderCounters.ProcessTimeNanos.getValue()
                - overall.PrimaryAssemblerCounters.InitialReadTimeNanos.getValue()
                - overall.AssemblyExtenderCounters.SubsequentReadTimeNanos.getValue()
                - overall.AssemblyExtenderCounters.DiscordantSearchTimeNanos.getValue();
        ioPieChart.add("Processing Time", processingTimeNanos, StringUtils.formatNanos(processingTimeNanos));
        append(ioPieChart, overall.AlignmentTime);
        ioPieChart.appendAsSVG(builder, "width: 400px;");

        final PieChart primaryPieChart = new PieChart("Primary Assembly");
        append(primaryPieChart, overall.PrimaryAssemblerCounters.InitialReadTimeNanos);
        append(primaryPieChart, overall.PrimaryAssemblerCounters.DiagramGenerationTimeNanos);
        append(primaryPieChart, overall.PrimaryAssemblerCounters.GraphSimplificationTimeNanos);
        append(primaryPieChart, overall.PrimaryAssemblerCounters.JunctionConstructionTimeNanos);
        append(primaryPieChart, overall.PrimaryAssemblerCounters.JunctionExtensionTimeNanos);
        append(primaryPieChart, overall.PrimaryAssemblerCounters.AnchorConstructionTimeNanos);
        primaryPieChart.appendAsSVG(builder, "width: 400px;");

        final PieChart extensionPieChart = new PieChart("Assembly Extension");
        append(extensionPieChart, overall.AssemblyExtenderCounters.SubsequentReadTimeNanos);
        append(extensionPieChart, overall.AssemblyExtenderCounters.DiscordantSearchTimeNanos);
        append(extensionPieChart, overall.AssemblyExtenderCounters.DiagramGenerationTimeNanos);
        append(extensionPieChart, overall.AssemblyExtenderCounters.GraphSimplificationTimeNanos);
        append(extensionPieChart, overall.AssemblyExtenderCounters.ExtendLeftTimeNanos);
        append(extensionPieChart, overall.AssemblyExtenderCounters.ExtendRightTimeNanos);
        append(extensionPieChart, overall.AssemblyExtenderCounters.CleanupTimeNanos);
        extensionPieChart.appendAsSVG(builder, "width: 400px;");

        final PieChart deduplicationPieChart = new PieChart("Deduplication");
        append(deduplicationPieChart, overall.VariantDeduplicationCounters.MatchedGroupWallTimeNanos);
        append(deduplicationPieChart, overall.VariantDeduplicationCounters.AssemblyDedupeWallTimeNanos);
        append(deduplicationPieChart, overall.VariantDeduplicationCounters.SinglesDedupeWallTimeNanos);
        append(deduplicationPieChart, overall.VariantDeduplicationCounters.NearbyDedupeWallTimeNanos);
        deduplicationPieChart.appendAsSVG(builder, "width: 400px;");

        builder.appendEndTag("</div>");
    }

    private static void append(final PieChart chart, final Counter counter)
    {
        chart.add(counter.Name, counter.getValue(), counter.formatValue());
    }
}
