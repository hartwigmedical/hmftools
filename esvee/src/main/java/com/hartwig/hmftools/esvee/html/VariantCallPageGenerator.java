package com.hartwig.hmftools.esvee.html;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Collection;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.Objects;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.esvee.models.AlignedAssembly;
import com.hartwig.hmftools.esvee.assembly.JunctionMetrics;
import com.hartwig.hmftools.esvee.assembly.SupportChecker;
import com.hartwig.hmftools.esvee.common.VariantAssembly;
import com.hartwig.hmftools.esvee.common.VariantCall;

import org.apache.commons.compress.utils.IOUtils;

public final class VariantCallPageGenerator
{
    public static void generatePage(
            final String folder, final RefGenomeInterface reference, final SupportChecker supportChecker, final VariantCall call)
    {
        new File(folder).mkdirs();

        final String title = call.isSingleSided()
                ? Objects.requireNonNullElse(call.LeftChromosome, call.RightChromosome) + ":" + (call.LeftPosition + call.RightPosition)
                : (call.LeftChromosome + ":" + call.LeftPosition + "&" + call.RightChromosome + ":" + call.RightPosition);

        final HTMLBuilder page = new HTMLBuilder();
        page.append("<!DOCTYPE html>\n");
        page.appendStartTag("<html>")
                .appendStartTag("<head>")
                .append("<title>%s - Summary</title>\n", title);

        page.appendStartTag("<style>");
        appendResource(page, "/summary.css");
        page.appendEndTag("</style>");
        page.appendEndTag("</head>");

        page.appendStartTag("<body>");
        generateSummarySection(page, call);
        generateAssemblySummaries(page, reference, supportChecker, call);

        page.appendStartTag("<script type=\"module\">");
        appendResource(page, "/summary.js");
        page.appendEndTag("</script>");
        page.appendEndTag("</body>").appendEndTag("</html>");

        try
        {
            final String filename = String.format("Variant_%s.html", call.name());
            Files.writeString(Path.of(folder, filename), page.toString());
        }
        catch(final IOException e)
        {
            throw new RuntimeException(e);
        }
    }

    private static void generateSummarySection(final HTMLBuilder builder, final VariantCall call)
    {
        builder.appendStartTag("<div class=\"summary\">");
        builder.append("<h1>Summary</h1>\n");
        builder.appendStartTag("<div class=\"force-horizontal-scroll\" style=\"display:flex;\">");

        builder.appendStartTag("<div>");
        builder.append("<center><h2>Details</h2></center>\n");

        Map<String,Object> attributes = mapOf(
                "Left", Objects.requireNonNullElse(call.LeftChromosome, "?") + ":" + call.LeftPosition,
                "Right", Objects.requireNonNullElse(call.RightChromosome, "?") + ":" + call.RightPosition,
                "Left Variant", call.LeftDescriptor,
                "Right Variant", call.RightDescriptor,
                "Support", call.somaticSupport() + "/" + call.germlineSupport(),
                "Phase", call.PhaseSets.size(),
                "Assemblies", call.associatedAssemblies().size(),
                "Quality", call.quality(),
                "Classification", call.Classification,
                "Overhang", call.overhang(),
                "Split Reads", call.splitReadSupport(),
                "Discordant Pairs", call.discordantSupport());

        builder.appendSimpleAttributeTable(attributes, 2);
        builder.appendEndTag("</div>");

        builder.appendEndTag("</div>"); // Flex div

        builder.append("<h2>Assemblies</h2>\n");
        builder.appendStartTag("<div class=\"summary-assemblies force-horizontal-scroll\">");
        builder.appendStartTag("<table style=\"font-family: monospace;\">");
        builder.append("<tr><th>#</th>")
                .append("<th>Length</th>")
                .append("<th>Support</th>")
                .append("<th>L/R Offset</th>")
                .append("<th>L/R Cigar</th>")
                .append("<th>L/R Overhang</th>")
                .append("<th>Source</th>")
                .append("<th>Assembly</th></tr>\n");
        int index = 0;
        for(VariantAssembly assembly : call.variantAssemblies())
            builder.append("<tr><td>").append(++index).append("</td>")
                    .append("<td>").append(assembly.Assembly.Assembly.length()).append("</td>")
                    .append("<td>").append("%s", assembly.Assembly.getSupportFragments().size()).append("</td>")
                    .append("<td>").append(assembly.LeftPosition).append("/").append(assembly.RightPosition).append("</td>")
                    .append("<td>").append(assembly.LeftAnchorCigar).append("/").append(assembly.RightAnchorCigar).append("</td>")
                    .append("<td>").append(assembly.LeftOverhang).append("/").append(assembly.RightOverhang).append("</td>")
                    .append("<td>").append(assembly.Assembly.getAllErrata(JunctionMetrics.class).stream()
                            .map(metrics -> String.format("%s:%s (%s)", metrics.JunctionChromosome, metrics.JunctionPosition, metrics.JunctionDirection))
                            .collect(Collectors.joining(", "))).append("</td>")
                    .append("<td>").append(assembly.Assembly.Assembly).append("</td>")
                    .append("</tr>\n");
        builder.appendEndTag("</table>");
        builder.appendEndTag("</div>");
        builder.appendStartTag("</div>");
    }

    private static Map<String,Object> mapOf(final Object... params)
    {
        final Map<String, Object> result = new LinkedHashMap<>();
        for(int i = 0; i < params.length; i += 2)
            result.put(params[i].toString(), params[i + 1]);
        return result;
    }

    private static void appendResource(final HTMLBuilder page, final String resourceName)
    {
        try
        {
            final String resource = new String(IOUtils.toByteArray(Objects.requireNonNull(
                    VariantCallPageGenerator.class.getResourceAsStream(resourceName))));
            page.append(resource);
        }
        catch(final Exception e)
        {
            throw new RuntimeException(e);
        }
    }

    private static String javascriptShowAllRows() {
        return "document.querySelectorAll('.read-summary tr.sequence').forEach(tr => {tr.style.display = 'table-row'; tr.style.filter = '';});";
    }

    private static String javascriptHideAllExcept(final String assemblyName, final Collection<String> fragmentNames) {
        final String set = "[" + fragmentNames.stream().map(name -> "'" + name + "'").collect(Collectors.joining(",")) + "]";
        final String highlightAssembly = "if(tr.getAttribute('name').includes('Assembly')) {"
                + " if(tr.getAttribute('name') == '" + assemblyName + "') {"
                + " tr.style.filter = ''"
                + " } else {"
                + " tr.style.filter = 'brightness(25%)'"
                + " } }";
        return "document.querySelectorAll('.read-summary tr.sequence').forEach(tr => {"
                + " if(!" + set + ".includes(tr.getAttribute('name'))"
                + " && !tr.getAttribute('name').includes('Assembly')"
                + " && !tr.getAttribute('name').includes('Reference')) { "
                + "tr.style.display = 'none'; "
                + "} else { "
                + "tr.style.display = 'table-row'; "
                + "} "
                + (assemblyName != null ? highlightAssembly : "")
                + "});";
    }

    private static void generateAssemblySummaries(
            final HTMLBuilder builder, final RefGenomeInterface reference,
            final SupportChecker supportChecker, final VariantCall call)
    {
        builder.appendStartTag("<div class=\"candidate-assemblies\">");
        builder.append("<h1>Assemblies</h1>\n");
        
        final AssemblyView assemblyView = new AssemblyView(reference, supportChecker);
        
        for(AlignedAssembly assembly : call.associatedAssemblies())
        {
            assemblyView.generate(builder, call, assembly, false);
        }
        
        builder.appendEndTag("</div>");
    }
}
