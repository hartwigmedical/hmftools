package com.hartwig.hmftools.esvee.html;

import java.util.AbstractMap;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashSet;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.Set;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.esvee.models.AlignedAssembly;
import com.hartwig.hmftools.esvee.models.AlignedSequence;
import com.hartwig.hmftools.esvee.models.Alignment;
import com.hartwig.hmftools.esvee.models.ExtendedAssembly;
import com.hartwig.hmftools.esvee.models.PrimaryAssembly;
import com.hartwig.hmftools.esvee.read.Read;
import com.hartwig.hmftools.esvee.models.Sequence;
import com.hartwig.hmftools.esvee.models.SimpleSequence;
import com.hartwig.hmftools.esvee.models.SupportedAssembly;
import com.hartwig.hmftools.esvee.util.Counter;
import com.hartwig.hmftools.esvee.assembly.AssemblyExtenderCounters;
import com.hartwig.hmftools.esvee.assembly.JunctionMetrics;
import com.hartwig.hmftools.esvee.assembly.PrimaryAssemblerCounters;
import com.hartwig.hmftools.esvee.assembly.SupportChecker;
import com.hartwig.hmftools.esvee.mermaid.Flowchart;
import com.hartwig.hmftools.esvee.common.VariantCall;

import org.apache.commons.lang3.tuple.Pair;
import org.jetbrains.annotations.Nullable;

public class AssemblyView
{
    private final RefGenomeInterface mRefGenomeSource;
    private final SupportChecker mSupportChecker;

    public AssemblyView(final RefGenomeInterface refGenomeSource, final SupportChecker supportChecker)
    {
        mRefGenomeSource = refGenomeSource;
        mSupportChecker = supportChecker;
    }

    public void generate(final HTMLBuilder builder, final VariantCall call, final AlignedAssembly assembly, final boolean basicInfoOnly)
    {
        builder.append("<h2>Assembly %s</h2>\n", assembly.Name);
        builder.appendStartTag("<div class=\"assembly-summary collapsible expanded\">");
        builder.appendStartTag("<div class=\"assembly-summary-inner\">");

        builder.append("<h3>Counters</h3>\n");
        builder.appendStartTag("<div class=\"assembly-counters\">");
        final PrimaryAssemblerCounters primaryCounters = assembly.getAllErrata(JunctionMetrics.class).stream()
                .map(metrics -> metrics.Counters)
                .reduce(new PrimaryAssemblerCounters(), (l, r) -> (PrimaryAssemblerCounters) l.add(r));

        final AssemblyExtenderCounters extension = assembly.getAllErrata(AssemblyExtenderCounters.class).stream()
                .reduce(new AssemblyExtenderCounters(), (l, r) -> (AssemblyExtenderCounters) l.add(r));
        final List<Counter> counters = new ArrayList<>();
        counters.addAll(primaryCounters.all());
        counters.addAll(extension.all());

        builder.appendCountersTable(counters, 3);

        builder.appendEndTag("</div>");

        builder.append("<h3>Alignment</h3>\n");
        builder.appendStartTag("<div class=\"assembly-alignment\">");
        generateAlignmentTable(builder, assembly);
        builder.appendEndTag("</div>");

        builder.append("<h3>Assembly & Support</h3>\n");
        builder.append("<p>Supporting Fragments: ").append(assembly.getSupportFragments().size()).append("</p>\n");

        generateSupportTable(builder, assembly);

        final SequenceView sequenceView = new SequenceView(mRefGenomeSource);
        sequenceView.addHeaderSequence(assembly);
        final Set<SupportedAssembly> supportingAssemblies = new LinkedHashSet<>();
        for(ExtendedAssembly source : assembly.Source.Sources)
        {
            @Nullable
            SupportedAssembly innerSource = source.Source;
            while(innerSource != null)
            {
                supportingAssemblies.add(innerSource);

                if(innerSource instanceof ExtendedAssembly)
                    innerSource = ((ExtendedAssembly) innerSource).Source;
                else
                    break;
            }
        }
        supportingAssemblies.addAll(assembly.getAllErrata(PrimaryAssembly.class));
        supportingAssemblies.addAll(assembly.getAllErrata(ExtendedAssembly.class));
        for(SupportedAssembly support : supportingAssemblies)
        {
            @Nullable
            final Integer index = determineBestSupportIndex(assembly, support);
            if(index != null)
                sequenceView.addSequence(support, index, "lost");
        }

        sequenceView.addPointOfInterest(call.LeftChromosome, call.LeftPosition);
        sequenceView.addPointOfInterest(call.RightChromosome, call.RightPosition);

        if(!basicInfoOnly)
        {
            // Sort fragments with their partner, ordered by the lowest start index.
            final Set<Read> finalSupport = new HashSet<>();
            final Set<Read> lostSupport = new HashSet<>();
            assembly.getSupportRecords().forEach(finalSupport::add);
            for(ExtendedAssembly source : assembly.Source.Sources)
            {
                source.getSupportRecords().forEach(support ->
                {
                    if(!finalSupport.contains(support))
                        lostSupport.add(support);
                });
                source.Source.getSupportRecords().forEach(support ->
                {
                    if(!finalSupport.contains(support))
                        lostSupport.add(support);
                });
            }

            final List<Map.Entry<Read, Integer>> rawSupport = assembly.getSupport().toList();
            for(Read read : lostSupport)
            {
                final int supportIndex = Objects.requireNonNullElse(determineBestSupportIndex(assembly, read), 0);
                rawSupport.add(new AbstractMap.SimpleEntry<>(read, supportIndex));
            }

            final Map<String, Integer> fragmentsByStartIndex = rawSupport.stream()
                    .collect(Collectors.toMap(kvp -> kvp.getKey().getName(), Map.Entry::getValue, Math::min));

            final List<Map.Entry<Read, Integer>> sortedSupport = rawSupport.stream()
                    .sorted(Comparator.<Map.Entry<Read, Integer>, Integer>comparing(kvp -> fragmentsByStartIndex.get(kvp.getKey()
                                    .getName()))
                            .thenComparing(kvp -> kvp.getKey().getName())
                            .thenComparingInt(Map.Entry::getValue))
                    .collect(Collectors.toList());

            for(int i = 0; i < sortedSupport.size(); i++)
            {
                final Map.Entry<Read, Integer> entry = sortedSupport.get(i);
                final Read read = entry.getKey();
                final int supportIndex = entry.getValue();

                final Map.Entry<Read, Integer> nextEntry = i + 1 < sortedSupport.size() ? sortedSupport.get(i + 1) : entry;
                final int thisEnd = supportIndex + read.getLength();
                final int nextBegin = nextEntry.getValue();
                final int gapSize = nextBegin - thisEnd;
                final boolean canCombine = read.getName().equals(nextEntry.getKey().getName()) && gapSize > 3;

                final boolean isReference = false; // TO-DO or just show by sampleId as per Sage  read.isGermline();
                final boolean isLost = !canCombine && lostSupport.contains(read);

                final Sequence sequence;
                if(canCombine)
                {
                    final String newName = read.getName() + (read.secondInPair() ? " 2&1" : " 1&2");
                    // FIXME: see comment below
                    // sequence = combine(newName, read, nextEntry.getKey(), gapSize);
                    i++; // Skip next
                }
                else
                    sequence = read;

                final String cssClass;
                if(isLost && isReference)
                    cssClass = "lost fromRef";
                else if(isLost)
                    cssClass = "lost";
                else if(isReference)
                    cssClass = "fromRef";
                else
                    cssClass = "default";

                // FIXME: requires Read to have Alignments but consider an alternative
                /*
                if(isLost)
                    sequenceView.addSequence(read, cssClass);
                else
                    sequenceView.addSequence(sequence, supportIndex, cssClass);
                 */
            }
        }

        builder.appendStartTag("<div class=\"read-summary\">");
        builder.appendStartTag("<div class=\"table-wrapper force-horizontal-scroll\">");
        sequenceView.generate(builder);
        builder.appendEndTag("</div>");
        builder.appendEndTag("</div>");

        for(ExtendedAssembly source : assembly.Source.Sources)
            generateAssemblyDiagrams(builder, source);

        builder.appendEndTag("</div>");
        builder.appendEndTag("</div>");
    }

    private static void generateSupportTable(final HTMLBuilder builder, final AlignedAssembly assembly)
    {
        final int columns = Math.max(1, Math.min(4, assembly.getSupport().size() / 10));
        final int rows = (int) Math.ceil(assembly.getSupport().size() / (double) columns);

        builder.appendStartTag("<div class=\"table-wrapper force-horizontal-scroll\">");
        builder.appendStartTag("<table>");
        builder.appendStartTag("<thead>");
        builder.append("<tr>");
        builder.append("<th>Name</th><th>Index</th>".repeat(columns));
        builder.append("</tr>");
        builder.appendEndTag("</thead>");
        builder.appendStartTag("<tbody>");

        final List<Map.Entry<Read, Integer>> entries = assembly.getSupport().toList();
        for(int i = 0; i < rows; i++)
        {
            builder.append("<tr>");
            for(int j = 0; j < columns; j++)
            {
                final int index = i * columns + j;
                if(index < entries.size())
                {
                    final Map.Entry<Read, Integer> entry = entries.get(index);
                    final Read read = entry.getKey();
                    final Integer startIndex = entry.getValue();
                    builder.append("<td>").append(read.getName()).append("</td>").append("<td>").append(startIndex).append("</td>");
                }
            }
            builder.append("</tr>");
        }
        builder.appendEndTag("</tbody>");
        builder.appendEndTag("</table>");
        builder.appendEndTag("</div>");
    }

    private static Sequence combine(final String newName, final AlignedSequence left, final AlignedSequence right, final int gap)
    {
        final byte[] combinedBases = (left.getBasesString() + "-".repeat(gap) + right.getBasesString()).getBytes();
        final byte[] combinedQuals = new byte[left.getLength() + right.getLength() + gap];
        System.arraycopy(left.getBaseQuality(), 0, combinedQuals, 0, left.getLength());
        System.arraycopy(right.getBaseQuality(), 0, combinedQuals, left.getLength() + gap, right.getLength());

        return new SimpleSequence(newName, combinedBases, combinedQuals);
    }

    @Nullable
    private Integer determineBestSupportIndex(final Sequence assembly, final Sequence record)
    {
        return mSupportChecker.determineBestSupportIndex(assembly, record, 30, 20);
    }

    private static void generateAssemblyDiagrams(final HTMLBuilder builder, final ExtendedAssembly assembly)
    {
        generateAssemblyDiagrams(builder, assembly.getDiagrams());
    }

    private static void generateAssemblyDiagrams(final HTMLBuilder builder, final List<DiagramSet> diagrams)
    {
        builder.appendStartTag("<div class=\"assemblyDiagrams\">");
        for(DiagramSet diagramSet : diagrams)
        {
            builder.appendStartTag("<div class=\"diagramSet\">");
            builder.append("<h3>").append(diagramSet.Label).append("</h3>\n");
            for(Pair<String, Flowchart> entry : diagramSet.Diagrams)
            {
                builder.appendStartTag("<div class=\"diagram\">");
                builder.append("<h4>").append(entry.getLeft()).append("</h4>\n");
                builder.appendStartTag("<div class=\"diagramInner force-horizontal-scroll\">");
                builder.append(entry.getRight().toSVG());
                builder.appendEndTag("</div>");
                builder.appendEndTag("</div>");
            }
            builder.appendEndTag("</div>");
        }
        builder.appendEndTag("</div>");
    }

    private static void generateAlignmentTable(final HTMLBuilder builder, final AlignedAssembly assembly)
    {
        if(assembly.getAlignmentBlocks() == null)
            return;

        builder.appendStartTag("<table>");
        builder.appendStartTag("<thead>");
        builder.append("<tr><th>Chromosome</th><th>Position</th><th>Length</th><th>Inverted</th><th>Qual</th></tr>\n");
        builder.appendEndTag("</thead>");
        builder.appendStartTag("<tbody>");
        for(Alignment alignment : assembly.getAlignmentBlocks())
        {
            builder.append("<tr><td>").append(alignment.Chromosome).append("</td>");

            builder.append("<td>");
            if(!alignment.isUnmapped())
                builder.append(alignment.ReferenceStartPosition).append("-")
                        .append(alignment.ReferenceStartPosition + alignment.Length - 1);
            builder.append("</td>");

            builder.append("<td>").append(alignment.Length).append("</td>");
            builder.append("<td>").append(alignment.Inverted).append("</td>");
            builder.append("<td>").append(alignment.Quality).append("</td>");
            builder.append("</tr>\n");
        }
        builder.appendEndTag("</tbody>");
        builder.appendEndTag("</table>");
    }
}
