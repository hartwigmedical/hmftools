package com.hartwig.hmftools.esvee.html;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.Optional;
import java.util.Set;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.esvee.models.AlignedSequence;
import com.hartwig.hmftools.esvee.models.Alignment;
import com.hartwig.hmftools.esvee.models.BasicAlignedSequence;
import com.hartwig.hmftools.esvee.models.Sequence;
import com.hartwig.hmftools.esvee.models.SimpleSequence;

import org.apache.commons.lang3.tuple.Pair;
import org.apache.commons.lang3.tuple.Triple;
import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.util.SequenceUtil;

public class SequenceView
{
    private final RefGenomeInterface mReferenceSource;
    private final List<AlignedSequence> mHeaderSequences = new ArrayList<>();
    private final List<Pair<String, AlignedSequence>> mSequences = new ArrayList<>();
    private final List<Triple<String, Integer, Sequence>> mRawSequences = new ArrayList<>();
    private final Map<String, Set<Integer>> mPointsOfInterest = new HashMap<>();
    private int mMinOffset = 0;

    public SequenceView(final RefGenomeInterface referenceSource)
    {
        mReferenceSource = referenceSource;
    }

    public void addPointOfInterest(final String chromosome, final int position)
    {
        mPointsOfInterest.computeIfAbsent(chromosome, c -> new HashSet<>()).add(position);
    }

    /**
     * Will be rendered as a pinned sequence at the top of the table. The table will only show the bounds of header sequences.
     * Header sequences should all at least partially overlap or rendering issues may result.
     */
    public void addHeaderSequence(final AlignedSequence sequence)
    {
        mHeaderSequences.add(sequence);
    }

    public void addSequence(final AlignedSequence sequence)
    {
        addSequence(sequence, "default");
    }

    public void addSequence(final AlignedSequence sequence, final String cssClass)
    {
        mSequences.add(Pair.of(cssClass, sequence));
    }

    public void addSequence(final Sequence sequence, final int offset, final String cssClass)
    {
        mRawSequences.add(Triple.of(cssClass, offset, sequence));
        mMinOffset = Math.min(mMinOffset, offset);
    }

    private List<Alignment> addToLegend(final List<Alignment> legend, final AlignedSequence sequence)
    {
        if(legend.isEmpty())
            return sequence.getAlignmentBlocks();

        // TODO: This

        return legend;
    }

    private List<Alignment> computeLegend()
    {
        List<Alignment> legend = List.of();
        for(final AlignedSequence sequence : mHeaderSequences)
            legend = addToLegend(legend, sequence);

        for(final Pair<String, AlignedSequence> pair : mSequences)
            legend = addToLegend(legend, pair.getRight());

        return legend;
    }

    private Sequence createReference(final List<Alignment> legend)
    {
        final StringBuilder reference = new StringBuilder();

        for(final Alignment block : legend)
        {
            if(block.isUnmapped())
            {
                reference.append("B".repeat(block.Length));
            }
            else
            {
                final int endPositionInclusive =
                        Math.min(mReferenceSource.getChromosomeLength(block.Chromosome), block.ReferenceStartPosition + block.Length - 1);
                String sequence = mReferenceSource.getBaseString(block.Chromosome, block.ReferenceStartPosition, endPositionInclusive);
                if(block.Inverted)
                    sequence = SequenceUtil.reverseComplement(sequence);

                if(sequence.length() != block.Length)
                    sequence += "X".repeat(block.Length - sequence.length());

                reference.append(sequence);
            }
        }

        final String referenceSequence = reference.toString();
        final byte[] qualities = new byte[referenceSequence.length()];
        Arrays.fill(qualities, (byte) 37);
        return new SimpleSequence("Reference", referenceSequence.getBytes(), qualities);
    }

    @VisibleForTesting
    AlignedSequence realignSequenceToLegend(final AlignedSequence sequence, final List<Alignment> legend)
    {
        if(sequence.getAlignmentBlocks().isEmpty() || !sequence.getAlignmentBlocks().get(0).isUnmapped())
            return sequence;

        // remap the unmapped portion against the legend
        final Optional<Alignment> firstMappedOpt = sequence.getAlignmentBlocks().stream()
                .filter(block -> !block.isUnmapped())
                .findFirst();
        if(firstMappedOpt.isEmpty()) // TODO: Allow a sequence to have a suggested alignment vs a header sequence
            return sequence; // No mapping at all

        final Alignment firstMapped = firstMappedOpt.get();
        final int firstMappedIndex = sequence.getAlignmentBlocks().indexOf(firstMapped);

        final int clippedBases = firstMapped.SequenceStartPosition - 1;
        if(clippedBases <= 0)
            return sequence;

        final List<Alignment> newSequenceAlignment = new ArrayList<>();

        for(int i = 0; i < legend.size(); i++)
        {
            final Alignment legendElement = legend.get(i);
            if(legendElement.isUnmapped() || !legendElement.includes(firstMapped.Chromosome, firstMapped.ReferenceStartPosition))
                continue;

            final int offsetVsLegendElement = firstMapped.ReferenceStartPosition - legendElement.ReferenceStartPosition;
            if(offsetVsLegendElement == 0)
            {
                if(i == 0)
                    return sequence; // Can't unclip

                final Alignment previousElement = legend.get(i - 1);
                if(previousElement.isUnmapped())
                    return sequence;
                final int basesToUnclip = Math.min(clippedBases, previousElement.Length);
                final int leftoverClipBases = clippedBases - basesToUnclip;
                if(leftoverClipBases > 0)
                    newSequenceAlignment.add(Alignment.unmapped(leftoverClipBases));

                final int referenceStartPosition = previousElement.ReferenceStartPosition - basesToUnclip;
                final int readStartPosition = firstMapped.SequenceStartPosition - basesToUnclip;

                newSequenceAlignment.add(new Alignment(previousElement.Chromosome,
                        referenceStartPosition, readStartPosition, basesToUnclip, previousElement.Inverted, 60));

                for(int j = firstMappedIndex; j < sequence.getAlignmentBlocks().size(); j++)
                    newSequenceAlignment.add(sequence.getAlignmentBlocks().get(j));

                final AlignedSequence newAlignment =
                        new BasicAlignedSequence(sequence.getName(), sequence.getBases(), sequence.getBaseQuality(), newSequenceAlignment);

                return leftoverClipBases > 0 ? realignSequenceToLegend(newAlignment, legend) : newAlignment;
            }
            else
            {
                final int basesToUnclip = Math.min(clippedBases, offsetVsLegendElement);
                final int leftoverClipBases = clippedBases - basesToUnclip;
                if(leftoverClipBases > 0)
                    newSequenceAlignment.add(Alignment.unmapped(leftoverClipBases));

                newSequenceAlignment.add(new Alignment(legendElement.Chromosome,
                        firstMapped.ReferenceStartPosition - (firstMapped.Inverted ? -basesToUnclip : basesToUnclip),
                        firstMapped.SequenceStartPosition - basesToUnclip,
                        firstMapped.Length + basesToUnclip, firstMapped.Inverted,
                        60));

                for(int j = firstMappedIndex + 1; j < sequence.getAlignmentBlocks().size(); j++)
                    newSequenceAlignment.add(sequence.getAlignmentBlocks().get(j));

                final AlignedSequence newAlignment =
                        new BasicAlignedSequence(sequence.getName(), sequence.getBases(), sequence.getBaseQuality(), newSequenceAlignment);

                return leftoverClipBases > 0 ? realignSequenceToLegend(newAlignment, legend) : newAlignment;
            }
        }

        // Didn't overlap legend.
        return sequence;
    }

    private int offsetInLegend(final AlignedSequence sequence, final List<Alignment> legend)
    {
        int unmappedLeft = 0;
        for(final Alignment alignment : sequence.getAlignmentBlocks())
        {
            if(alignment.isUnmapped())
            {
                unmappedLeft += alignment.Length;
                continue;
            }

            int legendLengthSoFar = 0;
            for(final Alignment legendAlignment : legend)
            {
                if(legendAlignment.includes(alignment.Chromosome, alignment.ReferenceStartPosition, alignment.Length))
                {
                    if(legendAlignment.Inverted)
                    {
                        final int legendEndPosition = legendAlignment.ReferenceStartPosition + legendAlignment.Length - 1;
                        final int sequenceEndPosition = alignment.ReferenceStartPosition + alignment.Length - 1;
                        return legendLengthSoFar + (legendEndPosition - sequenceEndPosition) - unmappedLeft;
                    }
                    else
                        return legendLengthSoFar + (alignment.ReferenceStartPosition - legendAlignment.ReferenceStartPosition) - (
                                alignment.SequenceStartPosition - 1);
                }

                legendLengthSoFar += legendAlignment.Length;
            }
        }
        return 0;
    }

    public void generate(final HTMLBuilder builder)
    {
        final List<Alignment> legend = computeLegend();
        final Sequence reference = createReference(legend);

        builder.appendStartTag("<table>");

        builder.appendStartTag("<thead>");
        appendLegend(builder, legend, null, null);
        appendSequence(builder, reference, 0, true, "default");

        for(final AlignedSequence sequence : mHeaderSequences)
            appendSequence(builder, legend, realignSequenceToLegend(sequence, legend), true, "default");
        builder.appendEndTag("</thead>");

        builder.appendStartTag("<tbody>");
        for(final Pair<String, AlignedSequence> pair : mSequences)
            appendSequence(builder, legend, realignSequenceToLegend(pair.getRight(), legend), false, pair.getLeft());
        for(final Triple<String, Integer, Sequence> pair : mRawSequences)
            appendSequence(builder, pair.getRight(), pair.getMiddle(), false, pair.getLeft());
        builder.appendEndTag("</tbody>");

        builder.appendEndTag("</table>");
    }

    private String wrap(final String contents, final boolean wrapInDiv)
    {
        return wrapInDiv ? "<div>" + contents + "</div>" : contents;
    }

    private void appendSequence(final HTMLBuilder builder, final List<Alignment> legend, final AlignedSequence sequence,
            final boolean wrapCellsInDiv, final String rowCSSClass)
    {
        builder.append("<tr class=\"sequence %s\" id=\"%s\">", rowCSSClass, sequence.getName());

        builder.append("<th class=\"label pinned\"><div>").append(sequence.getName()).append("</div></th>");
        int startOffset = offsetInLegend(sequence, legend);
        if(mMinOffset < 0)
            startOffset -= mMinOffset;
        if(startOffset > 0)
        {
            while(startOffset > 1000)
            {
                builder.append("<td class=\"B\" colspan=\"1000\"></td>");
                startOffset -= 1000;
            }
            builder.append("<td class=\"B\" colspan=\"").append(startOffset).append("\"></td>");
        }

        for(final Alignment alignment : sequence.getAlignmentBlocks())
        {
            if(-startOffset >= alignment.Length)
            {
                startOffset += alignment.Length;
                continue;
            }

            if(alignment.Chromosome.equals("-"))
            {
                builder.append("<td class=\"B\" colspan=\"%s\">", alignment.Length);
                continue;
            }

            if(alignment.SequenceStartPosition <= 0)
            {
                int length = alignment.Length;
                if(startOffset < 0)
                {
                    final int consume = Math.min(-startOffset, length);
                    startOffset += consume;
                    length -= consume;
                }
                if(length > 0)
                    builder.append(("<td class=\"D\">" + wrap("*", wrapCellsInDiv) + "</td>").repeat(length));
                continue;
            }

            for(int offset = 0; offset < alignment.Length; offset++)
            {
                if(startOffset < 0)
                {
                    startOffset++;
                    continue;
                }

                final int readIndex = alignment.SequenceStartPosition + offset - 1;
                if(readIndex >= sequence.getBases().length)
                    break;

                final char base = (char) sequence.getBases()[readIndex];
                final byte quality = sequence.getBaseQuality()[readIndex];

                final char cssClass = base == ' ' ? 'B' : base;
                final String qualityClass = quality < 11 ? " lowq" : quality < 30 ? " mediumq" : "";
                builder.append("<td class=\"").append(cssClass).append(qualityClass).append("\">")
                        .append(wrap(String.valueOf(base), wrapCellsInDiv))
                        .append("</td>");
            }
        }

        builder.append("</tr>\n");
    }

    private void appendSequence(final HTMLBuilder builder, final Sequence sequence, int startOffset,
            final boolean wrapCellsInDiv, final String rowCSSClass)
    {
        if(mMinOffset < 0)
            startOffset -= mMinOffset;

        builder.append("<tr class=\"sequence %s\" id=\"%s\">", rowCSSClass, sequence.getName());

        builder.append("<th class=\"label pinned\"><div>").append(sequence.getName()).append("</div></th>");
        if(startOffset > 0)
        {
            while(startOffset > 1000)
            {
                builder.append("<td class=\"B\" colspan=\"1000\">").append(wrap("B", wrapCellsInDiv)).append("</td>");
                startOffset -= 1000;
            }
            builder.append("<td class=\"B\" colspan=\"").append(startOffset).append("\">")
                    .append(wrap("B", wrapCellsInDiv)).append("</td>");
        }

        int blankLength = 0;
        for(int i = 0; i < sequence.getLength(); i++)
        {
            final char base = (char) sequence.getBases()[i];
            if(base == 'B' || base == '-')
            {
                blankLength++;
                continue;
            }

            if(blankLength > 0)
            {
                builder.append("<td class=\"B\" colspan=\"%s\">", blankLength);
                blankLength = 0;
            }
            if(base == '*')
            {
                builder.append(("<td class=\"D\">" + wrap("*", wrapCellsInDiv) + "</td>"));
                continue;
            }
            final byte quality = sequence.getBaseQuality()[i];
            final char cssClass = base == ' ' ? 'B' : base;
            final String qualityClass = quality < 11 ? " lowq" : quality < 30 ? " mediumq" : "";
            builder.append("<td class=\"").append(cssClass).append(qualityClass).append("\">")
                    .append(wrap(String.valueOf(base), wrapCellsInDiv))
                    .append("</td>");
        }

        builder.append("</tr>\n");
    }

    private void appendChromosomeLabelCell(final StringBuilder sb, final String chromosomeName, int span)
    {
        final String chromosomeLabel = span > 11 ? "Chromosome " + chromosomeName : chromosomeName;
        if(span > 100)
        {
            appendChromosomeLabelCellWithClass(sb, "norightborder", chromosomeLabel, 80);
            span -= 80;
            while(span > 100)
            {
                appendChromosomeLabelCellWithClass(sb, "noleftborder norightborder", chromosomeLabel, 80);
                span -= 80;
            }
            appendChromosomeLabelCellWithClass(sb, "noleftborder", chromosomeLabel, span);
        }
        else
            sb.append("<th colspan=").append(span).append("><div>").append(chromosomeLabel).append("</div></th>");
    }

    private void appendChromosomeLabelCellWithClass(final StringBuilder sb, final String className, final String label, final int span)
    {
        sb.append("<th class=\"")
                .append(className)
                .append("\" colspan=\"")
                .append(span)
                .append("\"><div>")
                .append(label)
                .append("</div></th>");
    }

    private void appendMajorNumberLabelCell(final StringBuilder sb, final int position, final int span)
    {
        String label = span <= 2 ? "*" : String.valueOf(position);
        if(label.length() > span)
            label = "*" + label.substring(label.length() - span + 1);

        sb.append("<th colspan=").append(span).append("><div>").append(label).append("</div></th>");
    }

    /**
     * The legend contains space for 2 actions, optionally included as raw-html.
     */
    private void appendLegend(final HTMLBuilder builder, final List<Alignment> legend, @Nullable final String action1,
            @Nullable final String action2)
    {
        final StringBuilder chromosomeRow = new StringBuilder();
        chromosomeRow.append("<tr class=\"legend chromosome\"><th class=\"pinned")
                .append(action1 == null ? " blank" : "").append("\"><div>");
        chromosomeRow.append(Objects.requireNonNullElse(action1, "X"));
        chromosomeRow.append("</div></th>");

        final StringBuilder majorNumberRow = new StringBuilder();
        majorNumberRow.append("<tr class=\"legend major\"><th class=\"pinned")
                .append(action2 == null ? " blank" : "").append("\"><div>");
        majorNumberRow.append(Objects.requireNonNullElse(action2, "X"));
        majorNumberRow.append("</div></th>");

        final StringBuilder minorNumberRow = new StringBuilder();
        minorNumberRow.append("<tr class=\"legend minor\"><th class=\"pinned blank\"><div>X</div></th>");

        if(mMinOffset < 0)
        {
            final String padding = String.format("<td colspan=\"%s\"></td>", -mMinOffset);
            chromosomeRow.append(padding);
            majorNumberRow.append(padding);
            minorNumberRow.append(padding);
        }

        @Nullable
        String currentChromosome = null;
        int chromosomeSpan = 0;
        @Nullable
        Integer currentMajor = null;
        int majorSpan = 0;
        for(final Alignment alignment : legend)
        {
            final int referenceLeftPosition = alignment.Inverted
                    ? alignment.ReferenceStartPosition + alignment.Length - 1
                    : alignment.ReferenceStartPosition;

            final boolean isNoReference = alignment.Chromosome.equals("*") || alignment.Chromosome.equals("?");
            if(!alignment.Chromosome.equals(currentChromosome) && (!isNoReference || currentChromosome == null))
            {
                if(currentChromosome != null && chromosomeSpan > 0)
                {
                    appendChromosomeLabelCell(chromosomeRow, currentChromosome, chromosomeSpan);
                    if(majorSpan > 0)
                        appendMajorNumberLabelCell(majorNumberRow, Objects.requireNonNullElse(currentMajor, 0), majorSpan);
                }

                currentChromosome = alignment.Chromosome;
                currentMajor = isNoReference ? null : referenceLeftPosition / 10 * 10;
                chromosomeSpan = 0;
                majorSpan = 0;
            }

            for(int i = 0; i < alignment.Length; i++)
            {
                if(isNoReference)
                {
                    minorNumberRow.append("<td><div>*</div></td>");
                    majorSpan++;
                    continue;
                }

                final int thisPosition = alignment.Inverted ? referenceLeftPosition - i : referenceLeftPosition + i;
                final int thisMajor = thisPosition / 10 * 10;
                if(thisMajor != currentMajor)
                {
                    appendMajorNumberLabelCell(majorNumberRow, currentMajor, majorSpan);
                    currentMajor = thisMajor;
                    majorSpan = 0;
                }
                majorSpan++;

                final boolean isPOI = mPointsOfInterest.getOrDefault(currentChromosome, Set.of()).contains(thisPosition);
                final String tdClass = isPOI ? " class=\"junction\"" : "";
                minorNumberRow.append("<td").append(tdClass).append("><div>").append(thisPosition % 10).append("</div></td>");
            }

            chromosomeSpan += alignment.Length;
        }

        if(chromosomeSpan > 0)
        {
            appendChromosomeLabelCell(chromosomeRow, currentChromosome, chromosomeSpan);
            assert majorSpan > 0;
            appendMajorNumberLabelCell(majorNumberRow, currentMajor == null ? 0 : currentMajor, majorSpan);
        }

        builder.append(chromosomeRow).append("</tr>\n");
        builder.append(majorNumberRow).append("</tr>\n");
        builder.append(minorNumberRow).append("</tr>\n");
    }
}
