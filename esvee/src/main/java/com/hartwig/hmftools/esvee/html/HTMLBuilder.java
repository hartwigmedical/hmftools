package com.hartwig.hmftools.esvee.html;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.hartwig.hmftools.esvee.util.Counter;

import org.apache.commons.lang3.tuple.Pair;

@SuppressWarnings({ "UnusedReturnValue", "unused" })
public class HTMLBuilder
{
    private final StringBuilder mSB;
    private String mIndent;
    private boolean mIsReadyForIndent;

    public HTMLBuilder()
    {
        this(new StringBuilder(), "");
    }

    public HTMLBuilder(final StringBuilder sb, final String indent)
    {
        mSB = sb;
        mIndent = indent;
        mIsReadyForIndent = true;
    }

    public StringBuilder getStringBuilder()
    {
        return mSB;
    }

    public HTMLBuilder append(final char value)
    {
        if(mIsReadyForIndent)
            mSB.append(mIndent);

        mIsReadyForIndent = false;
        mSB.append(value);
        return this;
    }

    public HTMLBuilder append(final int value)
    {
        if(mIsReadyForIndent)
            mSB.append(mIndent);

        mIsReadyForIndent = false;
        mSB.append(value);
        return this;
    }

    public HTMLBuilder append(final String html)
    {
        if(mIsReadyForIndent)
            mSB.append(mIndent);

        mSB.append(html);

        mIsReadyForIndent = html != null && html.endsWith("\n");
        return this;
    }

    public HTMLBuilder append(final Object object)
    {
        return append(String.valueOf(object));
    }

    public HTMLBuilder append(final String format, final Object... args)
    {
        return append(String.format(format, args));
    }

    public HTMLBuilder appendStartTag(final String format, final Object... args)
    {
        return appendStartTag(String.format(format, args));
    }

    public HTMLBuilder appendStartTag(final String tag)
    {
        append(tag).append("\n");
        increaseIndent();
        return this;
    }

    public HTMLBuilder appendEndTag(final String tag)
    {
        decreaseIndent();
        append(tag).append("\n");
        return this;
    }

    public HTMLBuilder increaseIndent()
    {
        mIndent += "\t";
        return this;
    }

    public HTMLBuilder decreaseIndent()
    {
        mIndent = mIndent.substring(1);
        return this;
    }

    public HTMLBuilder appendCountersTable(final List<Counter> counters, final int columns)
    {
        return appendSimpleAttributeTable(counters.stream()
                .map(counter -> Pair.of(counter.Name, counter.formatValue()))
                .collect(Collectors.toList()), columns);
    }

    public HTMLBuilder appendSimpleAttributeTable(final Map<String, ?> pairs)
    {
        return appendSimpleAttributeTable(pairs, 1);
    }

    public HTMLBuilder appendSimpleAttributeTable(final Map<String, ?> pairs, final int columns)
    {
        return appendSimpleAttributeTable(pairs.entrySet().stream()
                .map(entry -> Pair.of(entry.getKey(), entry.getValue()))
                .collect(Collectors.toList()), columns);
    }

    public HTMLBuilder appendSimpleAttributeTable(final List<Pair<String, ?>> pairs)
    {
        return appendSimpleAttributeTable(pairs, 1);
    }

    public HTMLBuilder appendSimpleAttributeTable(final List<Pair<String, ?>> pairs, final int columns)
    {
        appendStartTag("<table>");
        appendStartTag("<tbody>");

        final List<String> rowContents = new ArrayList<>();
        for(final Pair<String, ?> pair : pairs)
            rowContents.add(String.format("<th>%s</th><td>%s</td>", pair.getKey(), pair.getValue()));

        final int rowCount = (int) Math.ceil((double) rowContents.size() / columns);
        for(int i = 0; i < rowCount; i++)
        {
            final StringBuilder row = new StringBuilder();
            for(int j = 0; j < columns; j++)
            {
                final int index = i + (j * rowCount);
                if(index < rowContents.size())
                    row.append(rowContents.get(index));
            }

            append("<tr>").append(row.toString()).append("</tr>\n");
        }

        appendEndTag("</tbody>");
        appendEndTag("</table>");

        return this;
    }

    @Override
    public String toString()
    {
        return getStringBuilder().toString();
    }
}
