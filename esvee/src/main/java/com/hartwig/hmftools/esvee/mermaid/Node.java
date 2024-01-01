package com.hartwig.hmftools.esvee.mermaid;

import org.jetbrains.annotations.Nullable;

class Node
{
    public final int Id;
    public String Label;
    @Nullable
    public String Style;

    public Node(final int id, final String label, @Nullable final String style)
    {
        Id = id;
        Label = label;
        Style = style;
    }
}
