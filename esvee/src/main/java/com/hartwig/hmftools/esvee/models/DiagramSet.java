package com.hartwig.hmftools.esvee.models;

import java.util.ArrayList;
import java.util.List;

import com.hartwig.hmftools.esvee.mermaid.Flowchart;

import org.apache.commons.lang3.tuple.Pair;

public class DiagramSet
{
    public final String Label;
    public final List<Pair<String, Flowchart>> Diagrams;

    public DiagramSet(final String label)
    {
        Label = label;
        Diagrams = new ArrayList<>();
    }

    public void add(final String label, final Flowchart diagram)
    {
        Diagrams.add(Pair.of(label, diagram));
    }
}
