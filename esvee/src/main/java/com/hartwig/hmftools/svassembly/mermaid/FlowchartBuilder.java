package com.hartwig.hmftools.svassembly.mermaid;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.IdentityHashMap;
import java.util.Map;
import java.util.stream.Collectors;

import org.jetbrains.annotations.Nullable;

public class FlowchartBuilder
{
    final Map<Object, Node> nodes = new IdentityHashMap<>();
    final Map<String, String> registeredStyles = new HashMap<>();
    final Map<Link<Object>, Link<Object>> links = new HashMap<>();
    private int nextId = 0;

    public void addStyle(final String name, final String styleText)
    {
        registeredStyles.put(name, styleText);
    }

    /**
     * Updates the label associated with the given object identity
     *
     * @param identity user data associated with this node, identifies it to be used with links later.
     */
    public boolean addNode(final Object identity, final String label)
    {
        return addNode(identity, label, null);
    }

    public boolean addNode(final Object identity, final String label, @Nullable final String style)
    {
        final boolean[] existing = new boolean[1];
        nodes.compute(identity, (id, old) ->
        {
            if (old == null)
                return new Node(nextId++, label, style);
            else
            {
                existing[0] = true;
                old.Label = label;
                old.Style = style;
                return old;
            }
        });
        return !existing[0];
    }

    /** Uses the object identities registered using addNode */
    public void addLink(final Object node1, final Object node2)
    {
        addLink(node1, node2, null);
    }

    public void addLink(final Object node1, final Object node2, @Nullable final String label)
    {
        final var link = new Link<>(node1, node2, label);
        links.put(link, link);
    }

    public Flowchart build()
    {
        return new Flowchart(registeredStyles, new ArrayList<>(nodes.values()),
                links.values().stream()
                        .map(link -> new Link<>(nodes.get(link.Left), nodes.get(link.Right), link.Text))
                        .collect(Collectors.toList()));
    }
}
