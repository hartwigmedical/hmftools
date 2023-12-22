package com.hartwig.hmftools.svassembly.mermaid;

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Deque;
import java.util.HashMap;
import java.util.IdentityHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.hartwig.hmftools.svassembly.output.html.HTMLBuilder;

import org.apache.commons.lang3.tuple.Pair;

public enum FlowchartSVGGenerator
{
    ;

    private static final int NODE_RADIUS = 33;
    private static final int SCALE = 100;
    private static final int FONT_SIZE = 26;
    private static final int EDGE_WIDTH = 3;

    // Calculate the max depth for every node -- this is how far "right" this node is
    // Traverse the graph, track the "height" we're at.

    private static Map<Node, Integer> computeMaxDepth(final List<Node> startNodes,
            final Map<Node, List<Node>> linksBySource)
    {
        final Map<Node, Integer> maxDepth = new IdentityHashMap<>();

        int depth = 0;
        Collection<Node> toProcess = startNodes;
        while (!toProcess.isEmpty())
        {
            final Collection<Node> next = Collections.newSetFromMap(new IdentityHashMap<>());
            for(final Node node : toProcess)
            {
                maxDepth.put(node, depth);
                next.addAll(linksBySource.getOrDefault(node, List.of()));
            }
            depth++;
            toProcess = next;

            if (depth > 10000)
                throw new IllegalStateException("Graph contains a cycle or is too large");
        }

        return maxDepth;
    }

    private static List<Node> computeStartNodes(final List<Node> nodes, final Map<Node, List<Node>> linksBySource)
    {
        final Set<Node> destinationNodes = Collections.newSetFromMap(new IdentityHashMap<>());
        linksBySource.values().forEach(destinationNodes::addAll);

        final List<Node> startNodes = new ArrayList<>();
        for(final Node node : nodes)
            if (!destinationNodes.contains(node))
                startNodes.add(node);
        return startNodes;
    }

    private static Map<Node, List<Node>> createLinksBySource(final List<Link<Node>> links)
    {
        final Map<Node, List<Node>> linksBySource = new IdentityHashMap<>();
        for(final Link<Node> link : links)
        {
            final Node left = link.Left;
            final Node right = link.Right;
            linksBySource.computeIfAbsent(left, n -> new ArrayList<>()).add(right);
        }
        return linksBySource;
    }

    private static Map<Node, Pair<Integer, Integer>> computeNodePositions(
            final Map<Node, List<Node>> linksBySource,
            final List<Node> startNodes)
    {
        final Map<Node, Integer> maxDepth = computeMaxDepth(startNodes, linksBySource);

        final Map<Integer, Integer> heightAtDepth = new HashMap<>();
        final Set<Node> rendered = Collections.newSetFromMap(new IdentityHashMap<>());
        final Deque<Node> toProcess = new ArrayDeque<>();
        toProcess.addAll(startNodes);

        final Map<Node, Pair<Integer, Integer>> nodePositions = new IdentityHashMap<>();
        int minHeight = 0;
        while (!toProcess.isEmpty())
        {
            final Node node = toProcess.pollLast();
            if (!rendered.add(node))
                continue;

            final int currentMinHeight = minHeight;

            final int positionX = maxDepth.getOrDefault(node, 0);
            final int positionY = heightAtDepth.compute(positionX, (k, oldValue) -> oldValue == null
                    ? currentMinHeight
                    : Math.max(oldValue + 1, currentMinHeight));
            minHeight = Math.max(positionY, currentMinHeight);
            nodePositions.put(node, Pair.of(positionX * SCALE, positionY * SCALE));

            linksBySource.getOrDefault(node, List.of()).forEach(toProcess::addLast); // Depth-first traversal
        }
        return nodePositions;
    }

    private static Pair<Integer, Integer> generate(final HTMLBuilder nodesBuilder, final HTMLBuilder edgesBuilder, final Flowchart flowchart)
    {
        final List<Node> nodes = flowchart.createNodes();
        final List<Link<Node>> links = flowchart.createLinks(nodes);

        final Map<Node, List<Node>> linksBySource = createLinksBySource(links);
        final List<Node> startNodes = computeStartNodes(nodes, linksBySource);

        final Map<Node, Pair<Integer, Integer>> nodePositions = computeNodePositions(linksBySource, startNodes);

        renderEdges(edgesBuilder, linksBySource, nodePositions);
        renderNodes(nodesBuilder, nodePositions);

        final int width = nodePositions.values().stream().mapToInt(Pair::getLeft).max().orElseThrow() + SCALE;
        final int height = nodePositions.values().stream().mapToInt(Pair::getRight).max().orElseThrow() + SCALE;

        //noinspection SuspiciousNameCombination
        return Pair.of(width, height);
    }

    private static void renderEdges(final HTMLBuilder edgesBuilder,
            final Map<Node, List<Node>> linksBySource,
            final Map<Node, Pair<Integer, Integer>> nodePositions)
    {
        for(final Map.Entry<Node, List<Node>> entry : linksBySource.entrySet())
        {
            final Node source = entry.getKey();
            for (final Node destination : entry.getValue())
            {
                final Pair<Integer, Integer> startPosition = nodePositions.get(source);
                final Pair<Integer, Integer> endPosition = nodePositions.get(destination);

                final int startX = startPosition.getLeft();
                final int startY = startPosition.getRight();

                final int endX = endPosition.getLeft();
                final int endY = endPosition.getRight();

                edgesBuilder.appendStartTag("<g class=\"edge\">");
                if (startY == endY)
                    edgesBuilder.append("<line x1=\"%s\" y1=\"%s\" x2=\"%s\" y2=\"%s\"/>", startX, startY, endX, endY);
                else
                {
                    final double partialX;
                    final boolean isGoingDown = endY > startY;
                    if (isGoingDown)
                        partialX = ((double) endX - startX) / 4.0 + startX;
                    else
                        partialX = ((double) endX - startX) * 3.0 / 4.0 + startX;

                    edgesBuilder.append("<path d=\"M %s %s C %s %s, %s %s, %s %s\"/>",
                            startX, startY,
                            partialX, isGoingDown ? endY : startY,
                            partialX, isGoingDown ? endY : startY,
                            endX, endY);
                }

                edgesBuilder.appendEndTag("</g>");
            }
        }
    }

    private static void renderNodes(final HTMLBuilder nodesBuilder,
            final Map<Node, Pair<Integer, Integer>> nodePositions)
    {
        for(final Node node : nodePositions.keySet())
        {
            final var position = nodePositions.get(node);
            final int x = position.getLeft();
            final int y = position.getRight();

            nodesBuilder.appendStartTag("<g class=\"node" + (node.Style == null ? "" : " " + node.Style) + "\">");
            nodesBuilder.append("<circle cx=\"%s\" cy=\"%s\" r=\"%s\"></circle>\n", x, y, NODE_RADIUS);
            nodesBuilder.append("<text x=\"%s\" y=\"%s\">%s</text>\n", x, y, node.Label);
            nodesBuilder.appendEndTag("</g>");
        }
    }

    private static CharSequence createStyle(final Flowchart flowchart)
    {
        final StringBuilder styleBuilder = new StringBuilder();
        styleBuilder.append(".flowchart .edge {\n");
        styleBuilder.append("\tstroke: black;\n");
        styleBuilder.append("\tfill: none;\n");
        styleBuilder.append(String.format("\tstroke-width: %spx;\n", EDGE_WIDTH));
        styleBuilder.append("}\n");
        styleBuilder.append(".flowchart .node {\n");
        styleBuilder.append("}\n");
        styleBuilder.append(".flowchart text {\n");
        styleBuilder.append("\ttext-anchor: middle;\n");
        styleBuilder.append("\ttext-align: center;\n");
        styleBuilder.append("\tfill: black;\n");
        styleBuilder.append("\tdominant-baseline: central;\n");
        styleBuilder.append("\tfont-family: monospace;\n");
        styleBuilder.append(String.format("\tfont-size: %spx;\n", FONT_SIZE));
        styleBuilder.append("}\n");
        for(final Map.Entry<String, String> entry : flowchart.Styles.entrySet())
        {
            final String name = entry.getKey();
            final String text = entry.getValue();
            styleBuilder.append(".").append(name).append(" {\n");
            styleBuilder.append("\t").append(text.replace(",", ";")).append(";\n");
            styleBuilder.append("}\n");
        }

        return styleBuilder;
    }

    public static String generate(final Flowchart flowchart, final String customStyle, final int lineHeight)
    {
        final HTMLBuilder builder = new HTMLBuilder();

        final HTMLBuilder nodesBuilder = new HTMLBuilder();
        final HTMLBuilder edgesBuilder = new HTMLBuilder();

        final Pair<Integer, Integer> size = generate(nodesBuilder, edgesBuilder, flowchart);
        final int width = size.getLeft();
        final int height = size.getRight();

        builder.appendStartTag("<svg xmlns=\"http://www.w3.org/2000/svg\""
                + " role=\"graphics-document document\" class=\"flowchart\" style=\"height: %spx; %s\" viewBox=\"-%s -%s %s %s\">",
                ((double) height / SCALE) * lineHeight,
                customStyle,
                SCALE / 2, SCALE / 2, width, height);

        builder.appendStartTag("<style>");
        builder.append(createStyle(flowchart));
        builder.appendEndTag("</style>");

        builder.appendStartTag("<g class=\"root\">");

        builder.appendStartTag("<g class=\"edges\">");
        builder.append(edgesBuilder.getStringBuilder());
        builder.appendEndTag("</g>");

        builder.appendStartTag("<g class=\"nodes\">");
        builder.append(nodesBuilder.getStringBuilder());
        builder.appendEndTag("</g>");

        builder.appendEndTag("</g>");
        builder.appendEndTag("</svg>");

        return builder.getStringBuilder().toString();
    }
}
