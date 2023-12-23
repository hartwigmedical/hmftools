package com.hartwig.hmftools.esvee.mermaid;

import java.io.IOException;
import java.util.Arrays;
import java.util.Base64;
import java.util.List;
import java.util.Map;
import java.util.concurrent.CompletableFuture;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.zip.Deflater;

import com.hartwig.hmftools.esvee.util.StringCache;

import org.jetbrains.annotations.Nullable;

public class Flowchart
{
    public final Map<String, String> Styles;
    private final String[] mNodeLabels;
    private final String[] mNodeStyles;

    private final int[] mLinkLeft;
    private final int[] mLinkRight;
    private final String[] mLinkText;

    public Flowchart(final Map<String, String> styles, final List<Node> nodes, final List<Link<Node>> links)
    {
        Styles = styles;

        mNodeLabels = new String[nodes.size()];
        mNodeStyles = new String[nodes.size()];
        for(final Node node : nodes)
        {
            mNodeLabels[node.Id] = StringCache.tryDedupe(node.Label);
            mNodeStyles[node.Id] = StringCache.tryDedupe(node.Style);
        }

        mLinkLeft = new int[links.size()];
        mLinkRight = new int[links.size()];
        mLinkText = new String[links.size()];
        for(int i = 0; i < links.size(); i++)
        {
            final Link<Node> link = links.get(i);
            mLinkLeft[i] = link.Left.Id;
            mLinkRight[i] = link.Right.Id;
            mLinkText[i] = StringCache.tryDedupe(link.Text);
        }
    }

    public List<Node> createNodes()
    {
        return IntStream.range(0, mNodeLabels.length)
                .mapToObj(i -> new Node(i, mNodeLabels[i], mNodeStyles[i]))
                .collect(Collectors.toList());
    }

    public List<Link<Node>> createLinks(final List<Node> nodes)
    {
        return IntStream.range(0, mLinkLeft.length)
                .mapToObj(i -> new Link<>(nodes.get(mLinkLeft[i]), nodes.get(mLinkRight[i]), mLinkText[i]))
                .collect(Collectors.toList());
    }

    @Override
    public String toString()
    {
        return "Flowchart(" + mNodeLabels.length + " nodes, " + Styles.size() + " styles, " + mLinkLeft.length + " links)";
    }

    public StringBuilder build(final StringBuilder sb)
    {
        sb.append("flowchart LR\n");

        // Declare styles
        // https://mermaid.js.org/syntax/flowchart.html#styling-a-node
        for(final Map.Entry<String, String> styleEntry : Styles.entrySet())
            sb.append("    classDef ").append(styleEntry.getKey()).append(" ").append(styleEntry.getValue()).append("\n");

        final List<Node> nodes = createNodes();
        for(final Node node : nodes)
        {
            sb.append("    n").append(node.Id).append("((\"").append(node.Label).append("\"))");
            if(node.Style != null)
            {
                if(!Styles.containsKey(node.Style))
                    throw new IllegalStateException("Use of unregistered style " + node.Style);
                sb.append(":::").append(node.Style);
            }
            sb.append("\n");
        }

        for(final Link<Node> link : createLinks(nodes))
        {
            @Nullable
            final Node left = link.Left;
            @Nullable
            final Node right = link.Right;
            if(left == null || right == null)
                throw new IllegalStateException("Link uses unregistered nodes! From " + link.Left + " to " + link.Right);
            sb.append("    n").append(left.Id).append(" --- n").append(right.Id).append("\n");
        }

        return sb;
    }

    public String build()
    {
        return build(new StringBuilder()).toString();
    }

    public String toMarkdown()
    {
        final var sb = new StringBuilder();
        sb.append("```mermaid\n");
        build(sb);
        sb.append("\n```\n");
        return sb.toString();
    }

    /** Generates a compressed base64 representation of the diagram suitable for conversion via mermaid.live */
    public String toCompressedBase64()
    {
        final String forImage = "{\"code\":\""
                + build().replace("\"", "\\\"").replace("\n", "\\n")
                + "\",\"mermaid\":{\"theme\":\"default\"}}";

        return Base64.getUrlEncoder().encodeToString(compress(forImage.getBytes()));
    }

    /** Generates an SVG using the mmdc utility. */
    public CompletableFuture<Process> generateSVGUsingUtil(final String outputPath)
    {
        try
        {
            final Process process = new ProcessBuilder("/opt/homebrew/bin/mmdc", "-t", "dark", "-b", "transparent", "-o", outputPath)
                    .start();
            process.getOutputStream().write(build().getBytes());
            process.getOutputStream().flush();
            process.getOutputStream().close();
            return process.onExit();
        }
        catch(final IOException e)
        {
            throw new RuntimeException(e);
        }
    }

    public String toSVG()
    {
        return toSVG(null);
    }

    public String toSVG(@Nullable final Integer lineHeight)
    {
        return FlowchartSVGGenerator.generate(this, "", lineHeight == null ? 32 : lineHeight);
    }

    private static byte[] compress(final byte[] input)
    {
        final Deflater deflater = new Deflater();
        deflater.setInput(input);
        deflater.finish();

        final byte[] compressionBuffer = new byte[1024 * 1024];
        final int size = deflater.deflate(compressionBuffer);
        deflater.end();

        return Arrays.copyOfRange(compressionBuffer, 0, size);
    }
}
