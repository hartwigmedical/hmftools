package com.hartwig.hmftools.esvee.assembly;

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.List;
import java.util.Queue;

import com.hartwig.hmftools.esvee.mermaid.Flowchart;
import com.hartwig.hmftools.esvee.mermaid.FlowchartBuilder;
import com.hartwig.hmftools.esvee.models.Record;
import com.hartwig.hmftools.esvee.util.StringCache;

import org.jetbrains.annotations.Nullable;

public class Node
{
    public static class Support implements Comparable<Support>
    {
        public final com.hartwig.hmftools.esvee.models.Record Record;
        public final int ReadIndex;

        public Support(final Record record, final int readIndex)
        {
            Record = record;
            ReadIndex = readIndex;
        }

        @Override
        public boolean equals(final Object obj)
        {
            if(!(obj instanceof Node.Support))
                return false;

            final Node.Support other = (Node.Support) obj;
            return other.Record.equals(Record) && other.ReadIndex == ReadIndex;
        }

        @Override
        public int hashCode()
        {
            return Record.hashCode();
        }

        @Override
        public int compareTo(final Node.Support other)
        {
            if(Record == other.Record)
                return 0;

            int compare = Integer.compare(Record.getAlignmentStart(), other.Record.getAlignmentStart());
            if(compare != 0)
                return compare;

            compare = Boolean.compare(Record.isFirstOfPair(), other.Record.isFirstOfPair());
            if(compare != 0)
                return compare;

            return Record.getName().compareTo(other.Record.getName());
        }
    }

    public final char Base;
    public int Quality;
    public int MaxQuality;
    public List<Support> Support;
    @Nullable
    Node nextA, nextT, nextC, nextG, nextX;
    @Nullable
    private List<Node> cachedSuccessors;
    private int cachedSupportDepth = -1;

    public Node(final char base)
    {
        if(base == 'A' || base == 'T' || base == 'C' || base == 'G' || base == 'S' || base == 'X')
            Base = base;
        else
            Base = 'X';
    }

    public int successorCount()
    {
        return successors().size();
    }

    public List<Node> successors()
    {
        if(cachedSuccessors != null) return cachedSuccessors;

        final List<Node> nodes = new ArrayList<>(4);
        if(nextA != null)
            nodes.add(nextA);
        if(nextT != null)
            nodes.add(nextT);
        if(nextC != null)
            nodes.add(nextC);
        if(nextG != null)
            nodes.add(nextG);
        if(nextX != null)
            nodes.add(nextX);
        return cachedSuccessors = nodes;
    }

    public void removeSuccessor(final char base)
    {
        switch(base)
        {
            case 'A':
                nextA = null;
                break;
            case 'T':
                nextT = null;
                break;
            case 'C':
                nextC = null;
                break;
            case 'G':
                nextG = null;
                break;
            case 'X':
                nextX = null;
                break;
        }
        cachedSuccessors = null;
    }

    @Nullable
    public Node getNext(final char base)
    {
        switch(base)
        {
            case 'A':
                return nextA;
            case 'T':
                return nextT;
            case 'C':
                return nextC;
            case 'G':
                return nextG;
            case 'X':
                return nextX;
        }
        return null;
    }

    public void setNext(@Nullable final Node next)
    {
        if(next == null)
            return;

        cachedSuccessors = null;
        if(next.Base == 'A')
            nextA = next;
        else if(next.Base == 'T')
            nextT = next;
        else if(next.Base == 'C')
            nextC = next;
        else if(next.Base == 'G')
            nextG = next;
        else if(next.Base == 'X')
            nextX = next;
    }

    public void recomputeQuality()
    {
        int sum = 0;
        int max = 0;
        for(final Node.Support support : Support)
        {
            final int baseQuality = support.Record.getBaseQuality()[support.ReadIndex];
            sum += baseQuality;
            max = Math.max(max, baseQuality);
        }
        Quality = sum;
        MaxQuality = max;
    }

    public int supportDepth()
    {
        if(cachedSupportDepth != -1)
            return cachedSupportDepth;

        return cachedSupportDepth = Support.size();
    }

    public Node deepCopy()
    {
        final HeadNode headNode = new HeadNode();
        headNode.setNext(this);

        final HeadNode copy = headNode.deepCopy();
        return copy.successors().get(0);
    }

    public Flowchart toDiagram()
    {
        final FlowchartBuilder diagram = new FlowchartBuilder();
        diagram.addStyle("A", "fill:#00F000,color:#000000");
        diagram.addStyle("T", "fill:#F00000,color:#000000");
        diagram.addStyle("C", "fill:#003BFF,color:#000000");
        diagram.addStyle("G", "fill:#E0E000,color:#000000");
        diagram.addStyle("N", "fill:#FF00FF,color:#000000");
        diagram.addStyle("X", "fill:#FFFFFF,color:#000000");
        diagram.addStyle("S", "fill:#FFFFFF,color:#000000");

        final Queue<Node> nodes = new ArrayDeque<>();
        nodes.add(this);
        while(!nodes.isEmpty())
        {
            final Node node = nodes.poll();

            final String base = StringCache.of(node.Base);
            if(!diagram.addNode(node, StringCache.of(node.Quality), base))
                continue;

            for(final Node successor : node.successors())
            {
                nodes.add(successor);
                diagram.addLink(node, successor);
            }
        }

        return diagram.build();
    }
}
