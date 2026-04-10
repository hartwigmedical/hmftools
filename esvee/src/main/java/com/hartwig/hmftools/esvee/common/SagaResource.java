package com.hartwig.hmftools.esvee.common;

import static java.util.function.UnaryOperator.identity;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import com.hartwig.hmftools.common.genome.region.Orientation;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;

public class SagaResource
{
    private final String mFastaPath;
    private final List<AssemblyMetadata> mAssemblies;

    public SagaResource(final String fastaPath)
    {
        mFastaPath = fastaPath;
        mAssemblies = loadAssemblies();
    }

    public String bwaIndexImagePath()
    {
        return mFastaPath + ".img";
    }

    public Map<String, AssemblyMetadata> assembliesByFastaLabel()
    {
        return mAssemblies.stream().collect(Collectors.toMap(AssemblyMetadata::fastaLabel, identity()));
    }

    public Map<String, Variant> variantsById()
    {
        return mAssemblies.stream()
                .map(AssemblyMetadata::variant)
                .collect(Collectors.toMap(Variant::id, identity()));
    }

    public Map<String, List<IndexedBreakend>> breakendsByChromosome()
    {
        return mAssemblies.stream()
                .flatMap(assembly ->
                        assembly.variant.breakends().map(breakend ->
                                new IndexedBreakend(breakend, assembly.variant.id())))
                .collect(Collectors.groupingBy(IndexedBreakend::chromosome));
    }

    public record Breakend(
            String chromosome,
            int position,
            Orientation orientation
    )
    {
        public static Breakend fromString(final String string)
        {
            String[] parts = string.split(":");
            if(parts.length != 3)
            {
                throw new IllegalArgumentException("Invalid breakend string: " + string);
            }
            String chromosome = parts[0];
            int position = Integer.parseInt(parts[1]);
            Orientation orientation = Orientation.fromByteStr(parts[2]);
            return new Breakend(chromosome, position, orientation);
        }
    }

    public record Variant(
            String id,
            Breakend breakend1,
            Breakend breakend2
    )
    {
        public Stream<Breakend> breakends()
        {
            return Stream.of(breakend1, breakend2);
        }
    }

    public record AssemblyMetadata(
            String fastaLabel,
            Variant variant,
            int junctionOffset1,
            int junctionOffset2
    )
    {
        public static AssemblyMetadata fromFastaLabel(final String fastaLabel)
        {
            String[] parts = fastaLabel.split("\\|");
            if(parts.length != 5)
            {
                throw new IllegalArgumentException("Invalid fasta label: " + fastaLabel);
            }
            String id = parts[0];
            Breakend breakend1 = Breakend.fromString(parts[1]);
            Breakend breakend2 = Breakend.fromString(parts[2]);
            int junctionOffset1 = Integer.parseInt(parts[3]);
            int junctionOffset2 = Integer.parseInt(parts[4]);
            Variant variant = new Variant(id, breakend1, breakend2);
            return new AssemblyMetadata(fastaLabel, variant, junctionOffset1, junctionOffset2);
        }
    }

    public record IndexedBreakend(
            Breakend breakend,
            String variantId
    )
    {
        public String chromosome()
        {
            return breakend.chromosome();
        }

        public int position()
        {
            return breakend.position();
        }
    }

    private List<AssemblyMetadata> loadAssemblies()
    {
        List<AssemblyMetadata> assemblies;
        try(IndexedFastaSequenceFile fasta = loadFasta())
        {
            assemblies = fasta.getSequenceDictionary().getSequences().stream()
                    .map(seq -> AssemblyMetadata.fromFastaLabel(seq.getSequenceName()))
                    .toList();
        }
        catch(IOException e)
        {
            // Failed to load/read/close fasta file.
            throw new RuntimeException("Failed to load SAGA resource", e);
        }

        // Variants are keyed by ID later, so ensure the IDs are unique.
        long uniqueVariantIdCount = assemblies.stream().map(assembly -> assembly.variant().id()).distinct().count();
        if(uniqueVariantIdCount != assemblies.size())
        {
            throw new RuntimeException("Duplicate variant IDs in SAGA resource");
        }

        return assemblies;
    }

    private IndexedFastaSequenceFile loadFasta() throws IOException
    {
        return new IndexedFastaSequenceFile(new File(mFastaPath));
    }
}
