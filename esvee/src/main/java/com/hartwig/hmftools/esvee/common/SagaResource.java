package com.hartwig.hmftools.esvee.common;

import static java.util.function.UnaryOperator.identity;

import static com.hartwig.hmftools.esvee.assembly.AssemblyConfig.SV_LOGGER;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import com.hartwig.hmftools.common.genome.region.Orientation;
import com.hartwig.hmftools.common.region.BasePosition;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;

public class SagaResource
{
    private final String mFastaPath;
    private final List<AssemblyMetadata> mAssemblies;
    private final Map<String, Variant> mVariantsById;
    private final Map<String, List<IndexedBreakend>> mSearchableBreakends;

    public SagaResource(final String fastaPath)
    {
        mFastaPath = fastaPath;
        mAssemblies = loadAssemblies();
        mVariantsById = mAssemblies.stream()
                .map(AssemblyMetadata::variant)
                .collect(Collectors.toMap(Variant::id, identity()));
        mSearchableBreakends = mAssemblies.stream()
                .flatMap(assembly ->
                        assembly.variant.breakends().map(breakend ->
                                new IndexedBreakend(breakend, assembly.variant.id())))
                .collect(Collectors.groupingBy(IndexedBreakend::chromosome));
    }

//    public String bwaIndexImagePath()
//    {
//        return mFastaPath + ".img";
//    }

//    public Map<String, AssemblyMetadata> assembliesByFastaLabel()
//    {
//        return mAssemblies.stream().collect(Collectors.toMap(AssemblyMetadata::fastaLabel, identity()));
//    }

    public Variant getVariantById(final String variantId)
    {
        Variant variant = mVariantsById.get(variantId);
        if(variant == null)
        {
            throw new IllegalArgumentException("No SAGA variant with ID:" + variantId);
        }
        else
        {
            return variant;
        }
    }

    public Map<String, Variant> variantsById()
    {
        return mVariantsById;
    }

    public Map<String, List<IndexedBreakend>> searchableBreakends()
    {
        return mSearchableBreakends;
    }

    public record Breakend(
            BasePosition position,
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
            return new Breakend(new BasePosition(chromosome, position), orientation);
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
        public BasePosition position()
        {
            return breakend.position();
        }

        public String chromosome()
        {
            return position().Chromosome;
        }
    }

    private List<AssemblyMetadata> loadAssemblies()
    {
        SV_LOGGER.debug("Loading SAGA resource FASTA");
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

        SV_LOGGER.debug("Loaded {} SAGA variant assemblies", assemblies.size());

        return assemblies;
    }

    private IndexedFastaSequenceFile loadFasta() throws IOException
    {
        return new IndexedFastaSequenceFile(new File(mFastaPath));
    }
}
